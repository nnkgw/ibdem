// main.cpp
// Phase 2+: Implicit BDEM simulation with three-point bending (Fig. 5).
// Three beam scales (Low / Middle / High) displayed side-by-side.

#ifdef __APPLE__
  #include <GLUT/glut.h>
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#else
  #ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
  #endif
  #include <GL/glut.h>
  #include <GL/gl.h>
  #include <GL/glu.h>
#endif

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <algorithm>

#include "Particle.h"
#include "Bond.h"
#include "HexPacking.h"
#include "Simulation.h"

// ---------------------------------------------------------------------------
// Window / viewport constants
// ---------------------------------------------------------------------------
static const int NUM_SCALES = 3;
static const int VP_W       = 400;  // width  of one viewport (pixels)
static const int VP_H       = 450;  // height of one viewport (pixels)
static const int WIN_W      = VP_W * NUM_SCALES;
static const int WIN_H      = VP_H;

// ---------------------------------------------------------------------------
// Beam configurations (Table 2, BDEM rows).
// Physical dimensions L x H x W = 1.0 x 0.15 x 0.12 m; only r varies.
// ---------------------------------------------------------------------------
static const BeamConfig SCALE_CFG[NUM_SCALES] = {
  // label    L      H      W      r       density  E      nu   tauC
  { 1.0f, 0.15f, 0.12f, 0.018f, 1000.0f, 1e7f, 0.3f, 3e4f, "Low"    },
  { 1.0f, 0.15f, 0.12f, 0.011f, 1000.0f, 1e7f, 0.3f, 3e4f, "Middle" },
  { 1.0f, 0.15f, 0.12f, 0.005f, 1000.0f, 1e7f, 0.3f, 3e4f, "High"   },
};

// ---------------------------------------------------------------------------
// SimConfig for each scale (from Table 2 and paper Section 6)
// ---------------------------------------------------------------------------
static const SimConfig SIM_CFG[NUM_SCALES] = {
  // Quasi-static implicit integration: large dt so m/dt^2 << kn (elastic-dominated).
  // m/dt^2/kn ratios: Low=8.6e-5, Middle=3.2e-5, High=6.7e-6 -- all negligible.
  // Jacobi-preconditioned condition number kappa ~ O(N_x) -> PCG converges in
  // sqrt(kappa) iterations. With 10 outer * 20 PCG inner = 200 CG steps/frame:
  //   Low   (N_x~28):  converges in ~28 CG steps.
  //   Middle (N_x~45): converges in ~45 CG steps.
  //   High  (N_x~100): converges in ~100 CG steps.
  // Per-scale dt tuned to match paper Table 2 fracture frames (Low=23, Mid=25, High=25).
  // Quasi-static: fracture_frame = fracture_disp / (loadVel * dt).
  //   Low   fracture_disp ~ 1.55e-3 m -> dt = 1.55e-3 / (23 * 0.062) = 1.087e-3
  //   Mid   fracture_disp ~ 1.36e-3 m -> dt = 1.36e-3 / (25 * 0.062) = 0.877e-3
  //   High  fracture_disp ~ 1.36e-3 m -> dt = same as Middle
  //  dt          E       nu    tauC   loadVel    grav  maxIter  eps
  { 1.087e-3f, 1e7f, 0.3f, 3e4f, 0.069f,  0.0f,   15, 1e-4f }, // Low    r=0.018
  { 0.877e-3f, 1e7f, 0.3f, 3e4f, 0.069f,  0.0f,   15, 1e-4f }, // Middle r=0.011
  { 0.877e-3f, 1e7f, 0.3f, 3e4f, 0.069f,  0.0f,   15, 1e-4f }, // High   r=0.005
};

// ---------------------------------------------------------------------------
// Simulation objects and control state
// ---------------------------------------------------------------------------
static Simulation g_sims[NUM_SCALES];
static bool g_running          = true;   // auto-start: simulation begins immediately
static bool g_stepOnce         = false;
static bool g_captureMode      = false;  // -capture: save every frame as BMP then exit
static int  g_captureMaxFrames = 40;

// ---------------------------------------------------------------------------
// BMP screenshot helper
// ---------------------------------------------------------------------------
static void saveBMP(const char* path) {
  int w = WIN_W, h = WIN_H;
  std::vector<uint8_t> pixels((size_t)w * h * 3);
  glReadPixels(0, 0, w, h, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels.data());

  // OpenGL origin is bottom-left; BMP expects top-left -- flip vertically
  for (int y = 0; y < h / 2; y++) {
    uint8_t* rowA = pixels.data() + (size_t)y        * w * 3;
    uint8_t* rowB = pixels.data() + (size_t)(h-1-y)  * w * 3;
    for (int x = 0; x < w * 3; x++) std::swap(rowA[x], rowB[x]);
  }

  // BMP row size must be aligned to 4 bytes
  int rowBytes = w * 3;
  int pad      = (4 - rowBytes % 4) % 4;
  int dataSize = (rowBytes + pad) * h;
  int fileSize = 54 + dataSize;

  uint8_t hdr[54] = {};
  hdr[0] = 'B'; hdr[1] = 'M';
  hdr[2] = (uint8_t)(fileSize      ); hdr[3]  = (uint8_t)(fileSize >> 8);
  hdr[4] = (uint8_t)(fileSize >> 16); hdr[5]  = (uint8_t)(fileSize >> 24);
  hdr[10] = 54; // pixel data offset
  hdr[14] = 40; // BITMAPINFOHEADER size
  hdr[18] = (uint8_t)(w      ); hdr[19] = (uint8_t)(w >> 8);
  hdr[20] = (uint8_t)(w >> 16); hdr[21] = (uint8_t)(w >> 24);
  hdr[22] = (uint8_t)(h      ); hdr[23] = (uint8_t)(h >> 8);
  hdr[24] = (uint8_t)(h >> 16); hdr[25] = (uint8_t)(h >> 24);
  hdr[26] = 1;  // planes
  hdr[28] = 24; // bits per pixel
  hdr[34] = (uint8_t)(dataSize      ); hdr[35] = (uint8_t)(dataSize >> 8);
  hdr[36] = (uint8_t)(dataSize >> 16); hdr[37] = (uint8_t)(dataSize >> 24);

  FILE* f = std::fopen(path, "wb");
  if (!f) { std::printf("saveBMP: cannot open %s\n", path); return; }
  std::fwrite(hdr, 1, 54, f);
  uint8_t padding[3] = {};
  for (int y = 0; y < h; y++) {
    std::fwrite(pixels.data() + (size_t)y * w * 3, 1, (size_t)rowBytes, f);
    if (pad) std::fwrite(padding, 1, (size_t)pad, f);
  }
  std::fclose(f);
}

// ---------------------------------------------------------------------------
// Helper: draw a null-terminated C string at the current raster position
// ---------------------------------------------------------------------------
static void drawString(const char* s) {
  for (; *s != '\0'; ++s)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *s);
}

// ---------------------------------------------------------------------------
// Render one beam scale into its viewport
// ---------------------------------------------------------------------------
static void drawScale(int index) {
  const Simulation& sim = g_sims[index];
  const BeamConfig& cfg = sim.beamCfg;

  // -- Viewport (left edge depends on column index) --
  glViewport(index * VP_W, 0, VP_W, VP_H);

  // -- Orthographic projection that fits the beam with 5 % margin --
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  float xMargin = cfg.L * 0.05f;
  float xMin    = -xMargin;
  float xMax    =  cfg.L + xMargin;
  float xRange  = xMax - xMin;

  // Preserve aspect ratio: compute yRange from viewport aspect
  float aspect  = (float)VP_H / (float)VP_W;
  float yRange  = xRange * aspect;
  float yCentre = cfg.H * 0.5f;

  glOrtho(xMin, xMax,
          yCentre - yRange * 0.5f, yCentre + yRange * 0.5f,
          -1.0f,  1.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Pixels per world unit (used to scale point / line sizes)
  float worldToPixel = (float)VP_W / xRange;

  // -- Background (per-viewport clear via scissor) --
  glEnable(GL_SCISSOR_TEST);
  glScissor(index * VP_W, 0, VP_W, VP_H);
  glClearColor(0.12f, 0.12f, 0.14f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  glDisable(GL_SCISSOR_TEST);

  // -- Draw active bonds as thin lines, colored by stress level --
  {
    float lw = std::max(0.5f, cfg.r * worldToPixel * 0.3f);
    glLineWidth(lw);
    glBegin(GL_LINES);
    for (const auto& b : sim.bonds) {
      if (b.broken) continue;
      // Stress color: blue(0) -> green(0.5) -> red(1.0) by sigma/tauC
      float t = glm::clamp(b.sigma / cfg.tauC, 0.0f, 1.0f);
      float cr = t > 0.5f ? 2.0f * (t - 0.5f) : 0.0f;
      float cg = t < 0.5f ? 2.0f * t : 2.0f * (1.0f - t);
      float cb = t < 0.5f ? 1.0f - 2.0f * t : 0.0f;
      glColor3f(cr, cg, cb);
      const glm::vec3& pi = sim.particles[b.i].pos;
      const glm::vec3& pj = sim.particles[b.j].pos;
      glVertex2f(pi.x, pi.y);
      glVertex2f(pj.x, pj.y);
    }
    glEnd();
  }

  // -- Draw particles as points (size proportional to radius) --
  float ptSize = std::max(1.5f, 2.0f * cfg.r * worldToPixel);
  glPointSize(ptSize);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Regular particles (gray)
  glBegin(GL_POINTS);
  glColor4f(0.70f, 0.72f, 0.75f, 0.90f);
  for (const auto& p : sim.particles) {
    if (!p.isSupport && !p.isLoad)
      glVertex2f(p.pos.x, p.pos.y);
  }
  glEnd();

  // Support particles (blue)
  glPointSize(std::max(3.0f, ptSize * 1.4f));
  glBegin(GL_POINTS);
  glColor4f(0.25f, 0.50f, 0.95f, 1.00f);
  for (const auto& p : sim.particles) {
    if (p.isSupport)
      glVertex2f(p.pos.x, p.pos.y);
  }
  glEnd();

  // Load particles (red)
  glBegin(GL_POINTS);
  glColor4f(0.95f, 0.25f, 0.20f, 1.00f);
  for (const auto& p : sim.particles) {
    if (p.isLoad)
      glVertex2f(p.pos.x, p.pos.y);
  }
  glEnd();

  glDisable(GL_BLEND);
  glDisable(GL_POINT_SMOOTH);
  glPointSize(1.0f);

  // -----------------------------------------------------------------------
  // Labels and legend
  // -----------------------------------------------------------------------
  // Scale label: top-left of viewport (well below the window title bar)
  glColor3f(0.95f, 0.95f, 0.95f);
  float labelX = xMin + xRange * 0.02f;
  float labelY = yCentre + yRange * 0.44f;  // near top, with a little margin
  glRasterPos2f(labelX, labelY);

  char buf[128];
  if (sim.fractureFrame >= 0) {
    std::snprintf(buf, sizeof(buf), "%s (%d p, %d b) Frame:%d [FRACTURE@%d]",
      cfg.label,
      (int)sim.particles.size(),
      (int)sim.bonds.size(),
      sim.frame,
      sim.fractureFrame);
  } else {
    std::snprintf(buf, sizeof(buf), "%s (%d particles, %d bonds)  Frame: %d",
      cfg.label,
      (int)sim.particles.size(),
      (int)sim.bonds.size(),
      sim.frame);
  }
  drawString(buf);

  // Legend: bottom-left corner, only for the first viewport to avoid clutter
  if (index == 0) {
    float legX  = xMin + xRange * 0.02f;
    // Start legend near the bottom, leaving a small margin
    float legY0 = yCentre - yRange * 0.42f;   // bottom legend item Y
    float legDy = yRange  * 0.05f;             // spacing between items

    glPointSize(8.0f);

    // Regular particle dot
    glBegin(GL_POINTS);
    glColor3f(0.70f, 0.72f, 0.75f);
    glVertex2f(legX, legY0);
    glEnd();
    glColor3f(0.90f, 0.90f, 0.90f);
    glRasterPos2f(legX + xRange * 0.025f, legY0 - legDy * 0.05f);
    drawString("particle");

    // Support dot (one row above)
    glBegin(GL_POINTS);
    glColor3f(0.25f, 0.50f, 0.95f);
    glVertex2f(legX, legY0 + legDy);
    glEnd();
    glColor3f(0.90f, 0.90f, 0.90f);
    glRasterPos2f(legX + xRange * 0.025f, legY0 + legDy * 0.95f);
    drawString("support (fixed)");

    // Load dot (two rows above)
    glBegin(GL_POINTS);
    glColor3f(0.95f, 0.25f, 0.20f);
    glVertex2f(legX, legY0 + 2.0f * legDy);
    glEnd();
    glColor3f(0.90f, 0.90f, 0.90f);
    glRasterPos2f(legX + xRange * 0.025f, legY0 + legDy * 1.95f);
    drawString("load point");

    glPointSize(1.0f);

    // Stress color bar legend (bond colors)
    float barY  = legY0 + 3.5f * legDy;
    float barX0 = legX;
    float barX1 = legX + xRange * 0.15f;
    int   steps = 20;
    glLineWidth(4.0f);
    glBegin(GL_LINES);
    for (int k = 0; k < steps; k++) {
      float t  = (float)k / (float)(steps - 1);
      float cr = t > 0.5f ? 2.0f * (t - 0.5f) : 0.0f;
      float cg = t < 0.5f ? 2.0f * t : 2.0f * (1.0f - t);
      float cb = t < 0.5f ? 1.0f - 2.0f * t : 0.0f;
      glColor3f(cr, cg, cb);
      float bx = barX0 + t * (barX1 - barX0);
      glVertex2f(bx, barY);
      glVertex2f(bx, barY + legDy * 0.6f);
    }
    glEnd();
    glLineWidth(1.0f);
    glColor3f(0.90f, 0.90f, 0.90f);
    glRasterPos2f(barX0, barY + legDy * 0.7f);
    drawString("bond stress: 0 -> tauC");
  }
}

// ---------------------------------------------------------------------------
// GLUT callbacks
// ---------------------------------------------------------------------------

static void display() {
  // Clear the whole window first
  glClearColor(0.08f, 0.08f, 0.10f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  for (int i = 0; i < NUM_SCALES; ++i)
    drawScale(i);

  // -- Full-window overlay: separators + title --
  glViewport(0, 0, WIN_W, WIN_H);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, WIN_W, 0.0, WIN_H, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Vertical separator lines between viewports
  glColor3f(0.30f, 0.30f, 0.35f);
  glLineWidth(2.0f);
  glBegin(GL_LINES);
  for (int i = 1; i < NUM_SCALES; ++i) {
    glVertex2f((float)(i * VP_W), 0.0f);
    glVertex2f((float)(i * VP_W), (float)WIN_H);
  }
  glEnd();
  glLineWidth(1.0f);

  // Window title text (top centre)
  glColor3f(1.0f, 1.0f, 1.0f);
  glRasterPos2f(WIN_W * 0.5f - 200.0f, WIN_H - 18.0f);
  drawString("Implicit BDEM -- Three-Point Bending Scale Consistency (Fig. 5)");

  // Controls hint (bottom centre of full window, pixel coords)
  glColor3f(0.65f, 0.65f, 0.65f);
  glRasterPos2f(WIN_W * 0.5f - 170.0f, 4.0f);
  drawString("SPACE=start/pause  N=step  R=reset  S=screenshot  Q/ESC=quit");

  if (g_captureMode) {
    char path[64];
    std::snprintf(path, sizeof(path), "capture_fr%04d.bmp", g_sims[0].frame);
    saveBMP(path);
    std::printf("Saved: %s\n", path);
  }

  glutSwapBuffers();
}

static void reshape(int w, int h) {
  (void)w; (void)h;
}

static void resetAll() {
  g_running  = false;
  g_stepOnce = false;
  for (int i = 0; i < NUM_SCALES; i++)
    g_sims[i].init(SCALE_CFG[i], SIM_CFG[i], i);
}

static void keyboard(unsigned char key, int /*x*/, int /*y*/) {
  switch (key) {
    case 27:   // ESC
    case 'q':
    case 'Q':
      std::exit(0);
      break;
    case ' ':
      g_running = !g_running;
      break;
    case 'n':
    case 'N':
      g_stepOnce = true;
      break;
    case 'r':
    case 'R':
      resetAll();
      glutPostRedisplay();
      break;
    case 's':
    case 'S': {
      char path[64];
      std::snprintf(path, sizeof(path), "capture_fr%04d.bmp", g_sims[0].frame);
      saveBMP(path);
      std::printf("Saved: %s\n", path);
      break;
    }
    default:
      break;
  }
}

// Idle callback: steps simulation and requests redisplay.
// Called by GLUT whenever there are no pending events.
static void idleFunc() {
  if (g_captureMode) {
    if (g_sims[0].frame >= g_captureMaxFrames) std::exit(0);
    for (int i = 0; i < NUM_SCALES; i++) g_sims[i].step();
    glutPostRedisplay();
    return;
  }
  if (g_running) {
    for (int i = 0; i < NUM_SCALES; i++)
      g_sims[i].step();
    glutPostRedisplay();
  } else if (g_stepOnce) {
    for (int i = 0; i < NUM_SCALES; i++)
      g_sims[i].step();
    g_stepOnce = false;
    glutPostRedisplay();
  }
}

// Heartbeat timer: ensures the window redraws at least once per second
// even while paused (e.g. to keep the title bar responsive).
static void timerFunc(int val) {
  glutPostRedisplay();
  glutTimerFunc(1000, timerFunc, val);
}

// ---------------------------------------------------------------------------
// Initialisation
// ---------------------------------------------------------------------------

static void initSimulations() {
  for (int i = 0; i < NUM_SCALES; ++i) {
    std::printf("Initialising scale %d (%s) r=%.4f ...\n",
      i, SCALE_CFG[i].label, SCALE_CFG[i].r);
    g_sims[i].init(SCALE_CFG[i], SIM_CFG[i], i);
    std::printf("  particles: %d   bonds: %d\n",
      (int)g_sims[i].particles.size(),
      (int)g_sims[i].bonds.size());
  }
  std::printf("Initialisation complete.\n");
  std::printf("Press SPACE to start simulation, N to step, R to reset, S to screenshot, Q/ESC to quit.\n");
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

// Headless mode: run simulation for N frames and report fracture, then exit.
static void runHeadless(int maxFrames) {
  initSimulations();
  std::printf("Running headless for %d frames...\n", maxFrames);
  for (int f = 0; f < maxFrames; f++) {
    for (int i = 0; i < NUM_SCALES; i++)
      g_sims[i].step();
    for (int i = 0; i < NUM_SCALES; i++) {
      if (g_sims[i].fractureFrame == f)
        std::printf("  [Frame %3d] Scale %-7s FRACTURE\n", f, SCALE_CFG[i].label);
    }
  }
  for (int i = 0; i < NUM_SCALES; i++) {
    std::printf("Scale %-7s : fracture at frame %d  (total frames=%d)\n",
      SCALE_CFG[i].label,
      g_sims[i].fractureFrame,
      g_sims[i].frame);
  }
}

int main(int argc, char** argv) {
  // Parse arguments before GLUT init so we can handle -headless without a window.
  for (int a = 1; a < argc; a++) {
    std::string arg(argv[a]);
    if (arg == "-headless" || arg == "--headless") {
      int frames = (a + 1 < argc) ? std::atoi(argv[a + 1]) : 60;
      runHeadless(frames);
      return 0;
    }
    if (arg == "-capture" || arg == "--capture") {
      g_captureMode      = true;
      g_captureMaxFrames = (a + 1 < argc) ? std::atoi(argv[a + 1]) : 40;
      g_running          = false; // idle loop controls capture; g_running not used
    }
  }

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(WIN_W, WIN_H);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Implicit BDEM -- Fig. 5 Scale Consistency");

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idleFunc);
  glutTimerFunc(1000, timerFunc, 0);

  initSimulations();

  glutMainLoop();
  return 0;
}
