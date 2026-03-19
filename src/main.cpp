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
  // dt = 0.02975 * r  (from rho_Jacobi = 0.965 => m/dt^2 = 0.0363 * 6*(kn+kt))
  // loadVel = delta_crit / (25 * dt)  where delta_crit = 3.333mm (beam theory)
  // PCG converges in ~8 iterations; maxIter=20 is generous.
  //  dt          E       nu    tauC   loadVel  grav  maxIter  eps
  { 5.35e-4f, 1e7f, 0.3f, 3e4f, 0.2491f, 0.0f, 200, 1e-4f }, // Low    r=0.018
  { 3.27e-4f, 1e7f, 0.3f, 3e4f, 0.4075f, 0.0f, 200, 1e-4f }, // Middle r=0.011
  { 1.49e-4f, 1e7f, 0.3f, 3e4f, 0.8946f, 0.0f, 200, 1e-4f }, // High   r=0.005
};

// ---------------------------------------------------------------------------
// Simulation objects and control state
// ---------------------------------------------------------------------------
static Simulation g_sims[NUM_SCALES];
static bool       g_running   = false;
static bool       g_stepOnce  = false;

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

  // -- Draw active bonds as thin lines (XY projection) --
  {
    float lw = std::max(0.5f, cfg.r * worldToPixel * 0.3f);
    glLineWidth(lw);
    glColor3f(0.45f, 0.50f, 0.55f);
    glBegin(GL_LINES);
    for (const auto& b : sim.bonds) {
      if (b.broken) continue;
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
  glRasterPos2f(WIN_W * 0.5f - 140.0f, 4.0f);
  drawString("SPACE=start/pause  N=step  R=reset  Q/ESC=quit");

  glutSwapBuffers();
}

static void reshape(int w, int h) {
  (void)w; (void)h;
}

static void resetAll() {
  for (int i = 0; i < NUM_SCALES; i++) {
    g_sims[i].init(SCALE_CFG[i], SIM_CFG[i], i);
  }
  g_running  = false;
  g_stepOnce = false;
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
    default:
      break;
  }
}

// Timer fires ~60fps; advances simulation when running
static void timerFunc(int val) {
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
  glutTimerFunc(16, timerFunc, val);
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
  std::printf("Press SPACE to start simulation, N to step, R to reset, Q/ESC to quit.\n");
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
  // -headless [N]: run without window, print fracture frames, exit
  for (int a = 1; a < argc; a++) {
    std::string arg(argv[a]);
    if (arg == "-headless" || arg == "--headless") {
      int frames = (a + 1 < argc) ? std::atoi(argv[a + 1]) : 60;
      runHeadless(frames);
      return 0;
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
  glutTimerFunc(16, timerFunc, 0);

  initSimulations();

  glutMainLoop();
  return 0;
}
