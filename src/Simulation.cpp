#include "Simulation.h"
#include "BondForce.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>

// ---------------------------------------------------------------------------
// Helper: geodesic step on S3 from quaternion q along tangent dq by alpha
// ---------------------------------------------------------------------------
static glm::quat manifoldStep(const glm::quat& q, const glm::vec4& tangent, float alpha) {
  glm::vec4 dq = alpha * tangent;
  float norm_dq = glm::length(dq);
  if (norm_dq < 1e-14f) return q;
  glm::vec4 q_vec(q.x, q.y, q.z, q.w);
  glm::vec4 q_new = q_vec * std::cos(norm_dq) + (dq / norm_dq) * std::sin(norm_dq);
  // Normalize to stay on unit sphere
  float len = glm::length(q_new);
  if (len < 1e-12f) return q;
  q_new /= len;
  return glm::quat(q_new.w, q_new.x, q_new.y, q_new.z);
}

// ---------------------------------------------------------------------------
// Simulation::init
// ---------------------------------------------------------------------------
void Simulation::init(const BeamConfig& bcfg, const SimConfig& scfg, int idx) {
  beamCfg = bcfg;
  simCfg  = scfg;
  frame        = 0;
  fractureFrame = -1;
  scaleIdx     = idx;

  particles = generateBeam(bcfg);
  tagBoundaryParticles(particles, bcfg);
  bonds = generateBonds(particles);

  int N = (int)particles.size();
  p_hat.resize(N);
  q_hat.resize(N);
  pcg_r.resize(N);
  pcg_d.resize(N);
  pcg_z.resize(N);
  pcg_Hd.resize(N);
}

// ---------------------------------------------------------------------------
// Simulation::totalEnergy
// ---------------------------------------------------------------------------
float Simulation::totalEnergy() const {
  float E = 0.0f;
  for (const auto& b : bonds) {
    if (b.broken) continue;
    E += bondEnergy(particles[b.i], particles[b.j], b, simCfg.E, simCfg.nu);
  }
  return E;
}

// ---------------------------------------------------------------------------
// Simulation::solvePCG
// Preconditioned Conjugate Gradient for position DOFs (rotations held fixed).
// Hessian-vector product uses analytical stretch stiffness (dominant term).
// Converges in O(sqrt(kappa)) iterations where kappa = (1+rho)/(1-rho) ~ 56
// for rho=0.965, giving ~8 iterations vs 1000+ for Jacobi.
// ---------------------------------------------------------------------------
float Simulation::solvePCG(const std::vector<bool>& is_load,
                            const std::vector<glm::vec3>& load_target,
                            int maxIter, float epsilon)
{
  int N = (int)particles.size();
  float dt = simCfg.dt;
  float pi_f = glm::pi<float>();

  // ---- Diagonal preconditioner: stretch + shear transverse stiffness ----
  std::vector<float> h_p(N, 0.0f);
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed || is_load[i]) continue;
    h_p[i] = particles[i].mass / (dt * dt);
  }
  for (const auto& b : bonds) {
    if (b.broken) continue;
    float r0  = (particles[b.i].radius + particles[b.j].radius) * 0.5f;
    float S   = pi_f * r0 * r0;
    float kn  = simCfg.E * S / b.l0;
    float G   = simCfg.E / (2.0f * (1.0f + simCfg.nu));
    float kt  = G * S / b.l0;
    // Stretch diagonal contribution
    if (!particles[b.i].fixed && !is_load[b.i]) h_p[b.i] += kn;
    if (!particles[b.j].fixed && !is_load[b.j]) h_p[b.j] += kn;
    // Shear transverse stiffness: kt * l0^2 / dist^2 (two transverse DOFs)
    glm::vec3 dv = particles[b.j].pos - particles[b.i].pos;
    float dist2 = glm::dot(dv, dv);
    if (dist2 > 1e-20f) {
      float ksh = kt * b.l0 * b.l0 / dist2;
      if (!particles[b.i].fixed && !is_load[b.i]) h_p[b.i] += 2.0f * ksh;
      if (!particles[b.j].fixed && !is_load[b.j]) h_p[b.j] += 2.0f * ksh;
    }
  }

  // ---- Lambda: compute position gradient g = (m/dt^2)(p-p_hat) + bondGrad_p ----
  auto computeGrad = [&](std::vector<glm::vec3>& g) {
    std::fill(g.begin(), g.end(), glm::vec3(0.0f));
    for (int i = 0; i < N; i++) {
      if (particles[i].fixed || is_load[i]) continue;
      float m = particles[i].mass;
      g[i] = (m / (dt * dt)) * (particles[i].pos - p_hat[i]);
    }
    for (const auto& b : bonds) {
      if (b.broken) continue;
      BondGradient bg = bondGradient(particles[b.i], particles[b.j],
                                     b, simCfg.E, simCfg.nu);
      if (!particles[b.i].fixed && !is_load[b.i]) g[b.i] += bg.grad_pi;
      if (!particles[b.j].fixed && !is_load[b.j]) g[b.j] += bg.grad_pj;
    }
  };

  // ---- Lambda: Hessian-vector product Hd = (m/dt^2)d + (K_stretch + K_shear) * d ----
  // K_stretch[i,j] = -kn * dhat x dhat  (axial)
  // K_shear[i,j]   = -kt*l0^2/dist^2 * (I - dhat x dhat)  (transverse)
  auto hessVec = [&](const std::vector<glm::vec3>& d_vec,
                     std::vector<glm::vec3>&       Hd_vec) {
    std::fill(Hd_vec.begin(), Hd_vec.end(), glm::vec3(0.0f));
    for (int i = 0; i < N; i++) {
      if (particles[i].fixed || is_load[i]) continue;
      float m = particles[i].mass;
      Hd_vec[i] = (m / (dt * dt)) * d_vec[i];
    }
    for (const auto& b : bonds) {
      if (b.broken) continue;
      glm::vec3 dd = particles[b.j].pos - particles[b.i].pos;
      float dist = glm::length(dd);
      if (dist < 1e-12f) continue;
      glm::vec3 dhat = dd / dist;
      float r0  = (particles[b.i].radius + particles[b.j].radius) * 0.5f;
      float S   = pi_f * r0 * r0;
      float kn  = simCfg.E * S / b.l0;
      float G   = simCfg.E / (2.0f * (1.0f + simCfg.nu));
      float kt  = G * S / b.l0;
      float ksh = kt * b.l0 * b.l0 / (dist * dist);  // shear transverse stiffness
      glm::vec3 dvi = (particles[b.i].fixed || is_load[b.i]) ? glm::vec3(0.0f) : d_vec[b.i];
      glm::vec3 dvj = (particles[b.j].fixed || is_load[b.j]) ? glm::vec3(0.0f) : d_vec[b.j];
      glm::vec3 dv_rel = dvj - dvi;
      // Stretch: kn * dhat (x) dhat
      float proj_n = kn * glm::dot(dhat, dv_rel);
      // Shear: ksh * (I - dhat x dhat)  (transverse)
      glm::vec3 proj_t = ksh * (dv_rel - glm::dot(dhat, dv_rel) * dhat);
      if (!particles[b.i].fixed && !is_load[b.i]) {
        Hd_vec[b.i] -= proj_n * dhat + proj_t;
      }
      if (!particles[b.j].fixed && !is_load[b.j]) {
        Hd_vec[b.j] += proj_n * dhat + proj_t;
      }
    }
  };

  // ---- PCG iterations ----
  // r = -gradient (residual of the linear system)
  computeGrad(pcg_r);
  for (int i = 0; i < N; i++) pcg_r[i] = -pcg_r[i];

  // z = P^{-1} r
  for (int i = 0; i < N; i++) {
    if (h_p[i] > 1e-30f) pcg_z[i] = pcg_r[i] / h_p[i];
    else                  pcg_z[i] = glm::vec3(0.0f);
  }

  // d = z
  pcg_d = pcg_z;

  float rz = 0.0f;
  for (int i = 0; i < N; i++) rz += glm::dot(pcg_r[i], pcg_z[i]);

  float rz0 = rz;  // initial preconditioned residual norm squared (for relative check)
  int   pcg_iters_used = 0;

  for (int iter = 0; iter < maxIter; ++iter) {
    // Relative convergence: stop when ||r||_P / ||r0||_P < epsilon
    if (rz0 > 0.0f && rz < epsilon * epsilon * rz0) break;
    if (rz < 1e-30f) break;

    hessVec(pcg_d, pcg_Hd);

    float dHd = 0.0f;
    for (int i = 0; i < N; i++) dHd += glm::dot(pcg_d[i], pcg_Hd[i]);
    if (dHd < 1e-30f) break;

    float alpha_pcg = rz / dHd;

    // Update positions and residual
    for (int i = 0; i < N; i++) {
      if (particles[i].fixed || is_load[i]) continue;
      particles[i].pos += alpha_pcg * pcg_d[i];
      pcg_r[i]         -= alpha_pcg * pcg_Hd[i];
    }

    // Update preconditioned residual
    for (int i = 0; i < N; i++) {
      if (h_p[i] > 1e-30f) pcg_z[i] = pcg_r[i] / h_p[i];
      else                  pcg_z[i] = glm::vec3(0.0f);
    }

    float rz_new = 0.0f;
    for (int i = 0; i < N; i++) rz_new += glm::dot(pcg_r[i], pcg_z[i]);

    float beta = rz_new / rz;
    for (int i = 0; i < N; i++) pcg_d[i] = pcg_z[i] + beta * pcg_d[i];

    rz = rz_new;
    pcg_iters_used = iter + 1;
  }



  float gnorm_final = 0.0f;
  for (int i = 0; i < N; i++) gnorm_final += glm::dot(pcg_r[i], pcg_r[i]);
  return std::sqrt(gnorm_final);
}

// ---------------------------------------------------------------------------
// Simulation::gradientStep
// Compute full gradient and apply one preconditioned manifold step.
// Diagonal preconditioner: h_p[i] = m/dt^2 + sum_bonds(kn)
//                          h_q[i] = Iq/dt^2 + sum_bonds(EI/l0)
// This gives Jacobi-like convergence in O(1) iterations.
// alpha scales the preconditioned step (1.0 = full Jacobi step).
// ---------------------------------------------------------------------------
void Simulation::gradientStep(float alpha, bool update_pos) {
  int N = (int)particles.size();
  float dt = simCfg.dt;
  float pi_f = glm::pi<float>();

  // ---- Build diagonal preconditioner ----
  // fixed particles (support points): position pinned, rotation FREE (pin/roller support).
  // This simulates simply-supported conditions: rotation at supports can adapt to the beam.
  std::vector<float> h_p(N, 0.0f);
  std::vector<float> h_q(N, 0.0f);
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) continue;  // fixed particles get h_p from bonds only
    float m  = particles[i].mass;
    float r  = particles[i].radius;
    float Iq = (8.0f / 5.0f) * m * r * r;
    h_p[i] = m  / (dt * dt);
    h_q[i] = Iq / (dt * dt);  // inertia; fixed particles start at 0
  }
  for (const auto& b : bonds) {
    if (b.broken) continue;
    float r0  = (particles[b.i].radius + particles[b.j].radius) * 0.5f;
    float S   = pi_f * r0 * r0;
    float kn  = simCfg.E * S / b.l0;
    float G   = simCfg.E / (2.0f * (1.0f + simCfg.nu));
    float kt  = G * S / b.l0;
    float I   = pi_f * r0 * r0 * r0 * r0 / 4.0f;
    float EIl = simCfg.E * I / b.l0;
    float kt_ang = kt * b.l0 * b.l0;  // [Nm/rad]
    // Position preconditioner: only stretch (shear has no position gradient)
    if (!particles[b.i].fixed) h_p[b.i] += kn;
    if (!particles[b.j].fixed) h_p[b.j] += kn;
    // Rotation preconditioner: ALL particles (including fixed) can rotate
    h_q[b.i] += EIl + kt_ang;
    h_q[b.j] += EIl + kt_ang;
  }

  // ---- Accumulate gradients ----
  std::vector<glm::vec3> grad_p(N, glm::vec3(0.0f));
  std::vector<glm::vec4> grad_q(N, glm::vec4(0.0f));

  // Inertia terms: only for non-fixed particles (fixed have no positional inertia;
  // rotation inertia is also zeroed for fixed so supports are quasi-static in rotation)
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) continue;
    float m  = particles[i].mass;
    float r  = particles[i].radius;
    float Iq = (8.0f / 5.0f) * m * r * r;

    grad_p[i] += (m / (dt * dt)) * (particles[i].pos - p_hat[i]);

    glm::vec4 q_curr(particles[i].rot.x, particles[i].rot.y,
                     particles[i].rot.z, particles[i].rot.w);
    glm::vec4 q_hat_v(q_hat[i].x, q_hat[i].y, q_hat[i].z, q_hat[i].w);
    grad_q[i] += (Iq / (dt * dt)) * (q_curr - q_hat_v);
  }

  // Bond energy gradients (ALL particles, including fixed, contribute rotation gradient)
  for (const auto& b : bonds) {
    if (b.broken) continue;
    BondGradient bg = bondGradient(particles[b.i], particles[b.j],
                                   b, simCfg.E, simCfg.nu);
    if (!particles[b.i].fixed) grad_p[b.i] += bg.grad_pi;
    if (!particles[b.j].fixed) grad_p[b.j] += bg.grad_pj;
    grad_q[b.i] += bg.grad_qi;  // rotation gradient for ALL particles
    grad_q[b.j] += bg.grad_qj;
  }

  // Project q gradients onto tangent space of S3 (ALL particles including fixed)
  for (int i = 0; i < N; i++) {
    glm::vec4 q(particles[i].rot.x, particles[i].rot.y,
                particles[i].rot.z, particles[i].rot.w);
    grad_q[i] -= glm::dot(q, grad_q[i]) * q;  // (I - q*q^T) * grad_q
  }

  // ---- Apply preconditioned updates ----
  for (int i = 0; i < N; i++) {
    // Position update: only non-fixed particles (skip if update_pos=false)
    if (update_pos && !particles[i].fixed && h_p[i] > 1e-30f) {
      particles[i].pos -= (alpha / h_p[i]) * grad_p[i];
    }

    // Rotation update: ALL particles (including fixed/support) — simply-supported BC
    if (h_q[i] > 1e-30f) {
      glm::vec4 tang = -(alpha / h_q[i]) * grad_q[i];
      particles[i].rot = manifoldStep(particles[i].rot, tang, 1.0f);
    }
  }
}

// ---------------------------------------------------------------------------
// Simulation::checkFracture
// ---------------------------------------------------------------------------
int Simulation::checkFracture() {
  int count = 0;
  float maxSigma = 0.0f, maxTau = 0.0f;
  int maxTauBond = -1;
  for (int bi = 0; bi < (int)bonds.size(); bi++) {
    auto& b = bonds[bi];
    if (b.broken) continue;
    // Skip bonds directly adjacent to load/support particles (boundary artifacts)
    if (particles[b.i].isLoad || particles[b.j].isLoad) continue;
    if (particles[b.i].fixed || particles[b.j].fixed)   continue;
    float sigma = 0.0f, tau = 0.0f;
    bondStress(particles[b.i], particles[b.j], b, simCfg.E, simCfg.nu, sigma, tau);
    b.sigma = sigma;  // update for ALL interior bonds — used for rendering color
    if (sigma > maxSigma) { maxSigma = sigma; }
    if (tau   > maxTau)   { maxTau   = tau;   maxTauBond = bi; }

    // Before the initial fracture: restrict to the tension zone (bottom half of the beam).
    // This prevents premature shear fracture near the load-application point, which
    // creates a stress concentration at the top.  Once the first fracture has occurred
    // (fractureFrame >= 0), all bonds are eligible so the crack can propagate upward.
    if (fractureFrame < 0) {
      float y_mid = (particles[b.i].pos.y + particles[b.j].pos.y) * 0.5f;
      if (y_mid > beamCfg.H * 0.5f) continue;
    }
    if (sigma > simCfg.tauC || tau > simCfg.tauC) {
      b.broken = true;
      ++count;
    }
  }
  (void)maxTauBond;
  return count;
}

// ---------------------------------------------------------------------------
// Simulation::step
// One full implicit time step using manifold gradient descent.
// ---------------------------------------------------------------------------
void Simulation::step() {
  int   N  = (int)particles.size();
  float dt = simCfg.dt;

  // ---- 1. Predict: quasi-static (no velocity warm-start) ----
  // p_hat = current position.  Each frame finds the static equilibrium under
  // the new load increment.  With dt=1e-3 (large), m/dt^2 << kn so the
  // inertia term is a weak regularizer and x* ≈ argmin V(x) (quasi-static).
  // No velocity warm-start avoids dynamic oscillations across frames.
  for (int i = 0; i < N; i++) {
    p_hat[i] = particles[i].pos;
    q_hat[i] = particles[i].rot;
  }

  // ---- 2. Apply load displacement: move load particles downward each frame ----
  // Load particles receive a prescribed displacement; they are NOT fixed so
  // the optimizer can redistribute force, but we reset their position after
  // each solve to enforce the kinematic constraint (velocity-driven loading).
  // We store their target position here and clamp after each gradient step.
  std::vector<glm::vec3> load_target(N, glm::vec3(0.0f));
  std::vector<bool>      is_load(N, false);
  for (int i = 0; i < N; i++) {
    if (particles[i].isLoad && !particles[i].fixed) {
      is_load[i]     = true;
      load_target[i] = particles[i].pos;
      load_target[i].y -= simCfg.loadVel * dt;  // move downward
      // Initialize solver at target so inertia doesn't fight it
      p_hat[i]            = load_target[i];
      particles[i].pos    = load_target[i];
    }
  }

  // ---- 3. Set initial guess x0 = x_hat ----
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) continue;
    particles[i].pos = p_hat[i];
    particles[i].rot = q_hat[i];
  }

  // ---- 4. Solve positions + rotations with outer alternating loop ----
  // V_shear = (1/2)*kt*l0^2*angle(d_actual, qc⊙d0)^2 couples positions and orientations.
  // Block coordinate descent: alternate PCG for positions (stretch+shear) and
  // gradient descent for rotations (shear+bend-twist) until convergence.
  // Outer loop: 10 alternating iterations. Each iteration:
  //   (a) PCG for positions (with current orientations providing shear forces)
  //   (b) Rotation gradient descent (with current positions providing shear torques)
  int outerIters = simCfg.maxIter;  // typically 10
  int pcgPerOuter = 20;
  int rotPerOuter = 20;

  for (int outer = 0; outer < outerIters; ++outer) {
    // (a) PCG for positions — fixed orientations, solve positions
    solvePCG(is_load, load_target, pcgPerOuter, simCfg.epsilon);
    // Enforce position BCs after each PCG call
    for (int i = 0; i < N; i++) {
      if (particles[i].fixed) particles[i].pos = p_hat[i];
      if (is_load[i])         particles[i].pos = load_target[i];
    }

    // (b) Gradient descent for rotations — fixed positions, update rotations
    for (int r = 0; r < rotPerOuter; ++r) {
      gradientStep(1.0f, false);  // false = skip position update
    }
  }

  // ---- 5. No velocity update (quasi-static: no inertia accumulation) ----
  // Velocities stay at zero; p_hat = p_n each step so the solver always finds
  // the static equilibrium under the current load, not a dynamic trajectory.

  // ---- 6. Check fracture ----
  int broken = checkFracture();
  if (broken > 0 && fractureFrame < 0) {
    fractureFrame = frame;
  }

  ++frame;
}
