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

  // ---- Diagonal preconditioner (same as Jacobi) ----
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
    if (!particles[b.i].fixed && !is_load[b.i]) h_p[b.i] += kn + kt;
    if (!particles[b.j].fixed && !is_load[b.j]) h_p[b.j] += kn + kt;
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

  // ---- Lambda: Hessian-vector product Hd = (m/dt^2)d + K_stretch * d ----
  // Uses analytical stretch stiffness: K_stretch[i,j] = -kn * dhat (x) dhat
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
      // Stretch: kn * dhat (x) dhat
      glm::vec3 dvi = (particles[b.i].fixed || is_load[b.i]) ? glm::vec3(0.0f) : d_vec[b.i];
      glm::vec3 dvj = (particles[b.j].fixed || is_load[b.j]) ? glm::vec3(0.0f) : d_vec[b.j];
      float proj_n = kn * glm::dot(dhat, dvj - dvi);
      // Shear (transverse): kt * (I - dhat(x)dhat) treated as kt*I - kt*dhat(x)dhat
      glm::vec3 proj_t_i = kt * (dvj - dvi) - kt * glm::dot(dhat, dvj - dvi) * dhat;
      if (!particles[b.i].fixed && !is_load[b.i]) {
        Hd_vec[b.i] -= proj_n * dhat + proj_t_i;
      }
      if (!particles[b.j].fixed && !is_load[b.j]) {
        Hd_vec[b.j] += proj_n * dhat + proj_t_i;
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
void Simulation::gradientStep(float alpha) {
  int N = (int)particles.size();
  float dt = simCfg.dt;
  float pi_f = glm::pi<float>();

  // ---- Build diagonal preconditioner ----
  std::vector<float> h_p(N, 0.0f);
  std::vector<float> h_q(N, 0.0f);
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) continue;
    float m  = particles[i].mass;
    float r  = particles[i].radius;
    float Iq = (8.0f / 5.0f) * m * r * r;
    h_p[i] = m  / (dt * dt);
    h_q[i] = Iq / (dt * dt);
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
    // Shear angular stiffness: V_shear = 0.5*kt*l0^2*theta^2 -> d^2V/dq^2 ~ kt*l0^2
    float kt_ang = kt * b.l0 * b.l0;  // [Nm/rad]
    if (!particles[b.i].fixed) { h_p[b.i] += kn + kt; h_q[b.i] += EIl + kt_ang; }
    if (!particles[b.j].fixed) { h_p[b.j] += kn + kt; h_q[b.j] += EIl + kt_ang; }
  }

  // ---- Accumulate gradients ----
  std::vector<glm::vec3> grad_p(N, glm::vec3(0.0f));
  std::vector<glm::vec4> grad_q(N, glm::vec4(0.0f));

  // Inertia terms: (m/dt^2)*(p - p_hat) and (Iq/dt^2)*(q - q_hat)
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

  // Bond energy gradients
  for (const auto& b : bonds) {
    if (b.broken) continue;
    BondGradient bg = bondGradient(particles[b.i], particles[b.j],
                                   b, simCfg.E, simCfg.nu);
    grad_p[b.i] += bg.grad_pi;
    grad_p[b.j] += bg.grad_pj;
    grad_q[b.i] += bg.grad_qi;
    grad_q[b.j] += bg.grad_qj;
  }

  // Project q gradients onto tangent space of S3 at each quaternion
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) continue;
    glm::vec4 q(particles[i].rot.x, particles[i].rot.y,
                particles[i].rot.z, particles[i].rot.w);
    grad_q[i] -= glm::dot(q, grad_q[i]) * q;  // (I - q*q^T) * grad_q
  }

  // ---- Apply preconditioned updates ----
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) continue;

    // Position: p -= alpha * grad_p / h_p[i]
    particles[i].pos -= (alpha / h_p[i]) * grad_p[i];

    // Orientation: geodesic step on S3 with preconditioned tangent
    glm::vec4 tang = -(alpha / h_q[i]) * grad_q[i];
    particles[i].rot = manifoldStep(particles[i].rot, tang, 1.0f);
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
    // Skip bonds in the top half of the beam (compression zone; fracture initiates at bottom)
    if (particles[b.i].pos.y > beamCfg.H * 0.5f ||
        particles[b.j].pos.y > beamCfg.H * 0.5f) continue;
    float sigma = 0.0f, tau = 0.0f;
    bondStress(particles[b.i], particles[b.j], b, simCfg.E, simCfg.nu, sigma, tau);
    if (sigma > maxSigma) { maxSigma = sigma; }
    if (tau   > maxTau)   { maxTau   = tau;   maxTauBond = bi; }
    if (sigma > simCfg.tauC || tau > simCfg.tauC) {
      b.broken = true;
      ++count;
    }
  }
  (void)maxTauBond;
  if (frame <= 15) {
    std::printf("[sc%d frame%3d] maxSigma=%.3e maxTau=%.3e broken=%d\n",
      scaleIdx, frame, maxSigma, maxTau, count);
  }
  return count;
}

// ---------------------------------------------------------------------------
// Simulation::step
// One full implicit time step using manifold gradient descent.
// ---------------------------------------------------------------------------
void Simulation::step() {
  int   N  = (int)particles.size();
  float dt = simCfg.dt;

  // ---- 1. Predict: x_hat = x + v*dt, q_hat = geodesic step ----
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) {
      p_hat[i] = particles[i].pos;
      q_hat[i] = particles[i].rot;
      continue;
    }

    // Predicted position
    p_hat[i] = particles[i].pos + particles[i].vel * dt;

    // Predicted orientation via quaternion integration:
    // dq/dt = 0.5 * [angVel, 0] * q   (in pure-quaternion left-multiplication)
    glm::quat w_quat(0.0f,
                     particles[i].angVel.x * 0.5f,
                     particles[i].angVel.y * 0.5f,
                     particles[i].angVel.z * 0.5f);
    glm::quat q_dot = w_quat * particles[i].rot;
    glm::quat q_new(particles[i].rot.w + q_dot.w * dt,
                    particles[i].rot.x + q_dot.x * dt,
                    particles[i].rot.y + q_dot.y * dt,
                    particles[i].rot.z + q_dot.z * dt);
    q_hat[i] = glm::normalize(q_new);
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

  // ---- 4. Block coordinate descent: alternate PCG (positions) and Jacobi (rotations) ----
  // After each PCG block, rotations change, so positions must be re-solved.
  // After each Jacobi block, positions are re-solved. Repeat until joint convergence.
  const int outerLoops = 5;
  const int innerPCG   = std::max(1, simCfg.maxIter / outerLoops);
  const int innerJac   = 20;  // rotation Jacobi converges fast (rho_rot ~ 0.07)

  std::vector<glm::vec3> pos_save(N);

  for (int outer = 0; outer < outerLoops; ++outer) {
    // -- PCG: optimise positions with current rotations --
    solvePCG(is_load, load_target, innerPCG, simCfg.epsilon);

    // Enforce boundary conditions after PCG
    for (int i = 0; i < N; i++) {
      if (particles[i].fixed) {
        particles[i].pos = p_hat[i];
        particles[i].rot = q_hat[i];
      }
      if (is_load[i]) particles[i].pos = load_target[i];
    }

    // -- Jacobi: optimise rotations with current positions (save/restore positions) --
    for (int i = 0; i < N; i++) pos_save[i] = particles[i].pos;

    for (int iter = 0; iter < innerJac; ++iter) {
      gradientStep(1.0f);
      // Restore positions, keep rotation updates
      for (int i = 0; i < N; i++) particles[i].pos = pos_save[i];
      for (int i = 0; i < N; i++) {
        if (particles[i].fixed) particles[i].rot = q_hat[i];
      }
    }
  }

  // ---- 5. Update velocities from position change ----
  for (int i = 0; i < N; i++) {
    if (particles[i].fixed) {
      particles[i].vel    = glm::vec3(0.0f);
      particles[i].angVel = glm::vec3(0.0f);
      continue;
    }
    particles[i].vel = (particles[i].pos - p_hat[i]) / dt + particles[i].vel;
    // Approximate angular velocity from quaternion change
    // omega ~ 2 * q_dot * q_conj where q_dot ~ (q_new - q_hat) / dt
    glm::quat dq_quat(
      (particles[i].rot.w - q_hat[i].w) / dt,
      (particles[i].rot.x - q_hat[i].x) / dt,
      (particles[i].rot.y - q_hat[i].y) / dt,
      (particles[i].rot.z - q_hat[i].z) / dt
    );
    glm::quat q_conj(
       particles[i].rot.w,
      -particles[i].rot.x,
      -particles[i].rot.y,
      -particles[i].rot.z
    );
    glm::quat omega_quat = dq_quat * q_conj;
    particles[i].angVel = glm::vec3(2.0f * omega_quat.x,
                                    2.0f * omega_quat.y,
                                    2.0f * omega_quat.z);
  }

  // ---- 6. Check fracture ----
  int broken = checkFracture();
  if (broken > 0 && fractureFrame < 0) {
    fractureFrame = frame;
  }

  ++frame;
}
