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
void Simulation::init(const BeamConfig& bcfg, const SimConfig& scfg) {
  beamCfg = bcfg;
  simCfg  = scfg;
  frame        = 0;
  fractureFrame = -1;

  particles = generateBeam(bcfg);
  tagBoundaryParticles(particles, bcfg);
  bonds = generateBonds(particles);

  int N = (int)particles.size();
  p_hat.resize(N);
  q_hat.resize(N);
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
    // h_p gets both stretch (kn) and shear (kt) contributions
    if (!particles[b.i].fixed) { h_p[b.i] += kn + kt; h_q[b.i] += EIl; }
    if (!particles[b.j].fixed) { h_p[b.j] += kn + kt; h_q[b.j] += EIl; }
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
  int maxSigBond = -1, maxTauBond = -1;
  for (int bi = 0; bi < (int)bonds.size(); bi++) {
    auto& b = bonds[bi];
    if (b.broken) continue;
    float sigma = 0.0f, tau = 0.0f;
    bondStress(particles[b.i], particles[b.j], b, simCfg.E, simCfg.nu, sigma, tau);
    if (sigma > maxSigma) { maxSigma = sigma; maxSigBond = bi; }
    if (tau   > maxTau)   { maxTau   = tau;   maxTauBond = bi; }
    if (sigma > simCfg.tauC || tau > simCfg.tauC) {
      b.broken = true;
      ++count;
    }
  }
  if (frame <= 2) {
    std::printf("  [frame %d checkFracture] maxSigma=%.3e (bond %d)  maxTau=%.3e (bond %d)  tauC=%.3e  broken=%d\n",
      frame, maxSigma, maxSigBond, maxTau, maxTauBond, simCfg.tauC, count);
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

  // ---- 4. Iterative gradient descent ----
  int n_free = 0;
  for (int i = 0; i < N; i++) {
    if (!particles[i].fixed) ++n_free;
  }

  // With diagonal preconditioner, alpha=1 is the natural Jacobi step.
  // Line search will reduce it if needed for energy descent.
  float alpha_base = 1.0f;

  for (int iter = 0; iter < simCfg.maxIter; ++iter) {
    // Compute gradient norm to check convergence
    float gnorm = 0.0f;
    {
      // Quick gradient norm estimate (position part only for efficiency)
      std::vector<glm::vec3> gp(N, glm::vec3(0.0f));
      for (int i = 0; i < N; i++) {
        if (particles[i].fixed) continue;
        float m = particles[i].mass;
        gp[i] = (m / (dt * dt)) * (particles[i].pos - p_hat[i]);
      }
      for (const auto& b : bonds) {
        if (b.broken) continue;
        BondGradient bg = bondGradient(particles[b.i], particles[b.j],
                                       b, simCfg.E, simCfg.nu);
        gp[b.i] += bg.grad_pi;
        gp[b.j] += bg.grad_pj;
      }
      for (int i = 0; i < N; i++) {
        if (!particles[i].fixed) gnorm += glm::dot(gp[i], gp[i]);
      }
      gnorm = std::sqrt(gnorm);
    }

    if (n_free > 0 && gnorm < simCfg.epsilon * n_free) break;

    // Armijo backtracking line search
    float E0 = totalEnergy();
    // Add inertia contribution to energy (for line search)
    for (int i = 0; i < N; i++) {
      if (particles[i].fixed) continue;
      float m  = particles[i].mass;
      float r  = particles[i].radius;
      float Iq = (8.0f / 5.0f) * m * r * r;
      glm::vec3 dp = particles[i].pos - p_hat[i];
      E0 += 0.5f * (m / (dt * dt)) * glm::dot(dp, dp);
      glm::vec4 q_curr(particles[i].rot.x, particles[i].rot.y,
                       particles[i].rot.z, particles[i].rot.w);
      glm::vec4 q_hat_v(q_hat[i].x, q_hat[i].y, q_hat[i].z, q_hat[i].w);
      glm::vec4 dq = q_curr - q_hat_v;
      E0 += 0.5f * (Iq / (dt * dt)) * glm::dot(dq, dq);
    }

    // Save state before step
    std::vector<glm::vec3> pos_saved(N);
    std::vector<glm::quat> rot_saved(N);
    for (int i = 0; i < N; i++) {
      pos_saved[i] = particles[i].pos;
      rot_saved[i] = particles[i].rot;
    }

    float alpha = alpha_base;
    bool accepted = false;
    for (int bt = 0; bt < 10; ++bt) {
      // Apply gradient step
      gradientStep(alpha);

      // Clamp fixed particles
      for (int i = 0; i < N; i++) {
        if (particles[i].fixed) {
          particles[i].pos = pos_saved[i];
          particles[i].rot = rot_saved[i];
        }
        // Enforce load particle kinematic constraint
        if (is_load[i]) {
          particles[i].pos = load_target[i];
        }
      }

      // Check energy
      float E1 = totalEnergy();
      for (int i = 0; i < N; i++) {
        if (particles[i].fixed) continue;
        float m  = particles[i].mass;
        float r  = particles[i].radius;
        float Iq = (8.0f / 5.0f) * m * r * r;
        glm::vec3 dp = particles[i].pos - p_hat[i];
        E1 += 0.5f * (m / (dt * dt)) * glm::dot(dp, dp);
        glm::vec4 q_curr(particles[i].rot.x, particles[i].rot.y,
                         particles[i].rot.z, particles[i].rot.w);
        glm::vec4 q_hat_v(q_hat[i].x, q_hat[i].y, q_hat[i].z, q_hat[i].w);
        glm::vec4 dq = q_curr - q_hat_v;
        E1 += 0.5f * (Iq / (dt * dt)) * glm::dot(dq, dq);
      }

      if (E1 <= E0 + 1e-8f * std::abs(E0)) {
        accepted = true;
        break;
      }

      // Reject step: restore state and halve alpha
      for (int i = 0; i < N; i++) {
        particles[i].pos = pos_saved[i];
        particles[i].rot = rot_saved[i];
      }
      alpha *= 0.5f;
    }

    // If no alpha was accepted, accept the smallest step anyway to avoid stall
    if (!accepted) {
      gradientStep(alpha);
      for (int i = 0; i < N; i++) {
        if (particles[i].fixed) {
          particles[i].pos = pos_saved[i];
          particles[i].rot = rot_saved[i];
        }
        if (is_load[i]) {
          particles[i].pos = load_target[i];
        }
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
