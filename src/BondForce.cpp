#include "BondForce.h"

#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Gl(q): 3x4 nullspace operator (left quaternion multiplication helper).
// Returns a mat<4,3,float> in GLM column-major convention, i.e.
//   result[col][row]
// Mathematical form (3 rows x 4 cols):
//   [ w,  z, -y, -x ]
//   [-z,  w,  x, -y ]
//   [ y, -x,  w, -z ]
// Here q = (x,y,z,w) in GLM storage order.
static glm::mat<4, 3, float> Gl(const glm::quat& q) {
  float x = q.x, y = q.y, z = q.z, w = q.w;
  glm::mat<4, 3, float> m;
  // col 0 (d/dq.x): row0= w, row1=-z, row2= y
  m[0] = glm::vec<3, float>( w, -z,  y);
  // col 1 (d/dq.y): row0= z, row1= w, row2=-x
  m[1] = glm::vec<3, float>( z,  w, -x);
  // col 2 (d/dq.z): row0=-y, row1= x, row2= w
  m[2] = glm::vec<3, float>(-y,  x,  w);
  // col 3 (d/dq.w): row0=-x, row1=-y, row2=-z
  m[3] = glm::vec<3, float>(-x, -y, -z);
  return m;
}

// Gr(q): 3x4 nullspace operator (right quaternion multiplication helper).
// Mathematical form (3 rows x 4 cols):
//   [ w, -z,  y, -x ]
//   [ z,  w, -x, -y ]
//   [-y,  x,  w, -z ]
static glm::mat<4, 3, float> Gr(const glm::quat& q) {
  float x = q.x, y = q.y, z = q.z, w = q.w;
  glm::mat<4, 3, float> m;
  m[0] = glm::vec<3, float>( w,  z, -y);
  m[1] = glm::vec<3, float>(-z,  w,  x);
  m[2] = glm::vec<3, float>( y, -x,  w);
  m[3] = glm::vec<3, float>(-x, -y, -z);
  return m;
}

// Standard rotation matrix from unit quaternion q.
static glm::mat3 rotMat(const glm::quat& q) {
  return glm::mat3_cast(q);
}

// Rotate vector v by quaternion q (equivalent to R(q)*v).
static glm::vec3 rotVec(const glm::quat& q, const glm::vec3& v) {
  return q * v;
}

// Jacobian of Gl(q)*u with respect to q, for fixed u (vec4, (x,y,z,w) order).
// The result is a 3x4 matrix (GLM: mat<4,3,float>) whose columns are
// d(Gl(q)*u)/d(q.x), d(Gl(q)*u)/d(q.y), d(Gl(q)*u)/d(q.z), d(Gl(q)*u)/d(q.w).
//
// From the definition of Gl:
//   Gl(q)*u = [  w*u1 + z*u2 - y*u3 - x*u4,
//               -z*u1 + w*u2 + x*u3 - y*u4,
//                y*u1 - x*u2 + w*u3 - z*u4 ]
// where u = (u1, u2, u3, u4) = (u.x, u.y, u.z, u.w).
// Differentiating:
//   d/d(q.x): [-u4, u3, -u2]
//   d/d(q.y): [-u3,-u4,  u1]
//   d/d(q.z): [ u2,-u1, -u4]
//   d/d(q.w): [ u1, u2,  u3]
static glm::mat<4, 3, float> jacGl(const glm::vec4& u) {
  glm::mat<4, 3, float> m;
  // col 0 = d/d(q.x)
  m[0] = glm::vec<3, float>(-u.w,  u.z, -u.y);
  // col 1 = d/d(q.y)
  m[1] = glm::vec<3, float>(-u.z, -u.w,  u.x);
  // col 2 = d/d(q.z)
  m[2] = glm::vec<3, float>( u.y, -u.x, -u.w);
  // col 3 = d/d(q.w)
  m[3] = glm::vec<3, float>( u.x,  u.y,  u.z);
  return m;
}

// ---------------------------------------------------------------------------
// Stiffness parameters for a bond between particles pi and pj
// ---------------------------------------------------------------------------
struct BondStiffness {
  float S;   // cross-section area
  float kn;  // normal (stretch) stiffness
  float kt;  // shear stiffness
  float r0;  // mean radius
  float I;   // second moment of area
  float J;   // polar moment of area
  glm::mat3 K;  // 3x3 bending+torsion stiffness diagonal
};

static BondStiffness computeStiffness(const Particle& pi, const Particle& pj,
                                      const Bond& b, float E, float nu) {
  BondStiffness st;
  st.r0 = (pi.radius + pj.radius) * 0.5f;
  float pi_const = glm::pi<float>();
  st.S  = pi_const * st.r0 * st.r0;
  st.kn = E * st.S / b.l0;
  float G = E / (2.0f * (1.0f + nu));
  st.kt = G * st.S / b.l0;
  st.I  = pi_const * st.r0 * st.r0 * st.r0 * st.r0 / 4.0f;
  st.J  = pi_const * st.r0 * st.r0 * st.r0 * st.r0 / 2.0f;
  // K = diag(G*J/l0, E*I/l0, E*I/l0)
  float Kt_val  = G * st.J / b.l0;
  float Kb_val  = E * st.I / b.l0;
  st.K = glm::mat3(0.0f);
  st.K[0][0] = Kt_val;
  st.K[1][1] = Kb_val;
  st.K[2][2] = Kb_val;
  return st;
}

// ---------------------------------------------------------------------------
// Stretch energy
// ---------------------------------------------------------------------------

static float stretchEnergy(const Particle& pi, const Particle& pj,
                            const Bond& b, float kn) {
  float dist  = glm::length(pj.pos - pi.pos);
  float delta = dist - b.l0;
  return 0.5f * kn * delta * delta;
}

// ---------------------------------------------------------------------------
// Shear energy (Eq. 9 in paper)
// V_shear = (1/2) kt l0^2 theta^2
// theta = angle(d_actual, qc⊙d0)
// d_actual = (pj - pi) / |pj - pi|  is the CURRENT physical bond direction.
// qc = (qi+qj)/||qi+qj|| is the mean orientation.
// d_exp = qc⊙d0 is the "expected direction" tracked by the mean orientation.
// Shear energy is zero when orientations correctly predict the bond direction.
// This creates position-orientation coupling: bending positions drive orientations.
// ---------------------------------------------------------------------------
static float shearEnergy(const Particle& pi, const Particle& pj,
                          const Bond& b, float kt) {
  glm::vec3 dv   = pj.pos - pi.pos;
  float     dist = glm::length(dv);
  if (dist < 1e-12f) return 0.0f;
  glm::vec3 d_actual = dv / dist;

  glm::vec4 qi_vec(pi.rot.x, pi.rot.y, pi.rot.z, pi.rot.w);
  glm::vec4 qj_vec(pj.rot.x, pj.rot.y, pj.rot.z, pj.rot.w);
  glm::vec4 s_sum  = qi_vec + qj_vec;
  float     norm_s = glm::length(s_sum);
  if (norm_s < 1e-12f) return 0.0f;
  glm::quat qc(s_sum.w / norm_s, s_sum.x / norm_s,
               s_sum.y / norm_s, s_sum.z / norm_s);

  // d_exp = qc⊙d0; theta = angle(d_actual, d_exp)
  glm::vec3 d_exp    = rotVec(qc, b.d0);
  float cos_theta    = glm::clamp(glm::dot(d_actual, d_exp), -1.0f, 1.0f);
  float theta        = std::acos(cos_theta);
  return 0.5f * kt * b.l0 * b.l0 * theta * theta;
}

// ---------------------------------------------------------------------------
// Bend-twist energy (Eq. 10-11 in paper)
// ---------------------------------------------------------------------------
static float bendTwistEnergy(const Particle& pi, const Particle& pj,
                               const Bond& b, const BondStiffness& st) {
  // R0 = rotation matrix from q0 (maps (1,0,0) -> d0)
  glm::mat3 R0 = rotMat(b.q0);

  // Gl(qj) applied to vec4(qi)
  glm::vec4 qi_vec(pi.rot.x, pi.rot.y, pi.rot.z, pi.rot.w);
  glm::mat<4, 3, float> Glqj = Gl(pj.rot);
  // t = R0 * (Gl(qj) * qi_vec)   result is 3-vector
  // Gl(qj) is (4 cols) x (3 rows), qi_vec is vec4 -> result is vec3
  glm::vec3 Glqj_times_qi = Glqj * qi_vec;  // mat<4,3> * vec4 = vec3
  glm::vec3 t = R0 * Glqj_times_qi;

  return 0.5f * glm::dot(t, st.K * t);
}

// ---------------------------------------------------------------------------
// bondEnergy
// ---------------------------------------------------------------------------
float bondEnergy(const Particle& pi, const Particle& pj, const Bond& b,
                 float E, float nu) {
  BondStiffness st = computeStiffness(pi, pj, b, E, nu);
  float Vs = stretchEnergy(pi, pj, b, st.kn);
  float Vsh = shearEnergy(pi, pj, b, st.kt);
  float Vbt = bendTwistEnergy(pi, pj, b, st);
  return Vs + Vsh + Vbt;
}

// ---------------------------------------------------------------------------
// bondGradient
// ---------------------------------------------------------------------------
BondGradient bondGradient(const Particle& pi, const Particle& pj,
                           const Bond& b, float E, float nu) {
  BondStiffness st = computeStiffness(pi, pj, b, E, nu);
  BondGradient bg;
  bg.grad_pi = glm::vec3(0.0f);
  bg.grad_pj = glm::vec3(0.0f);
  bg.grad_qi = glm::vec4(0.0f);
  bg.grad_qj = glm::vec4(0.0f);

  // ---- Stretch gradient ----
  glm::vec3 d    = pj.pos - pi.pos;
  float     dist = glm::length(d);
  if (dist > 1e-12f) {
    glm::vec3 dhat  = d / dist;
    float     delta = dist - b.l0;
    bg.grad_pi += -st.kn * delta * dhat;
    bg.grad_pj +=  st.kn * delta * dhat;
  }

  // ---- Shear gradient ----
  // V_shear = (1/2)*kt*l0^2*theta^2, theta = angle(d_actual, qc⊙d0)   [Paper Eq. 9]
  // d_actual = (pj-pi)/|pj-pi| is the current physical bond direction.
  // Gradient w.r.t. positions (transverse shear force) AND orientations.
  if (dist > 1e-12f) {
    glm::vec3 d_actual = d / dist;  // d = pj - pi computed above

    glm::vec4 qi_vec(pi.rot.x, pi.rot.y, pi.rot.z, pi.rot.w);
    glm::vec4 qj_vec(pj.rot.x, pj.rot.y, pj.rot.z, pj.rot.w);
    glm::vec4 s_sum  = qi_vec + qj_vec;
    float     norm_s = glm::length(s_sum);
    if (norm_s > 1e-12f) {
      glm::quat qc(s_sum.w / norm_s, s_sum.x / norm_s,
                   s_sum.y / norm_s, s_sum.z / norm_s);
      glm::vec3 d_exp    = rotVec(qc, b.d0);
      float cos_theta    = glm::clamp(glm::dot(d_actual, d_exp), -1.0f, 1.0f);
      float theta        = std::acos(cos_theta);
      float sin_theta    = std::sin(theta);

      if (sin_theta > 1e-6f) {
        float factor = st.kt * b.l0 * b.l0 * theta / sin_theta;

        // ---- Position gradient ----
        // dV/d(pj) = factor * (-1/sin_theta) * d(d_actual . d_exp)/d(pj)
        // d(d_actual)/d(pj) = (I - d_actual x d_actual) / dist
        // => dV/d(pj) = -factor * (d_exp - cos_theta * d_actual) / dist
        glm::vec3 transverse = d_exp - cos_theta * d_actual;
        bg.grad_pj += -factor * transverse / dist;
        bg.grad_pi +=  factor * transverse / dist;

        // ---- Orientation gradient ----
        // d(cos_theta)/d(qc) where cos_theta = dot(d_actual, R(qc)*d0).
        // Using d/dq [v^T R(q) u] with v = d_actual, u = d0:
        //   vector part: 2(v.t3)u + 2(t3.u)v - 2(v.u)t3 + 2s(u x v)
        //   scalar part: 2s(v.u) + 2t3.(u x v)
        glm::vec3 t3(qc.x, qc.y, qc.z);
        float s_qc   = qc.w;
        glm::vec3 uxv    = glm::cross(b.d0, d_actual);  // d0 x d_actual
        float tdotu0     = glm::dot(t3, b.d0);
        float vdotu0     = glm::dot(d_actual, b.d0);

        glm::vec3 dcos_dt3 = 2.0f * glm::dot(d_actual, t3) * b.d0
                           + 2.0f * tdotu0 * d_actual
                           - 2.0f * vdotu0 * t3
                           + 2.0f * s_qc * uxv;
        float dcos_ds = 2.0f * s_qc * vdotu0 + 2.0f * glm::dot(t3, uxv);

        glm::vec4 dcos_dqc(dcos_dt3.x, dcos_dt3.y, dcos_dt3.z, dcos_ds);
        glm::vec4 dV_dqc = -factor * dcos_dqc;

        // Project onto tangent space of S3 at qc
        glm::vec4 qc_vec(qc.x, qc.y, qc.z, qc.w);
        glm::vec4 dV_proj = dV_dqc - glm::dot(qc_vec, dV_dqc) * qc_vec;

        // Chain rule: dqc/dqi = dqc/dqj = (I - qc*qc^T)/norm_s
        glm::vec4 grad_q_shear = dV_proj / norm_s;
        bg.grad_qi += grad_q_shear;
        bg.grad_qj += grad_q_shear;
      }
    }
  }

  // ---- Bend-twist gradient ----
  // V_bt = 0.5 * t^T * K * t
  // t = R0 * Gl(qj) * vec4(qi)
  {
    glm::mat3  R0      = rotMat(b.q0);
    glm::vec4  qi_vec(pi.rot.x, pi.rot.y, pi.rot.z, pi.rot.w);
    glm::mat<4, 3, float> Glqj = Gl(pj.rot);

    glm::vec3 Glqj_qi = Glqj * qi_vec;   // 3-vector
    glm::vec3 t       = R0 * Glqj_qi;    // 3-vector
    glm::vec3 Kt      = st.K * t;         // 3-vector

    // grad_qi V_bt = Gl(qj)^T * R0^T * (K*t)
    // Gl(qj) is mat<4,3,float> -> Gl(qj)^T is mat<3,4,float> in GLM = mat<3,4,float>
    // But GLM transpose(mat<4,3>) = mat<3,4>
    glm::mat<3, 4, float> GlT = glm::transpose(Glqj);
    glm::vec3 R0T_Kt    = glm::transpose(R0) * Kt;    // 3-vector
    bg.grad_qi         += GlT * R0T_Kt;               // vec4

    // grad_qj V_bt = jacGl(vec4(qi))^T * R0^T * (K*t)
    glm::mat<4, 3, float> jGl    = jacGl(qi_vec);
    glm::mat<3, 4, float> jGlT   = glm::transpose(jGl);
    bg.grad_qj                   += jGlT * R0T_Kt;   // vec4
  }

  return bg;
}

// ---------------------------------------------------------------------------
// bondStress
// ---------------------------------------------------------------------------
void bondStress(const Particle& pi, const Particle& pj, const Bond& b,
                float E, float nu, float& sigma, float& tau) {
  BondStiffness st = computeStiffness(pi, pj, b, E, nu);

  // Normal force from stretch gradient
  glm::vec3 d    = pj.pos - pi.pos;
  float     dist = glm::length(d);
  glm::vec3 Fn(0.0f);
  if (dist > 1e-12f) {
    glm::vec3 dhat  = d / dist;
    float     delta = dist - b.l0;
    Fn = st.kn * delta * dhat;
  }

  // Bending moment from bend-twist energy gradient
  glm::mat3  R0     = rotMat(b.q0);
  glm::vec4  qi_vec(pi.rot.x, pi.rot.y, pi.rot.z, pi.rot.w);
  glm::mat<4, 3, float> Glqj = Gl(pj.rot);
  glm::vec3 Glqj_qi = Glqj * qi_vec;
  glm::vec3 t       = R0 * Glqj_qi;
  glm::vec3 Kt      = st.K * t;

  // Kt[0] = torsion component, Kt[1], Kt[2] = bending components
  float Mt = std::abs(Kt[0]);                                          // torsion magnitude
  float Mb = std::sqrt(Kt[1] * Kt[1] + Kt[2] * Kt[2]);               // bending moment magnitude

  // Shear force: Fs = nabla_p V_shear (transverse component)
  // theta = angle(d_actual, qc⊙d0); Fs = kt*l0^2*theta/sin(theta)*(d_exp - cos*d_actual)/dist
  glm::vec3 Fs(0.0f);
  {
    glm::vec4 qi_vec_s(pi.rot.x, pi.rot.y, pi.rot.z, pi.rot.w);
    glm::vec4 qj_vec_s(pj.rot.x, pj.rot.y, pj.rot.z, pj.rot.w);
    glm::vec4 s_sum_s  = qi_vec_s + qj_vec_s;
    float     ns = glm::length(s_sum_s);
    if (ns > 1e-12f && dist > 1e-12f) {
      glm::quat qcs(s_sum_s.w/ns, s_sum_s.x/ns, s_sum_s.y/ns, s_sum_s.z/ns);
      glm::vec3 d_actual_s = (dist > 1e-12f) ? (d / dist) : glm::vec3(0.0f);
      glm::vec3 d_exp_s = rotVec(qcs, b.d0);
      float cth = glm::clamp(glm::dot(d_actual_s, d_exp_s), -1.0f, 1.0f);
      float tth = std::acos(cth);
      float sth = std::sin(tth);
      if (sth > 1e-6f) {
        float fac = st.kt * b.l0 * b.l0 * tth / sth;
        Fs = fac * (d_exp_s - cth * d_actual_s) / dist;
      }
    }
  }

  // Normal stress: maximum tensile stress at the outer fibre of the bond cross-section.
  // sigma = max(0, Fn/S + Mb*r0/I) where Fn is the *signed* axial force.
  // For compressive bonds (delta < 0) the axial compression cancels the bending term;
  // since r0 << H/2, top-half bonds remain near zero and do not fracture prematurely.
  // Crack-tip bonds enter tension as the crack opens, enabling natural upward propagation.
  float delta_axial = dist > 1e-12f ? (dist - b.l0) : 0.0f;
  float Fn_signed   = st.kn * delta_axial;   // positive = tension, negative = compression
  sigma = std::max(0.0f, Fn_signed / st.S + Mb * st.r0 / st.I);

  // Shear stress: transverse shear force + torsion
  tau = glm::length(Fs) / st.S + Mt * st.r0 / st.J;
}
