#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include "Particle.h"
#include "Bond.h"

// Gradient of bond energy w.r.t. the four DOFs of particles i and j.
// q gradient is in (x,y,z,w) component order matching glm::quat storage.
struct BondGradient {
  glm::vec3 grad_pi;   // d(V)/d(p_i)
  glm::vec3 grad_pj;   // d(V)/d(p_j)
  glm::vec4 grad_qi;   // d(V)/d(q_i)  (x,y,z,w order)
  glm::vec4 grad_qj;   // d(V)/d(q_j)  (x,y,z,w order)
};

// Total bond potential energy V = V_stretch + V_shear + V_bendtwist
float bondEnergy(const Particle& pi, const Particle& pj, const Bond& b,
                 float E, float nu);

// Gradient of bond energy w.r.t. particle positions and orientations
BondGradient bondGradient(const Particle& pi, const Particle& pj, const Bond& b,
                          float E, float nu);

// Compute fracture stress measures sigma (normal) and tau (shear) using Eq. 12-17
void bondStress(const Particle& pi, const Particle& pj, const Bond& b,
                float E, float nu, float& sigma, float& tau);
