#pragma once
#include <vector>
#include "Particle.h"
#include "Bond.h"
#include "HexPacking.h"

// Simulation configuration parameters
struct SimConfig {
  float dt;        // timestep (s)
  float E;         // Young's modulus (Pa)
  float nu;        // Poisson's ratio
  float tauC;      // fracture stress threshold (Pa)
  float loadVel;   // load point downward velocity (m/s)
  float gravity;   // gravitational acceleration (m/s^2, 0 for quasi-static bending test)
  int   maxIter;   // max gradient descent iterations per time step
  float epsilon;   // convergence threshold for gradient norm (per DOF)
};

// Main simulation class implementing implicit manifold gradient descent (BDEM)
class Simulation {
public:
  std::vector<Particle> particles;
  std::vector<Bond>     bonds;
  BeamConfig            beamCfg;
  SimConfig             simCfg;
  int                   frame;          // current frame number
  int                   fractureFrame;  // frame when first bond broke (-1 = not yet)

  // Initialize with beam and simulation configurations
  void init(const BeamConfig& bcfg, const SimConfig& scfg);

  // Advance simulation by one timestep (calls gradient descent + fracture check)
  void step();

  // Compute total potential energy of all active bonds
  float totalEnergy() const;

  // Apply one gradient descent step on SO(3) x R^3 manifold.
  // alpha: step size multiplier
  void gradientStep(float alpha);

  // Check all bonds for fracture; break those exceeding tauC.
  // Returns number of bonds broken in this call.
  int checkFracture();

private:
  // Predicted positions and orientations (inertia terms)
  std::vector<glm::vec3> p_hat;
  std::vector<glm::quat> q_hat;
};
