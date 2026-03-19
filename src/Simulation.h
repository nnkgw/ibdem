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
  int                   scaleIdx;       // 0=Low,1=Middle,2=High (for diagnostics)

  // Initialize with beam and simulation configurations
  void init(const BeamConfig& bcfg, const SimConfig& scfg, int idx = 0);

  // Advance simulation by one timestep (calls gradient descent + fracture check)
  void step();

  // Compute total potential energy of all active bonds
  float totalEnergy() const;

  // Apply one gradient descent step on SO(3) x R^3 manifold.
  // alpha: step size multiplier
  void gradientStep(float alpha);

  // PCG solver for position DOFs (rotations updated by Jacobi afterwards).
  // Uses analytical stretch-only Hessian-vector product.
  // Returns residual norm after solving.
  float solvePCG(const std::vector<bool>& is_load,
                 const std::vector<glm::vec3>& load_target,
                 int maxIter, float epsilon);

  // Check all bonds for fracture; break those exceeding tauC.
  // Returns number of bonds broken in this call.
  int checkFracture();

private:
  // Predicted positions and orientations (inertia terms)
  std::vector<glm::vec3> p_hat;
  std::vector<glm::quat> q_hat;

  // PCG working vectors (allocated once in init)
  std::vector<glm::vec3> pcg_r, pcg_d, pcg_z, pcg_Hd;
};
