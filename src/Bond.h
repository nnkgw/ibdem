#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

// A bond connecting two discrete elements i and j.
// Stores the rest-state geometry needed to evaluate stretch, shear,
// bend and twist energies (Section 3.2 of the paper).
struct Bond {
  int i, j;        // indices into the particle array
  float l0;        // rest length  ||p_j - p_i|| at creation (m)
  glm::vec3 d0;    // unit rest direction from i to j
  glm::quat q0;    // rotation that maps (1,0,0) to d0;
                   // used to transform the bend/twist vector to the bond-local frame
  bool  broken;    // true after the bond has fractured
  float sigma = 0.0f;  // cached normal stress (Pa) — updated by checkFracture(), used for rendering
};
