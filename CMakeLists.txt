cmake_minimum_required(VERSION 3.16)
project(ibdem CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect all source files under src/
file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE src)

# ---------- GLM (header-only) ----------
find_package(glm CONFIG QUIET)
if(glm_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE glm::glm)
  message(STATUS "GLM: found via config package")
else()
  find_path(GLM_INCLUDE_DIR
    NAMES glm/glm.hpp
    PATHS
      /usr/include
      /usr/local/include
      /opt/homebrew/include
      "$ENV{VCPKG_ROOT}/installed/x64-windows/include"
      "${CMAKE_SOURCE_DIR}/third_party"
  )
  if(GLM_INCLUDE_DIR)
    target_include_directories(${PROJECT_NAME} PRIVATE ${GLM_INCLUDE_DIR})
    message(STATUS "GLM: found headers at ${GLM_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR
      "GLM not found.\n"
      "  macOS  : brew install glm\n"
      "  Ubuntu : sudo apt install libglm-dev\n"
      "  Windows: vcpkg install glm\n")
  endif()
endif()

# ---------- OpenGL + GLUT ----------
if(APPLE)
  # macOS: use Apple system frameworks (no freeglut needed)
  find_library(OPENGL_FRAMEWORK OpenGL REQUIRED)
  find_library(GLUT_FRAMEWORK   GLUT   REQUIRED)
  target_link_libraries(${PROJECT_NAME} PRIVATE
    ${OPENGL_FRAMEWORK}
    ${GLUT_FRAMEWORK}
  )
  # Suppress deprecated OpenGL API warnings on macOS 10.14+
  target_compile_definitions(${PROJECT_NAME} PRIVATE GL_SILENCE_DEPRECATION)
  message(STATUS "GLUT: using Apple GLUT framework (macOS)")
else()
  # Windows / Ubuntu: use freeglut
  find_package(OpenGL REQUIRED)
  find_package(GLUT   REQUIRED)

  # OpenGL
  if(TARGET OpenGL::GL)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL)
  else()
    target_include_directories(${PROJECT_NAME} PRIVATE ${OPENGL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${OPENGL_LIBRARIES})
  endif()

  # GLUT (freeglut)
  if(TARGET GLUT::GLUT)
    target_link_libraries(${PROJECT_NAME} PRIVATE GLUT::GLUT)
  else()
    target_include_directories(${PROJECT_NAME} PRIVATE ${GLUT_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${GLUT_LIBRARIES})
  endif()

  if(WIN32)
    message(STATUS "GLUT: using freeglut (Windows) -- install via: vcpkg install freeglut")
  else()
    message(STATUS "GLUT: using freeglut (Linux)   -- install via: sudo apt install freeglut3-dev")
  endif()
endif()

# ---------- MSVC: treat source files as UTF-8 (suppresses C4819) ----------
if(MSVC)
  target_compile_options(${PROJECT_NAME} PRIVATE /utf-8)
endif()
