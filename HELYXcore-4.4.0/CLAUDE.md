# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

HELYXcore 4.4.0 is an open-source CFD (Computational Fluid Dynamics) library based on OpenFOAM, developed by ENGYS Ltd. It provides solvers and utilities for various CFD applications including incompressible, compressible, multiphase, heat transfer, and reacting flows.

## Build Commands

### Initial Configuration
```bash
# Standard configuration
cmake -S . -B build

# Release build (optimized)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Debug build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

# With custom compiler
cmake -S . -B build -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
```

### Building
```bash
# Build all targets with parallel compilation
cmake --build build -j
cmake --build build -j$(nproc)

# Build specific target
cmake --build build --target <target_name>

# List available targets
cmake --build build --target help

# Clean build
cmake --build build --target clean

# Complete rebuild
rm -rf build
cmake -S . -B build
cmake --build build -j
```

### Environment Setup
```bash
# Source the environment after building (required for running)
source build/bin/helyxcore-bashrc

# Verify environment
echo $HELYXCORE_DIR  # Should show the build directory path
```

## Testing Commands

### Running Examples
```bash
# Run all examples (comprehensive test)
cd examples
./Allrun

# Clean example results
./Allclean

# Run specific solver examples
cd examples/helyxSolve/<category>
./Allrun

# Manual execution
cd examples/helyxSolve/<category>/<case>
helyxMesh      # Generate mesh
helyxSolve     # Run solver
helyxReport    # Generate reports
```

### Test Categories
- **helyxSolve/**: Main solver test cases
  - incompressible/: Steady and transient incompressible flows
  - compressible/: Compressible flow problems
  - multiphase/: VOF and mixture model cases
  - heatTransfer/: Conjugate heat transfer
  - turbulence/: Various turbulence models

### Quick Validation
```bash
# After sourcing environment
helyxSolve -help  # Should display solver options

# Run a simple test case
cd examples/helyxSolve/incompressible/cavity
helyxMesh
helyxSolve
```

## Code Architecture

### Directory Structure
```
HELYXcore-4.4.0/
├── src/                  # Core libraries
│   ├── HELYXcore/       # Base classes and utilities
│   ├── finiteVolume/    # FV discretization
│   ├── meshTools/       # Mesh manipulation
│   ├── solvers/         # Solver implementations
│   └── turbulence/      # Turbulence models
├── applications/        # Executable programs
│   ├── solvers/         # CFD solver applications
│   │   └── helyxSolve/  # Main multi-physics solver
│   └── utilities/       # Pre/post-processing tools
│       ├── helyxMesh/   # Mesh generation
│       └── helyxReport/ # Result analysis
├── modules/             # Optional functionality
├── thirdParty/          # External dependencies
├── examples/            # Test cases and tutorials
├── bin/                 # Scripts and utilities
├── doc/                 # Documentation
└── etc/                 # Configuration files
```

### Key Components

**helyxSolve**: Unified multi-physics solver
- Incompressible/compressible flows
- Heat transfer and conjugate heat transfer
- Multiphase flows (VOF, mixture model)
- Reacting flows and combustion
- Moving mesh and dynamic mesh refinement

**helyxMesh**: Advanced meshing utility
- Cartesian cut-cell meshing
- Automatic refinement regions
- Boundary layer generation
- STL surface handling

**helyxReport**: Post-processing and reporting
- Force coefficients
- Flow rate calculations
- Average field values
- Custom function objects

### CMake Build System
```cmake
# Key CMake files
CMakeLists.txt           # Main configuration
src/CMakeLists.txt       # Library definitions
applications/CMakeLists.txt  # Application targets

# Build artifacts location
build/
├── bin/                 # Executables and scripts
├── lib/                 # Compiled libraries
└── helyxcore-bashrc     # Environment setup script
```

### Development Patterns

**Solver Structure**:
```cpp
// Typical solver loop in helyxSolve
while (runTime.run())
{
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"
    
    runTime++;
    
    // Solver-specific equations
    #include "UEqn.H"
    #include "pEqn.H"
    
    runTime.write();
}
```

**Field Definitions**:
```cpp
volScalarField p(...);     // Pressure
volVectorField U(...);     // Velocity
surfaceScalarField phi(...); // Face flux
```

**Custom Boundary Conditions**:
- Located in `src/finiteVolume/fields/fvPatchFields/`
- Inherit from base fvPatchField classes
- Registered via runtime selection

## Development Workflow

### Adding New Features
1. Create source files in appropriate `src/` subdirectory
2. Update `CMakeLists.txt` to include new sources
3. Rebuild: `cmake --build build -j`
4. Create test case in `examples/`
5. Validate with example run

### Debugging
```bash
# Debug build
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j

# Run with debugger
source build-debug/bin/helyxcore-bashrc
gdb helyxSolve
run -case <case-directory>

# Enable detailed output
helyxSolve -case <dir> -verbose

# Check for floating point errors
export FOAM_SIGFPE=true
helyxSolve -case <dir>
```

### Running Parallel Cases
```bash
# Decompose domain
helyxDecompose -case <dir> -n 4

# Run in parallel
mpirun -np 4 helyxSolve -parallel -case <dir>

# Reconstruct results
helyxReconstruct -case <dir>
```

### Testing Changes
```bash
# Run specific example to test changes
cd examples/helyxSolve/<relevant-category>/<test-case>
./Allclean  # Clean previous results
./Allrun    # Run test

# Compare with reference solution
diff -r <results> <reference-results>
```

## Important Notes

1. **No environment conflict**: HELYXcore uses CMake, so it doesn't conflict with OpenFOAM environment variables
2. **Build directory**: Always use out-of-source builds (build/ directory)
3. **Environment setup**: Must source `build/bin/helyxcore-bashrc` before running applications
4. **Parallel builds**: Use `-j` flag with cmake for faster compilation
5. **Code style**: Follow OpenFOAM conventions (4 spaces, CamelCase classes, .H/.C file pairs)
6. **Testing**: Always validate changes with relevant examples
7. **Third-party deps**: Located in `thirdParty/`, built automatically with main build
8. **Module loading**: Optional modules in `modules/` can be enabled/disabled via CMake options
9. **STL files**: Place geometry files in case `constant/geometry/` directory
10. **Clean rebuilds**: Remove entire `build/` directory for complete clean rebuild

## Cross-References

- **[Root CLAUDE.md](../CLAUDE.md)**: Multi-project overview and comparison
- **[OpenFOAM-13](../OpenFOAM-13/CLAUDE.md)**: Foundation version with modular solvers
- **[OpenFOAM-v2506](../OpenFOAM-v2506/CLAUDE.md)**: ESI-OpenCFD version with extended features
- **[Root AGENTS.md](../AGENTS.md)**: Development guidelines and commit conventions
- **[Project README](README.md)**: HELYXcore specific information
- **[Examples](examples/)**: Test cases and validation examples