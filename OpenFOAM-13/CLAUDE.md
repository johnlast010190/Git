# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

OpenFOAM-13 is the Foundation release of the OpenFOAM open-source CFD (Computational Fluid Dynamics) toolbox. Key features include:
- **Modular solver framework**: `foamRun`/`foamMultiRun` with runtime-loadable solver modules
- **wmake build system**: Custom make system with dependency resolution
- **Comprehensive test suite**: Unit tests in `test/` and functional tests via `tutorials/`

## Build Commands

### Environment Setup (Required First)
```bash
# Source environment - REQUIRED for every new shell
source etc/bashrc

# Check environment is set
echo $WM_PROJECT_DIR  # Should output the OpenFOAM-13 directory path
foamVersion           # Display version information
```

### Building
```bash
# Full build (may need 2 passes for circular dependencies)
./Allwmake -j
./Allwmake -j  # Run again if first pass has unresolved dependencies

# Build with specific number of cores
./Allwmake -j$(nproc)
./Allwmake -j8

# Build specific library
cd src/<library>
wmake
wmake libso  # Explicitly build shared library

# Build specific solver/utility
cd applications/solvers/<solver>
wmake

# Build solver module
cd applications/modules/<module>
wmake

# Force rebuild
wmake -a  # Rebuild all dependencies
wclean && wmake  # Clean then build
```

### Cleaning
```bash
# Clean current directory
wclean

# Clean all (removes all object files and binaries)
wclean all

# Clean specific platform
wcleanPlatform

# Clean libraries only
wclean libso

# Clean and rebuild from scratch
wclean all && ./Allwmake -j

# Clean tutorials
cd tutorials && ./Allclean
```

## Test Commands

### Unit Tests
```bash
# Run all tests
cd test
./Allrun

# Run specific test category
cd test/<test-name>
./Allrun

# Available test categories:
# - dictionary: Dictionary and I/O tests
# - fvMeshTools: Mesh manipulation tests
# - meshTools: Mesh utility tests
# - parallel: Parallel communication tests

# Clean test artifacts
cd test
./Allclean
```

### Tutorial Cases (Functional Tests)
```bash
# Run all tutorials (comprehensive but time-consuming)
cd tutorials
./Alltest

# Run specific tutorial
cd tutorials/<category>/<solver>/<case>
./Allrun

# Manual tutorial execution
cd tutorials/<category>/<solver>/<case>
blockMesh           # Generate mesh
checkMesh           # Verify mesh quality
<solver>            # Run solver
foamRun             # For modular solver cases
paraFoam            # Visualize results
foamLog log         # Extract residuals

# Quick validation test
cd tutorials/incompressibleFluid/pitzDaily
./Allrun

# Clean tutorial results
cd tutorials/<category>/<solver>/<case>
./Allclean
foamCleanTutorials  # Clean all subdirectories
```

### System Validation
```bash
# Check OpenFOAM installation
foamSystemCheck
foamInstallationTest

# Display system information
foamVersion -full
foamVersion -api

# Create user run directory
mkdir -p $FOAM_RUN
run  # Changes to ~/OpenFOAM/<user>-13/run
```

## Code Architecture

### Core Structure
```
OpenFOAM-13/
├── src/                          # Core libraries
│   ├── OpenFOAM/                # Base classes, containers, I/O
│   ├── finiteVolume/            # FV discretization, fields, BCs
│   ├── meshTools/               # Mesh manipulation
│   ├── thermophysicalModels/    # Thermophysics, equations of state
│   ├── MomentumTransportModels/ # Turbulence and transport
│   ├── lagrangian/              # Particle tracking
│   └── functionObjects/         # Runtime post-processing
├── applications/                 # Executables
│   ├── modules/                 # Modular solver components
│   │   ├── fluid/              # Fluid flow modules
│   │   ├── solid/              # Solid mechanics
│   │   └── multiphaseEuler/    # Multiphase modules
│   ├── solvers/                 # Standalone solvers
│   │   ├── foamRun             # Modular solver executor
│   │   └── foamMultiRun        # Multi-region modular
│   └── utilities/               # Pre/post tools
│       ├── mesh/               # Mesh utilities
│       └── postProcessing/     # Analysis tools
├── wmake/                       # Build system
│   ├── rules/                  # Platform-specific rules
│   ├── scripts/                # Build scripts
│   └── wmake                   # Main build script
├── test/                        # Unit tests
└── tutorials/                   # Example cases
```

### Key Patterns

**Field Operations**:
```cpp
// Field declarations
volScalarField p(...);           // Pressure field
volVectorField U(...);           // Velocity field
surfaceScalarField phi(...);     // Face flux field

// Field operations
fvScalarMatrix pEqn(...);        // Pressure equation
solve(UEqn == -fvc::grad(p));    // Solve momentum
p.correctBoundaryConditions();   // Update BCs
```

**Boundary Conditions**:
```cpp
// In field dictionary
boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0);
    }
    wall
    {
        type            noSlip;
    }
}
```

**Runtime Selection Tables (RTS)**:
```cpp
// Usage example
autoPtr<turbulenceModel> turbulence
(
    turbulenceModel::New(U, phi, transport)
);
```

**Solver Modules** (OpenFOAM-13 specific):
```cpp
// In controlDict
solver          fluid;  // Loads fluid module

// Module structure
class fluid : public solver
{
    // Implementation of solver phases
    virtual void preSolve();
    virtual void momentumPredictor();
    virtual void pressureCorrector();
    virtual void postSolve();
};
```

### wmake Build System

**Make/files**: Lists source files
```make
source1.C
source2.C

LIB = $(FOAM_LIBBIN)/libmyLibrary
```

**Make/options**: Compiler flags and dependencies
```make
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude

LIB_LIBS = \
    -lfiniteVolume \
    -lmeshTools
```

## Development Workflow

### Creating a New Solver Module
```bash
# Copy template
cd applications/modules
cp -r fluid myModule
cd myModule

# Rename files and class
mv fluid.C myModule.C
mv fluid.H myModule.H

# Update Make/files
sed -i 's/fluid/myModule/g' Make/files

# Edit source files to implement functionality
# Build
wmake

# Use in case via controlDict:
# solver myModule;
```

### Adding a Boundary Condition
```bash
# Copy template BC
cd src/finiteVolume/fields/fvPatchFields
cp -r derived/fixedValue myBC

# Implement in myBC.C and myBC.H
# Add to Make/files
# Register with TypeName macro
wmake libso
```

### Running Parallel Cases
```bash
# Setup parallel decomposition
cat > system/decomposeParDict << EOF
numberOfSubdomains 4;
method          scotch;
EOF

# Decompose mesh
decomposePar

# Run in parallel
mpirun -np 4 <solver> -parallel
# Or for modular solver:
mpirun -np 4 foamRun -parallel

# Monitor progress
tail -f log

# Reconstruct results
reconstructPar
reconstructPar -latestTime  # Only latest time

# Clean parallel files
rm -rf processor*
```

### Debugging
```bash
# Debug build
export WM_COMPILE_OPTION=Debug
source etc/bashrc
./Allwmake -j

# Run with debugger
gdb <solver>
(gdb) run -case <dir>

# Enable floating point exceptions
export FOAM_SIGFPE=true

# Increase verbosity
<solver> -case <dir> -v

# Enable debug switches in controlDict
DebugSwitches
{
    <className> 1;  // Enable debug output
}

# Profile performance
<solver> -case <dir> -profiling
```

### Code Style Guidelines
```cpp
// Header guard format
#ifndef myClass_H
#define myClass_H

// Namespace
namespace Foam
{

// Class documentation
/*---------------------------------------------------------------------------*\
                         Class myClass Declaration
\*---------------------------------------------------------------------------*/

class myClass
:
    public baseClass
{
    // Private data
        scalar myData_;

public:

    // Constructors
        myClass();

    // Member Functions
        scalar calculate() const;
};

} // End namespace Foam

#endif
```

## Common Operations

### Finding Examples
```bash
# Find tutorial cases using specific solver
find tutorials -name controlDict -exec grep -l "application.*<solver>" {} \;

# Find usage of specific boundary condition
grep -r "type.*<BC-name>" tutorials/

# List available solvers
ls -1 $FOAM_APPBIN

# List available libraries
ls -1 $FOAM_LIBBIN
```

### Environment Variables
```bash
# Key variables
$WM_PROJECT_DIR    # OpenFOAM-13 root
$FOAM_SRC          # Source code
$FOAM_APP          # Applications
$FOAM_TUTORIALS    # Tutorial cases
$FOAM_RUN          # User run directory
$FOAM_USER_APPBIN  # User applications
$FOAM_USER_LIBBIN  # User libraries

# Compilation settings
$WM_COMPILE_OPTION # Debug/Opt
$WM_COMPILER       # Gcc/Clang/Icc
$WM_PRECISION_OPTION # SP/DP
$WM_MPLIB          # SYSTEMOPENMPI/OPENMPI/MPICH
```

### Function Objects
```cpp
// In controlDict
functions
{
    forces
    {
        type            forces;
        libs            ("libforces.so");
        patches         (wall);
        rho             rhoInf;
        rhoInf          1.225;
        CofR            (0 0 0);
    }
    
    fieldAverage
    {
        type            fieldAverage;
        libs            ("libfieldFunctionObjects.so");
        fields
        (
            U
            {
                mean        on;
                prime2Mean  on;
                base        time;
            }
        );
    }
}
```

## Important Notes

1. **Environment isolation**: Never source multiple OpenFOAM environments in same shell
2. **Build dependencies**: May need 2 `./Allwmake` passes for circular dependencies
3. **Testing changes**: Always validate with relevant tutorial case
4. **Code style**: Follow existing patterns - 4 spaces, no tabs, Foam:: namespace
5. **Parallel builds**: Use `-j` flag for faster compilation
6. **Clean builds**: Use `wclean all` when switching configurations
7. **Debug builds**: Set `WM_COMPILE_OPTION=Debug` before sourcing bashrc
8. **Modular solvers**: Preferred over monolithic solvers for new developments
9. **Function objects**: Use for runtime post-processing instead of separate utilities
10. **Version control**: Never commit platforms/, lnInclude/, or *.dep files

## Cross-References

- **[Root CLAUDE.md](../CLAUDE.md)**: Multi-project overview and comparison
- **[OpenFOAM-v2506](../OpenFOAM-v2506/CLAUDE.md)**: ESI-OpenCFD version with extended features
- **[HELYXcore](../HELYXcore-4.4.0/CLAUDE.md)**: CMake-based unified solver framework
- **[Root AGENTS.md](../AGENTS.md)**: Development guidelines and commit conventions
- **[Project README](README.md)**: OpenFOAM-13 specific information