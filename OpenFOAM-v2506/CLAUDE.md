# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

OpenFOAM-v2506 is the ESI-OpenCFD distribution of OpenFOAM (v2506 release, June 2025) with extended functionality through modules and plugins. This computational fluid dynamics (CFD) framework includes additional features not found in the Foundation version.

## Build Commands

### Environment Setup (Required First)
```bash
# Source environment - CRITICAL: only one OpenFOAM per shell
source etc/bashrc

# Verify environment
echo $WM_PROJECT_VERSION  # Should show v2506
foamVersion               # Display version details
foamSystemCheck           # Check system readiness
```

### Full Build
```bash
# Main build with options:
# -j: parallel compilation
# -s: silent (reduced output)
# -q: queue jobs properly
# -l: log output to file
./Allwmake -j -s -q -l

# Check build log
tail -f log.Allwmake

# Optional: Build extended modules (after main build)
./Allwmake-modules

# Optional: Build third-party plugins (after main build)
./Allwmake-plugins

# Complete build (all components)
./Allwmake -j -s -q -l && ./Allwmake-modules && ./Allwmake-plugins
```

### Component Building
```bash
# Build specific library
cd src/<library_name>
wmake
wmake libso  # Explicitly shared library

# Build specific application
cd applications/<category>/<name>
wmake

# Build module
cd modules/<module_name>
./Allwmake

# Force rebuild
wmake -a
wclean && wmake
```

### Cleaning
```bash
# Clean current directory
wclean

# Clean all recursively
wclean all

# Clean platform-specific builds
wcleanPlatform

# Clean modules
cd modules && ./Allwclean

# Clean plugins
cd plugins && ./Allwclean

# Complete clean rebuild
wclean all
./Allwmake -j -s -q -l
```

## Test Commands

### Quick Validation
```bash
# Installation test
foamInstallationTest

# Test a simple case
mkdir -p $FOAM_RUN && cd $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily ./
cd pitzDaily
blockMesh
simpleFoam
paraFoam  # Visualize results
```

### Tutorial Cases
```bash
# Run all tutorials (time-consuming)
cd tutorials
./Alltest

# Run specific tutorial
cd tutorials/<category>/<solver>/<case>
./Allrun

# Manual execution
blockMesh           # Generate mesh
checkMesh           # Verify mesh quality
setFields           # Initialize fields if needed
<solver>            # Run solver
postProcess -func <function>  # Post-processing

# Clean tutorial
./Allclean
foamCleanTutorials  # Clean all subdirectories
```

### Module Testing
```bash
# Test specific module
cd modules/<module>/tutorials
./Alltest

# Example: Test OpenQBMM module
cd modules/OpenQBMM/tutorials
./Allrun
```

## Code Architecture

### Directory Structure
```
OpenFOAM-v2506/
├── src/                        # Core libraries
│   ├── OpenFOAM/              # Base classes, I/O, containers
│   ├── finiteVolume/          # FV discretization, fields
│   ├── meshTools/             # Mesh utilities
│   ├── turbulenceModels/      # Turbulence modeling
│   ├── thermophysicalModels/  # Thermophysics
│   ├── functionObjects/       # Runtime post-processing
│   └── sampling/              # Data sampling
├── applications/              # Executable programs
│   ├── solvers/              # CFD solvers
│   │   ├── incompressible/   # Incompressible flow
│   │   ├── compressible/     # Compressible flow
│   │   ├── multiphase/       # Multiphase solvers
│   │   └── heatTransfer/     # Heat transfer
│   ├── utilities/            # Pre/post tools
│   │   ├── mesh/            # Mesh manipulation
│   │   ├── preProcessing/   # Case setup
│   │   └── postProcessing/  # Analysis tools
│   └── test/                # Test applications
├── modules/                  # Optional extensions
│   ├── OpenQBMM/            # Quadrature-based moments
│   ├── avalanche/           # Granular flows
│   ├── adios/               # ADIOS I/O
│   └── visualization/       # Advanced visualization
├── plugins/                  # Third-party integrations
│   ├── paraview/            # ParaView reader
│   └── cfmesh/              # cfMesh integration
├── wmake/                    # Build system
│   ├── rules/               # Platform rules
│   ├── scripts/             # Build scripts
│   └── wmake                # Main build tool
├── tutorials/                # Example cases
└── META-INFO/                # Package metadata
```

### ESI-Specific Features

**Modules** (v2506 exclusive):
- OpenQBMM: Quadrature-based moment methods
- avalanche: Granular and particle-laden flows
- external-solver: Coupling with external solvers
- adios: ADIOS2 I/O for large-scale simulations
- visualization: Advanced post-processing

**Enhanced Capabilities**:
```cpp
// Adaptive mesh refinement
dynamicRefineFvMesh
{
    type            dynamicRefineFvMesh;
    refineInterval  1;
    field           alpha.water;
    lowerRefineLevel 0.001;
    upperRefineLevel 0.999;
}

// Overset mesh support
oversetMesh
{
    type            overset;
    active          true;
}
```

**Function Objects** (extended):
```cpp
functions
{
    #includeFunc Q        // Q-criterion
    #includeFunc vorticity
    #includeFunc yPlus
    #includeFunc MachNo
    #includeFunc CourantNo
    
    surfaces
    {
        type            surfaces;
        libs            (sampling);
        writeControl    writeTime;
        surfaces
        {
            plane1
            {
                type        plane;
                planeType   pointAndNormal;
                point       (0 0 0);
                normal      (0 0 1);
            }
        }
        fields          (U p);
    }
}
```

### Build System (wmake)

**Make/files**:
```make
mySource1.C
mySource2.C

LIB = $(FOAM_USER_LIBBIN)/libmyLibrary
```

**Make/options**:
```make
sinclude $(GENERAL_RULES)/module

EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/turbulenceModels/incompressible/lnInclude

LIB_LIBS = \
    -lfiniteVolume \
    -lmeshTools \
    -lincompressibleTurbulenceModels
```

## Development Workflow

### Creating New Solver
```bash
# Copy template
cd $WM_PROJECT_USER_DIR
cp -r $FOAM_SOLVERS/incompressible/simpleFoam mySolver
cd mySolver

# Modify Make/files
sed -i "s|FOAM_APPBIN|FOAM_USER_APPBIN|g" Make/files
sed -i "s|simpleFoam|mySolver|g" Make/files

# Rename main file
mv simpleFoam.C mySolver.C

# Build
wmake
```

### Adding Custom Boundary Condition
```bash
# Create BC in user directory
mkdir -p $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/myBC
cd $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/myBC

# Copy template and modify
cp -r $FOAM_SRC/finiteVolume/fields/fvPatchFields/derived/fixedValue/* .
# Edit files...

# Build
wmake libso
```

### Parallel Execution
```bash
# Decomposition setup
cat > system/decomposeParDict << EOF
numberOfSubdomains  4;
method             scotch;

scotchCoeffs
{
    processorWeights (1 1 1 1);
}
EOF

# Decompose
decomposePar

# Run parallel
mpirun -np 4 simpleFoam -parallel

# Monitor
tail -f log.simpleFoam

# Reconstruct
reconstructPar
reconstructPar -latestTime

# Clean
rm -rf processor*
```

### Debugging
```bash
# Debug build
export WM_COMPILE_OPTION=Debug
source etc/bashrc
./Allwmake -j

# Run with debugger
gdb simpleFoam
(gdb) run -case <case>

# Memory debugging
valgrind --leak-check=full simpleFoam -case <case>

# Enable floating point exceptions
export FOAM_SIGFPE=true
export FOAM_SETNAN=true

# Debug output in controlDict
DebugSwitches
{
    lduMatrix       2;
    GAMGSolver      1;
}

# Profiling
simpleFoam -case <case> -profiling
```

### Using Modules
```bash
# Build specific module
cd modules/OpenQBMM
./Allwmake

# Use in case
# Add to controlDict:
libs ("libOpenQBMM.so");

# Run module tutorials
cd tutorials
./Allrun
```

## Common Operations

### Environment Variables
```bash
# Core variables
$WM_PROJECT         # OpenFOAM
$WM_PROJECT_VERSION # v2506
$WM_PROJECT_DIR     # Installation root
$FOAM_SRC           # Source directory
$FOAM_APP           # Applications
$FOAM_TUTORIALS     # Tutorial cases
$FOAM_RUN           # User run directory

# User directories
$WM_PROJECT_USER_DIR # User's OpenFOAM directory
$FOAM_USER_APPBIN    # User applications
$FOAM_USER_LIBBIN    # User libraries

# Build settings
$WM_COMPILER        # Gcc/Clang/Icc
$WM_COMPILE_OPTION  # Debug/Opt
$WM_PRECISION_OPTION # SP/DP/LP
$WM_LABEL_SIZE      # 32/64
$WM_MPLIB           # MPI implementation
```

### Useful Utilities
```bash
# Mesh quality
checkMesh -allGeometry -allTopology

# Convert mesh formats
fluent3DMeshToFoam <mesh.msh>
ideasUnvToFoam <mesh.unv>
gmshToFoam <mesh.msh>

# Field operations
postProcess -func mag(U)
postProcess -func 'div(U)'

# Sampling
postProcess -func sampleDict

# Force calculation
postProcess -func forces

# Time averaging
postProcess -func fieldAverage
```

### Performance Optimization
```bash
# Solver performance
simpleFoam -case <case> -profiling

# Check decomposition quality
decomposePar -cellDist

# Optimize parallel communication
# In decomposeParDict:
preservePatches (wall inlet outlet);
preserveFaceZones (zone1 zone2);

# Cache optimization
export OMP_NUM_THREADS=1
export FOAM_SIGFPE=false  # Disable for production runs
```

## Important Notes

1. **Environment conflicts**: NEVER source multiple OpenFOAM versions in the same shell
2. **Module dependencies**: Build main OpenFOAM before modules/plugins
3. **Version compatibility**: Check module compatibility with v2506
4. **Parallel scalability**: Test decomposition strategies for large cases
5. **Memory management**: Use `valgrind` for memory leak detection
6. **Code standards**: Follow OpenFOAM coding style guide
7. **Testing protocol**: Validate with tutorials before production use
8. **Build optimization**: Use `-j$(nproc)` for parallel compilation
9. **Debug vs Release**: Use Debug builds only for development
10. **File cleanup**: Never commit `lnInclude/`, `*.dep`, or `platforms/` directories

## Cross-References

- **[Root CLAUDE.md](../CLAUDE.md)**: Multi-project overview and comparison
- **[OpenFOAM-13](../OpenFOAM-13/CLAUDE.md)**: Foundation version with modular solvers
- **[HELYXcore](../HELYXcore-4.4.0/CLAUDE.md)**: CMake-based unified solver framework
- **[Root AGENTS.md](../AGENTS.md)**: Development guidelines and commit conventions
- **[Project README](README.md)**: OpenFOAM-v2506 specific information
- **[Modules README](modules/README.md)**: Extended modules documentation
- **[Plugins README](plugins/README.md)**: Third-party plugin information