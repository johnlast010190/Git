# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains three OpenFOAM-based CFD projects:

1. **[OpenFOAM-13](OpenFOAM-13/CLAUDE.md)**: Foundation core version with modular solver framework
   - Focus: Modular solvers (`foamRun`/`foamMultiRun`)
   - Build: wmake system
   - [README](OpenFOAM-13/README.md) • [AGENTS](OpenFOAM-13/AGENTS.md) • [CLAUDE](OpenFOAM-13/CLAUDE.md)

2. **[OpenFOAM-v2506](OpenFOAM-v2506/CLAUDE.md)**: ESI-OpenCFD fork with extended capabilities
   - Focus: Extended modules, plugins, advanced features
   - Build: wmake + optional modules/plugins
   - [README](OpenFOAM-v2506/README.md) • [AGENTS](OpenFOAM-v2506/AGENTS.md) • [CLAUDE](OpenFOAM-v2506/CLAUDE.md)

3. **[HELYXcore-4.4.0](HELYXcore-4.4.0/CLAUDE.md)**: CMake-based unified solver framework
   - Focus: Unified multi-physics solver, GUI/tools
   - Build: CMake out-of-source
   - [README](HELYXcore-4.4.0/README.md) • [AGENTS](HELYXcore-4.4.0/AGENTS.md) • [CLAUDE](HELYXcore-4.4.0/CLAUDE.md)

## Quick Start Guide

### Environment Setup
```bash
# CRITICAL: Only source ONE OpenFOAM environment per shell
# Option 1: OpenFOAM-13
source OpenFOAM-13/etc/bashrc

# Option 2: OpenFOAM-v2506
source OpenFOAM-v2506/etc/bashrc

# Option 3: HELYXcore (after building)
source HELYXcore-4.4.0/build/bin/helyxcore-bashrc

# Verify environment
echo $WM_PROJECT_DIR    # OpenFOAM projects
echo $HELYXCORE_DIR     # HELYXcore
foamVersion            # Check version info
```

### Build Commands Summary

| Project | Initial Build | Clean | Rebuild |
|---------|--------------|-------|---------|
| **OpenFOAM-13** | `./Allwmake -j`<br>`./Allwmake -j` (2nd pass) | `wclean all` | `wclean all && ./Allwmake -j` |
| **OpenFOAM-v2506** | `./Allwmake -j -s -q -l`<br>`./Allwmake-modules`<br>`./Allwmake-plugins` | `wclean all` | `wclean all && ./Allwmake -j -s -q -l` |
| **HELYXcore** | `cmake -S . -B build`<br>`cmake --build build -j` | `rm -rf build` | `rm -rf build && cmake -S . -B build && cmake --build build -j` |

### Testing Quick Reference

| Test Type | OpenFOAM-13 | OpenFOAM-v2506 | HELYXcore |
|-----------|-------------|----------------|-----------|
| **Unit Tests** | `cd test && ./Allrun` | `cd test && ./Allrun` | N/A |
| **Tutorials** | `cd tutorials && ./Alltest` | `cd tutorials && ./Alltest` | `cd examples && ./Allrun` |
| **Quick Test** | `tutorials/incompressibleFluid/pitzDaily` | `tutorials/incompressible/simpleFoam/pitzDaily` | `examples/helyxSolve/incompressible/cavity` |

## Project Comparison Matrix

| Feature | OpenFOAM-13 | OpenFOAM-v2506 | HELYXcore-4.4.0 |
|---------|-------------|----------------|-----------------|
| **Version** | Foundation v13 | ESI-OpenCFD v2506 | 4.4.0 (OpenFOAM-based) |
| **Build System** | wmake | wmake | CMake |
| **Modular Solvers** | ✅ foamRun/foamMultiRun | ❌ Traditional | ✅ helyxSolve unified |
| **Extended Modules** | ❌ | ✅ OpenQBMM, avalanche, etc. | ✅ Integrated |
| **Plugins** | ❌ | ✅ cfmesh, turbulence-community | ❌ |
| **Overset Mesh** | ❌ | ✅ | ❌ |
| **ADIOS2 I/O** | ❌ | ✅ | ❌ |
| **GUI Tools** | ❌ | ❌ | ✅ helyxMesh, helyxReport |
| **Runtime Post-processing** | ✅ Function objects | ✅ Extended function objects | ✅ Integrated reporting |

## Code Architecture Overview

### Common Structure
```
<project>/
├── src/                    # Core libraries
│   ├── OpenFOAM/          # Base classes, containers, I/O
│   ├── finiteVolume/      # FV discretization
│   ├── meshTools/         # Mesh manipulation
│   └── ...                # Domain-specific libraries
├── applications/          # Executables
│   ├── solvers/          # CFD solvers
│   └── utilities/        # Pre/post tools
├── tutorials/            # Example cases (OpenFOAM)
├── examples/             # Example cases (HELYXcore)
└── test/                 # Unit tests
```

### Key Differences

**OpenFOAM-13 (Modular Focus)**:
- `applications/modules/`: Modular solver components
- `foamRun`: Runtime module loading
- Emphasis on code reusability

**OpenFOAM-v2506 (Extended Features)**:
- `modules/`: Optional advanced capabilities
- `plugins/`: Third-party integrations
- Extended boundary conditions and models

**HELYXcore (Unified Framework)**:
- `helyxSolve`: Single multi-physics solver
- `helyxMesh`: Advanced meshing with GUI
- CMake-based for easier integration

## Development Workflow

### 1. Environment Management
```bash
# Create separate terminals for different projects
# Terminal 1: OpenFOAM-13
source OpenFOAM-13/etc/bashrc

# Terminal 2: OpenFOAM-v2506
source OpenFOAM-v2506/etc/bashrc

# Terminal 3: HELYXcore
source HELYXcore-4.4.0/build/bin/helyxcore-bashrc
```

### 2. Code Navigation
```bash
# Find implementations
grep -r "className" src/
find . -name "*.H" -exec grep -l "className" {} \;

# Locate tutorials/examples
find tutorials -name "controlDict" -exec grep -l "solver.*<name>" {} \;
find examples -type d -name "<case-name>"

# Check dependencies
wmake -list-dependencies  # OpenFOAM
cmake --build build --target help  # HELYXcore
```

### 3. Building Components

**OpenFOAM Projects**:
```bash
cd src/<library>
wmake libso        # Build library
wmake -a          # Force rebuild

cd applications/<solver>
wmake             # Build application
```

**HELYXcore**:
```bash
cmake --build build --target <target>
cmake --build build -j$(nproc)  # Parallel build
```

### 4. Testing Changes
```bash
# Quick validation
foamSystemCheck
foamInstallationTest

# Run specific test case
cd <test-case>
./Allrun
./Allclean  # Clean after test

# Compare results
diff -r <results> <reference>
```

### 5. Debugging
```bash
# Debug builds
export WM_COMPILE_OPTION=Debug  # OpenFOAM
cmake -DCMAKE_BUILD_TYPE=Debug  # HELYXcore

# Runtime debugging
export FOAM_SIGFPE=true        # Floating point exceptions
export FOAM_SETNAN=true        # Initialize with NaN
gdb <solver>                   # Debug with gdb
valgrind --leak-check=full <solver>  # Memory checks
```

## Cross-Project Guidelines

### File Naming Conventions
- Headers: `.H` (all projects)
- Sources: `.C` (all projects)
- Dictionaries: `Dict` suffix
- Patches: `.patch` files
- CMake: `CMakeLists.txt` (HELYXcore only)

### Code Style (Universal)
```cpp
namespace Foam
{

class myClass
:
    public baseClass
{
    // Private data
        scalar data_;

public:

    // Constructors
        myClass();

    // Member Functions
        scalar calculate() const;
};

} // End namespace Foam
```

### Testing Protocol
1. **Unit tests**: For core library changes
2. **Tutorial validation**: For solver modifications
3. **Parallel testing**: For decomposition-sensitive changes
4. **Performance profiling**: For optimization work

## Important Notes & Best Practices

1. **Environment Isolation**: NEVER source multiple OpenFOAM environments in the same shell
2. **Build Order**: 
   - OpenFOAM: May need 2 passes for circular dependencies
   - HELYXcore: Always use out-of-source builds
3. **Version Control**:
   - Never commit: `platforms/`, `lnInclude/`, `*.dep`, `build/`
   - Use `.gitignore` appropriately
4. **Path Variables**:
   - Use `$WM_PROJECT_DIR`, `$FOAM_*` for OpenFOAM
   - Use `$HELYXCORE_DIR` for HELYXcore
5. **Parallel Builds**: Always use `-j` or `-j$(nproc)` for faster compilation
6. **Documentation**:
   - Update CLAUDE.md for AI assistance guidance
   - Follow AGENTS.md for development guidelines
   - Include README for new features/modules
7. **Testing Coverage**:
   - Minimum: Run relevant tutorial case
   - Preferred: Add test case demonstrating feature
   - Required: Ensure no regression in existing tests
8. **Cross-Compatibility**:
   - Check if feature exists in other projects
   - Consider portability when adding new code
   - Document project-specific features clearly

## Quick Reference Links

### Documentation
- [Root Guidelines](AGENTS.md)
- [OpenFOAM-13 Details](OpenFOAM-13/CLAUDE.md)
- [OpenFOAM-v2506 Details](OpenFOAM-v2506/CLAUDE.md)
- [HELYXcore Details](HELYXcore-4.4.0/CLAUDE.md)

### Key Commands
- Build: `./Allwmake -j` (OpenFOAM) | `cmake --build build -j` (HELYXcore)
- Clean: `wclean all` (OpenFOAM) | `rm -rf build` (HELYXcore)
- Test: `./Alltest` (tutorials) | `./Allrun` (examples)
- Debug: `export WM_COMPILE_OPTION=Debug` | `-DCMAKE_BUILD_TYPE=Debug`

### Environment Variables
- `$WM_PROJECT_DIR`: OpenFOAM root directory
- `$FOAM_SRC`: Source code location
- `$FOAM_APP`: Applications directory
- `$FOAM_TUTORIALS`: Tutorial cases
- `$HELYXCORE_DIR`: HELYXcore root
- `$FOAM_USER_APPBIN`: User applications
- `$FOAM_USER_LIBBIN`: User libraries