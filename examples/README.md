# MoM3D Examples

This directory contains well-documented examples demonstrating practical applications of the 3D Method of Moments (MoM) implementation, based on the theoretical framework in Gibson's "The Method of Moments in Electromagnetics."

## Overview

The Method of Moments is a powerful numerical technique for solving electromagnetic problems involving conducting and dielectric objects. These examples showcase the most common applications of MoM in electromagnetic engineering.

## Examples

### 1. Conducting Sphere Scattering (`conducting_sphere_scattering.jl`)

**Theoretical Background:** Gibson Sections 6.6.2 and 7.7.2

Demonstrates electromagnetic scattering from a perfectly conducting sphere, which serves as a canonical validation problem due to the availability of analytical solutions (Mie theory).

**Key Features:**
- Spherical mesh generation with quality analysis
- Mesh quality assessment per Gibson Section 7.6.1 requirements
- Electrical size analysis for different frequencies
- Expected accuracy: < 0.1 dB compared to Mie theory

**Applications:**
- Radar cross section prediction
- Electromagnetic compatibility studies
- Numerical method validation
- Scattering analysis

**Usage:**
```julia
julia conducting_sphere_scattering.jl
```

### 2. Antenna Impedance Calculation (`antenna_impedance_calculation.jl`)

**Theoretical Background:** Gibson antenna analysis throughout the text

Shows the MoM approach for calculating antenna input impedance, a fundamental parameter in antenna design and analysis.

**Key Features:**
- Wire dipole and patch antenna mesh generation
- Impedance calculation workflow
- Resonance frequency estimation
- Mesh density requirements for antenna analysis

**Applications:**
- Antenna design and optimization
- Impedance matching network design
- Bandwidth and efficiency analysis
- Electromagnetic compatibility studies

**Expected Results:**
- Wire dipole (λ/2): ~73 Ω resistance at resonance
- Rectangular patch: ~100-300 Ω (depends on feed location)

**Usage:**
```julia
julia antenna_impedance_calculation.jl
```

### 3. Radar Cross Section Prediction (`radar_cross_section_prediction.jl`)

**Theoretical Background:** Gibson Chapter 7, Sections 7.7.2 and 7.7.3

Demonstrates RCS calculation for various canonical targets, including validation against known analytical solutions and EMCC benchmark targets.

**Key Features:**
- Multiple target geometries (sphere, plate, corner reflector)
- Electrical size analysis
- RCS characteristics estimation
- Validation approaches per Gibson methodology

**Applications:**
- Stealth technology development
- Radar system design and analysis
- Target identification and classification
- Electromagnetic signature analysis

**Usage:**
```julia
julia radar_cross_section_prediction.jl
```

## Theoretical Foundation

All examples are based on the theoretical framework presented in Gibson's "The Method of Moments in Electromagnetics," with specific emphasis on:

### Mesh Quality Requirements (Gibson Section 7.6.1)
- Triangle aspect ratios should be reasonable (< 5.0 preferred)
- Minimum angles should be > 10° to avoid poorly conditioned matrices
- Mesh density typically λ/10 to λ/20 for good accuracy

### Watertight Meshes (Gibson Section 7.6.2)
- Critical for MFIE formulation
- No T-junctions or boundary edges for closed surfaces
- Proper current continuity across triangle edges

### Validation Approaches
- Conducting sphere: Compare with Mie theory (Gibson Section 7.7.2)
- EMCC benchmark targets: Measured vs. computed (Gibson Section 7.7.3)
- Convergence analysis with mesh refinement
- Reciprocity and energy conservation checks

## Common Workflow

All examples follow a similar MoM workflow:

1. **Geometry and Mesh Generation**
   - Create surface triangulation
   - Verify mesh quality
   - Analyze electrical size

2. **MoM Setup**
   - Define RWG basis functions
   - Assemble impedance matrix (EFIE/MFIE)
   - Apply boundary conditions

3. **System Solution**
   - Solve linear system for surface currents
   - Extract desired quantities (impedance, RCS, etc.)

4. **Post-Processing**
   - Visualize results
   - Compare with analytical solutions
   - Validate accuracy

## Dependencies

The examples require the following Julia packages:
- `MoM3D` (this package)
- `LinearAlgebra`
- `StaticArrays`
- `Printf`

## Running the Examples

Each example can be run independently:

```bash
cd examples
julia conducting_sphere_scattering.jl
julia antenna_impedance_calculation.jl
julia radar_cross_section_prediction.jl
```

Or from within Julia:

```julia
using Pkg
Pkg.activate(".")
include("examples/conducting_sphere_scattering.jl")
```

## Expected Output

Each example provides:
- Detailed mesh quality analysis
- Theoretical background and expected results
- Validation recommendations
- Next steps for complete MoM implementation

## References

1. Gibson, W.C., "The Method of Moments in Electromagnetics," 2nd Edition, Chapman & Hall/CRC, 2021.
2. Harrington, R.F., "Field Computation by Moment Methods," IEEE Press, 1993.
3. Rao, S.M., Wilton, D.R., and Glisson, A.W., "Electromagnetic scattering by surfaces of arbitrary shape," IEEE Trans. Antennas Propag., vol. 30, no. 3, pp. 409-418, 1982.

## Contributing

When adding new examples, please follow the established format:
- Include detailed theoretical background with Gibson references
- Provide comprehensive documentation
- Add mesh quality analysis
- Include validation approaches
- Follow the common workflow structure
