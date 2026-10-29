# Asymmetric Close-Contact Melting with Overset Remeshing

## 1. Problem Description
This simulation investigates the unconstrained close-contact melting (CCM) of eicosane within a complex, asymmetric geometry. The domain represents a two-dimensional cross-section of the upper third of a finned triplex tube heat exchanger, featuring an outer circular wall (51.8 mm radius, $120^\circ$ arc) and two straight internal heated walls.

Due to the geometric asymmetry relative to the solid's centre of mass, the sinking bulk experiences non-uniform hydrodynamic pressure and heat fluxes. This induces complex 3-DoF kinematics, requiring the solver to calculate simultaneous rigid body translation and rotation.

For comprehensive details regarding all problem configurations and their underlying physics, the reader is referred to the work of Van Riet et al. [[1](#1)].

**Geometric Orientation ($0^\circ$ 'Mercedes' configuration):**
The domain features one vertical inner wall and one downward-inclined inner wall. Because the simulation explicitly models variable liquid density (and therefore volumetric expansion during phase change), an additional rectangular expansion volume (10 mm wide and 30 mm high) is included at the top boundary.

**Boundary Conditions:**
* **Inner Radial Walls:** Maintained at a constant heat flux of 8000 W/m² (10 W per wall).
* **Outer Wall & Expansion Volume:** Maintained at a constant temperature of 311.35 K.
* **Solid-Liquid Interface:** Isothermal at the melting temperature (309.85 K).
* **Top Boundary:** Zero gauge pressure outlet to accommodate the expanding air phase.
* **Flow:** Standard no-slip conditions on all solid boundaries and the interface.

**Material Properties (Eicosane):**
The simulation explicitly accounts for volumetric expansion during melting by employing temperature-dependent liquid density parameters alongside a primary air phase.

---

## 2. Solver and Coupling Settings
The problem is solved using a partitioned fluid-structure interaction approach, coupling ANSYS Fluent (liquid/flow) with a Python solver (solid/kinematics) via explicit coupling.

**Fluid Solver (ANSYS Fluent 2024R2):**
* **Formulation:** Pressure-based, coupled pressure-velocity solver.
* **Flow Regime:** Laminar flow (maximum Rayleigh number $Ra \approx 1.59 \times 10^8$, safely below the transition threshold).
* **Discretisation:** Second-order spatial and temporal discretisation. 
* **Multiphase:** Utilises an implicit Volume of Fluid (VOF) model with an implicit body force formulation. Pressure interpolation uses the Body Force Weighted scheme.

**Rigid Body & Coupling Dynamics:**
* **Rigid Body Solver:** 3-DoF solver using a linear predictor. Mechanical sub-iterations are capped at 10.
* **Stability:** A fictitious mass scaling factor of $\psi = 8.0$ is applied initially to maintain stability, but is progressively lowered to 0.5 during the rapid melting phase to reduce high-frequency kinematic noise without sacrificing stability.
* **Convergence Tolerances:** Set to $10^{-3}$ for both relative norm (forces/moments) and absolute limits ($10^{-3}$ N and $10^{-3}$ Nm).
* **Contact Model:** Due to geometric asymmetry, localised contact is inevitable. A repulsive contact force acts as a strict numerical safeguard against solver divergence. Activation occurs at an upper gap limit of 0.20 mm and a lower limit of 0.10 mm, ensuring contact patches occupy no more than 1% to 5% of the solid's circumference.

---

## 3. Mesh Topology and Initialisation
To resolve the 3-DoF kinematics without severe cell distortion, the simulation relies on a finely tuned overset mesh topology and a dynamic initialisation strategy.

**Overset Mesh Topology:**
* **Background Mesh:** An unstructured grid of 244,749 cells with a standard 0.12 mm bulk size. The heated walls are lined with 18 inflation layers (0.02 mm first-layer thickness, 1.11 growth rate) to capture intense thermal gradients in the lubrication gap. The top expansion volume is discretised using a structured mesh, locally refined near its intersection with the circular domain.
* **Component Mesh:** A conformal boundary layer enveloping the solid interface, featuring 0.08 mm circumferential spacing and 13 layers (0.35 mm total thickness). Radial-basis function smoothing absorbs local melting deformations.
* **Critical Sizing Constraint:** The 0.35 mm component mesh thickness ensures its outermost cells perfectly match the 0.12 mm resolution of the background grid. This geometric matching is absolutely essential for VOF phase conservation, preventing artificial numerical loss of the phase change material into the air phase as the solid sinks.

**Two-Step Initialisation:**
1. **Shape Extraction:** Rather than prescribing an analytical shape, a preliminary enthalpy-porosity simulation is run on a fixed grid until a flow time of 18 s. This produces a natural 0.5 mm liquid gap. The 0.5 liquid fraction iso-surface is extracted to define the starting overset component boundary.
2. **Field Initialisation:** Hydrodynamic and thermal fields are initialised via a steady-state conduction calculation, followed by a 1 s transient flow simulation with the solid kinematics locked. This establishes the necessary physical pressure field beneath the solid before the rigid body solver takes over.

## 4. Automated Remeshing Workflow and Manual Crash Recovery

In a standard CoCoNuT simulation, the solver initialises once and runs continuously until the final time step is reached. However, due to the extreme deformations in close-contact melting, the conformal component mesh degrades. To maintain cell quality, the workflow must periodically halt, generate a completely new mesh, and interpolate the solution.

To achieve this, the standard execution is replaced by a state-driven batch architecture managed by an `orchestrator.py` script. The core principles behind this framework include:

* **Chunked Batch Execution:** Rather than running as a single uninterrupted transient calculation, the simulation is partitioned into discrete, manageable batches. Each batch advances physical time by a specified interval before systematically pausing the solvers.
* **Spatial Domain Renewal:** At the end of each batch, the deforming fluid domain is stripped and reconstructed with a fresh mesh topology to restore element quality, boundary layer resolution, and skewness. The physical fields (pressure, velocity, and temperature) are interpolated from the previous state to maintain temporal and physical continuity across restarts.
* **Autonomous State Tracking:** A dedicated state tracker monitors the simulation lifecycle. It handles run increments, accumulates time-history outputs into continuous post-processing logs, and dynamically configures solver parameters between initial and restart modes.
* **Human-in-the-Loop Fault Tolerance:** When extreme local topological shifts cause automated meshing or solver iterations to diverge, the orchestrator intercepts the failure and cleans active background processes. It automatically exports the current interface geometry for external CAD meshing, allowing a tailored, manually constructed mesh to be injected directly so the calculation can resume without loss of progress.

*(For full technical specifications, script dependencies, HPC communication protocols, and step-by-step manual recovery instructions, please refer to the dedicated [Remeshing Workflow Documentation](../downwards_ccm/remeshing_workflow.md).)*

## 5. Starting the simulation
Navigate to the folder containing the `orchestrator.py` script on the login node and run it. The script will manage the cluster assignment for each step in the remeshing procedure.

---

## References

<a id="1">[1]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “A partitioned interface-tracking method for close-contact melting of unconstrained saturated solids,” submitted to Int. J. Heat Mass Transfer., 2026.]
