# Unidirectional Close-Contact Melting with Overset Remeshing

## 1. Problem Description
This simulation models the unconstrained, close-contact melting of eicosane within an isothermal square cavity heated from all sides. The focus is on resolving the continuous hydrodynamic pressure within the thin liquid film beneath the solid, which naturally balances gravitational forces without artificially relying on a contact model.

For comprehensive details regarding all problem configurations and their underlying physics, the reader is referred to the work of Van Riet et al. [[1](#1)].

**Geometry and Boundary Conditions:**
* Domain: 52 mm by 52 mm square cavity.
* Walls: No-slip boundaries maintained at a constant 329.85 K (a 20 K superheat above the melting temperature).
* Solid Bulk: Initialised at the melting temperature, Tm = 309.85 K, with zero initial velocity.
* Interface: The solid-liquid interface acts as a moving isothermal wall at Tm with a no-slip condition. 

**Material Properties (Eicosane):**
Volume change is neglected by setting the solid density equal to the liquid density. Natural convection is modelled using the Boussinesq approximation ($\beta = 8.61 \cdot 10^{-4}$ K$^{-1}$).

| Phase | Dynamic Viscosity | Specific Heat | Thermal Conductivity | Density | Latent Heat |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Solid** | N/A | 1926 J/(kg K) | 0.4 W/(m K) | 856 kg/m³ | 248 kJ/kg |
| **Liquid** | 0.00401 kg/(m s) | 2400 J/(kg K) | 0.15 W/(m K) | 775.13 kg/m³ | 248 kJ/kg |

## 2. Solver and Coupling Settings
The problem is solved using a partitioned fluid-structure interaction approach, coupling ANSYS Fluent (liquid/flow) with a Python solver (solid/kinematics).

**Fluid Solver (ANSYS Fluent 2024R2):**
* Formulation: Pressure-based, coupled pressure-velocity solver.
* Discretisation: PRESTO! scheme for pressure, second-order upwind for momentum and energy, second-order implicit for time.
* Flow Regime: Laminar flow (maximum Rayleigh number is $Ra \approx 5.7 \cdot 10^7$).

**Rigid Body & Coupling Dynamics:**
* Rigid Body Solver: 3-DoF solver utilising a linear predictor.
* Stability: A fictitious mass scaling factor of $\psi = 8.0$ is applied. Mechanical sub-iterations are capped at 10 per time step.
* Convergence Tolerances: Absolute tolerances are set to $10^{-3}$ N for forces and $10^{-3}$ Nm for moments. 
* Contact Model: Acts strictly as a numerical safeguard against divergence during asymmetric melting. It is not intended to support the solid bulk under normal conditions.

## 3. Mesh Topology and Initialisation
To handle the extensive rigid body translation and local melting deformations, an overset mesh topology is utilised alongside a two-step initialisation process.

**Overset Mesh Topology:**
* Background Mesh: Static 832 by 837 grid. Uniform 0.0625 mm square cells in the bulk, with a refined boundary layer at the bottom wall (10 layers, 0.015 mm initial thickness, 1.16 growth rate) to capture intense lubrication gap gradients.
* Component Mesh: Dynamic, conformal boundary layer attached to the interface (9 layers, 1.16 growth rate, 0.263 mm total thickness). Radial-basis function smoothing accommodates local shape deformation.
* Overset Interpolation: Least-squares interpolation and a cell volume-based donor-priority method. The extremely thin component cells spanning the close-contact gap are dynamically disabled as dead cells.

**Two-Step Initialisation:**
1. Shape Extraction: The initial interface morphology (liquid fraction of 0.5) is extracted from a preliminary enthalpy-porosity simulation at a flow time of 2.3 s (yielding a 0.25 mm bottom film).
2. Field Initialisation: Hydrodynamic and thermal fields are initialised via a steady-state conduction calculation, followed by a 1 s transient flow simulation with solid kinematics locked. This ensures a fully developed pressure field before activating the rigid body solver.

## 4. The Automated Remeshing Workflow
In a standard CoCoNuT simulation, the solver initialises once and runs continuously until the final time step is reached. However, due to the extreme deformations in close-contact melting, the conformal component mesh degrades. To maintain cell quality, the workflow must periodically halt, generate a completely new mesh, and interpolate the solution.

To achieve this, the standard execution is replaced by a state-driven batch architecture managed by an `orchestrator.py` script.

**How the Auto-Remesh Script Works:**
* Centralised Configuration: A single `config.yaml` file controls cluster names, directory paths, sync targets, and parameter toggles. Users no longer need to edit Python or JSON files directly.
* State Tracking: A `state.json` file tracks the current run number. If the simulation stops, it automatically resumes from the last successful batch.
* Dynamic Parameter Injection: The orchestrator dynamically rewrites `parameters.json` before execution, injecting initial or restart configurations (like `case_restart.cas.h5`) based on the current loop phase.
* The Loop Cycle: Every 100 time steps, the orchestrator triggers the following remote pipeline:
    1. Collect Data: Consolidates and offsets `report-file.out` and rigid body data.
    2. Remesh: Archives the old folders, triggers `remesh.sh` (Gambit/Python background meshing), and runs Fluent journal files to interpolate the fields onto the new grid.
    3. Update Interface: Translates the new Fluent boundary nodes into a fresh `interface_nodes.csv` for the solid Python solver.
    4. Restart: Launches CoCoNuT with the newly configured parameters.

## 5. Manual Crash Recovery
If the automated remeshing fails or generates an overly distorted mesh, the simulation will crash (or generate unphysical results). The orchestrator is designed to handle these interruptions gracefully.

**The Recovery Process:**
1. Catching the Interrupt: If the solver crashes remotely, or if the user spots a drop in accuracy and presses `Ctrl+C` locally, the orchestrator catches the exception.
2. Safe Pause: It immediately logs into the cluster, kills any orphaned Fluent/API processes to save memory and licenses, and pauses the script with a prompt.
3. Manual Injection: The user selects to pause and wait for a manual remesh. They generate a high-quality mesh in Workbench, save it as a `.msh` file in the `manual_interventions` folder, and press ENTER.
4. Resume: The orchestrator bypasses the automated Gambit scripting, injects the user's manual mesh directly into the workflow, interpolates the data, and resumes the coupled simulation loop perfectly from where it left off.

---

## References

<a id="1">[1]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “A partitioned interface-tracking method for close-contact melting of unconstrained saturated solids,” submitted to Int. J. Heat Mass Transfer., 2026.]
