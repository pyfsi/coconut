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

---

## 4. The Automated Remeshing Workflow
Because the component mesh slowly degrades due to the unconstrained melting deformations, the simulation cannot run continuously. It must periodically halt, generate a new mesh, and interpolate the flow field. This is managed by a state-driven Python orchestration framework.

**How the Auto-Remesh Script Works:**
* **Centralised Configuration:** A single `config.yaml` file controls cluster environments (including HPC toggles), directory paths, and parameter settings. 
* **Dynamic Parameter Injection:** The orchestrator reads the base `parameters.json` file and dynamically injects the correct case files (`pc_fluid.cas.h5` vs `case_restart.cas.h5`) and restart values depending on whether it is initiating a new run or looping a restart.
* **The Hybrid Loop Cycle:** Every 100 time steps, the orchestrator triggers the pipeline:
    1. **Sync & Retrieve:** Automatically retrieves the latest rigid-body and Fluent report data from the HPC to the local processing node.
    2. **Collect Data:** Consolidates and offsets the time series data for continuous post-processing.
    3. **Remesh:** Archives the old output folders, executes `remesh.sh` (Gambit/Python background meshing), and uses Fluent journal files to interpolate the previous fields onto the updated grid.
    4. **Update Interface:** Translates the newly generated Fluent boundary nodes into a fresh `interface_nodes.csv` for the solid Python solver.
    5. **HPC Submission:** Syncs the fresh mesh back to the HPC, submits the restart PBS job, and monitors the queue.

---

## 5. Manual Crash Recovery
If the automated remeshing generates an overly distorted mesh, or if the Fluent solver diverges due to complex asymmetric contact, the orchestrator is designed to handle the interruption gracefully without breaking the workflow loop.

**The Recovery Process:**
1. **Catching the Interrupt:** If the HPC job crashes, or if the user spots a drop in accuracy and presses `Ctrl+C` locally, the orchestrator catches the exception.
2. **Safe Pause:** It immediately logs into the remote cluster, cancels any lingering PBS jobs (`qdel`), kills any orphaned Fluent/API processes to save memory and licensing tokens, and pauses the local script with a prompt.
3. **Manual Injection:** The user selects to pause and wait for a manual remesh. The user can then generate a high-quality, stable mesh manually in Workbench, save it as a `.msh` file in the `manual_interventions` folder, and press ENTER.
4. **Resume:** The orchestrator bypasses the automated Gambit scripting, injects the user's manual mesh directly into the workflow, interpolates the data, and resumes the coupled simulation loop seamlessly.

## 6. Starting the simulation
Navigate to the folder containing the `orchestrator.py` script on the login node and run it. The script will manage the cluster assignment for each step in the remeshing procedure.

---

## References

<a id="1">[1]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “A partitioned interface-tracking method for close-contact melting of unconstrained saturated solids,” submitted to Int. J. Heat Mass Transfer., 2026.]
