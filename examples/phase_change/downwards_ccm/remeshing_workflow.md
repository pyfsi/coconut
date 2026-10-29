# External remeshing workflow
The automated remeshing framework manages multi-run, coupled fluid-structure interaction (FSI) and 
close-contact melting (CCM) simulations within the CoCoNuT coupling environment. 
The workflow couples the liquid phase of a phase change material (PCM) in ANSYS Fluent with the solid phase, which is 
governed by rigid-body kinematics and phase change. During longer simulations, 
continuous boundary deformation due to melting distorts the computational cells of the liquid component mesh, 
necessitating regular geometry extraction, remeshing, field interpolation, 
and solver restart cycles without loss of temporal or physical continuity.
A more detailed explanation about component mesh degradation can be found in Van Riet et al. [[1](#1), [2](#2)].

The general workflow looks like this:

```
+----------------------------------------------------------------------------------------------------+
|                                           ORCHESTRATOR                                             |
|                                         (orchestrator.py)                                          |
+-------------------------------------------------+--------------------------------------------------+
                                                  |
           +--------------------------------------+--------------------------------------+
           |                                                                             |
           v                                                                             v
+--------------------------------------+                       +--------------------------------------+
|        AUTOMATED REMESHING           |                       |           MANUAL REMESHING           |
|            (Phase B Loop)            |                       |        (Crash Recovery Mode)         |
+--------------------------------------+                       +--------------------------------------+
| 1. collect_data.py                   |                       | 1. Solver divergence / cell crush    |
|    Appends continuous histories      |                       | 2. handle_crash() cleans processes   |
| 2. remesh_pipeline.py (Automated)    |                       | 3. export_designmodeler.py           |
|    - Archives completed run folder   |                       |    Adaptive RDP node decimation      |
|    - Runs Fluent step_1.jou          |                       | 4. User crafts new_mesh.msh in DM    |
|    - Obtains GAMBIT license lock     |                       | 5. remesh_pipeline.py (--manual-mesh)|
|    - Parses topology & runs GAMBIT   |                       |    - Pulls clean archive restart     |
|    - Runs Fluent step_2.jou          |                       |    - Injects manual mesh             |
| 3. solid_interface.py                |                       |    - Runs Fluent step_2.jou          |
|    Extracts boundary coordinates     |                       | 6. solid_interface.py                |
| 4. prepare_parameters(restart=True)  |                       | 7. prepare_parameters(restart=True)  |
| 5. run_simulation.py (CoCoNuT FSI)   |                       | 8. run_simulation.py (Retry Run i)   |
+--------------------------------------+                       +--------------------------------------+
                                                  |
                                                  v
                               +--------------------------------------+
                               |          PHASE C: FINALISATION       |
                               |          (Run N == total_runs)       |
                               +--------------------------------------+
                               | 1. collect_data.py (Final logs)      |
                               | 2. archive_run.py (Final archive)    |
                               +--------------------------------------+

```

---

## 1. System Architecture and Distributed Execution

## 1. System Architecture and Distributed Execution

The workflow coordinates tasks across distinct physical nodes or clusters to separate CFD calculations, geometry generation, and post-processing tasks, depending on machine performance and software licence availability. While initialisation and dedicated meshing operations typically run across local cluster nodes, the core coupled fluid-structure interaction calculations can run locally or be offloaded directly to an HPC cluster.

```
+----------------------------------------------------------------------------------------------------+
|                                      DISTRIBUTED NODE TOPOLOGY                                     |
+----------------------------------+----------------------------------+------------------------------+
|       local_setup (cfdclu32)     |      local_process (cfdclu32/33) |        remesh (cfdclu13)     |
+----------------------------------+----------------------------------+------------------------------+
| - Initial geometry build         | - CoCoNuT coupled execution      | - GAMBIT 2.4.6 remeshing     |
| - setup_case.py execution        | - Time-history collection        | - Fluent step_1 stripping    |
| - Fluent baseline initialisation | - Interface node extraction      | - Fluent step_2 injection    |
| - Local process sanitisation     | - DesignModeler curve export     | - Licence lock management    |
+----------------------------------+----------------------------------+------------------------------+
| (Alternatively: CoCoNuT and Fluent solve stages offloaded to HPC via automated PBS submission)      |
+----------------------------------------------------------------------------------------------------+

```

The node labels listed in the topology above (`cfdclu32`, `cfdclu33`, and `cfdclu13`) serve as an illustrative example of a hardware mapping defined in `config.yaml`. These roles are fully configurable and can be reassigned to different physical machines, virtual environments, or alternative compute clusters depending on resource requirements and network mounts.

When HPC mode is enabled in the configuration (`hpc.enabled: true`), the solver execution stage on `local_process` is automatically substituted. The orchestrator synchronises the active case folders to the supercomputer using `rsync`, submits the coupled solver script to the batch scheduler using PBS (`qsub`), polls job status until completion, and synchronises critical result files back to the local filesystem before advancing to the remeshing phase.

### Configuration Model (`config.yaml`)

Central execution settings reside in `config.yaml`, specifying machine hostnames, working directories, High Performance Computing (HPC) switches, and file synchronisation manifests:

* **`machines`**: Defines the target hostnames for individual tasks. `local_setup` handles case initialization, `local_process` runs solver execution and data extraction, and `remesh` hosts legacy meshing utilities (GAMBIT).
* **`coconut`**: Details parameter files and specifies the exact case file and rigid-body restart flags (`restart_rb`) for initial runs versus restart runs.
* **`directories`**: Establishes absolute filesystem paths on shared network storage (`/cfdfile1`), identifying subdirectories for the solid solver (`CFD_1`), fluid solver (`CFD_2`), template files (`setup_files`), and manual meshes (`manual_interventions`).
* **`simulation`**: Sets `total_runs`, defining the macro-loop execution horizon.

### Remote Execution Layer (`src/utils.py`)

All execution across machines is abstracted through helper functions in `src/utils.py`:

* **`execute_local_node(script_name, cluster_name, base_dir, extra_args)`**: Dispatches non-interactive Python commands to target cluster nodes via OpenSSH using pseudo-terminal allocation (`ssh -tt`). The environment loads Anaconda Python via login shells before running scripts with unbuffered standard output (`python3 -u`).
* **`clean_local_fluent(cluster_name, base_dir, cfd2_dir)`**: Connects to the simulation execution node to evaluate whether orphaned Fluent processes or license communication daemons (`apip-standalone`) remain active. It triggers any generated `cleanup-fluent*` scripts and force-kills lingering API wrappers.
* **`prepare_parameters(config, is_restart)`**: Reads `parameters.json`, navigates through the coupled solver wrapper array to locate the liquid solver wrapper (`pc_fluent_liquid`), dynamically injects the appropriate case file name (`pc_fluid.cas.h5` versus `case_restart.cas.h5`), updates the rigid-body restart flag, and writes the updated JSON back to disk.
* **HPC Primitives**: Houses functions (`rsync_to_hpc`, `rsync_from_hpc`, `submit_hpc_job`, `monitor_job_completion`, `cancel_hpc_job`) to handle PBS job submission (`qsub`), remote polling (`qstat`), status tracking, and file retrieval when running on external supercomputing clusters.

### High-Performance Computing (HPC) Subsystem and Remote Offloading

When HPC execution is activated (`hpc.enabled: true` in `config.yaml`), the orchestrator shifts solver execution from the local processing node to the Tier-2 supercomputing infrastructure (such as the VSC clusters at Ghent University). This remote execution layer is driven by non-interactive SSH channels, PBS batch scheduling, and robust network synchronisation routines implemented in `src/utils.py`.

```
+----------------------------------------------------------------------------------------------------+
|                                      HPC COMMUNICATION PIPELINE                                    |
+----------------------------------------------------------------------------------------------------+
|                                                                                                    |
|  Local Controller                          SSH / rsync Control Channel               HPC Cluster   |
|  (cfdfile1 storage)                                                             (VSC Scratch Tier) |
|                                                                                                    |
|  1. rsync_to_hpc()          ─── Upload CFD_1 & CFD_2 folders (rsync) ────────>  Scratch Working    |
|                                                                                 Directory          |
|  2. submit_hpc_job()        ─── qsub job_shinx.pbs via SSH ──────────────────>  PBS Batch Queue    |
|                                                                                                    |
|  3. monitor_job_completion() <── Poll qstat status & verify critical files ───  Lustre File Check  |
|                                                                                                    |
|  4. rsync_from_hpc()        <── Download updated CFD folders & results (rsync)  Completed Run Data |
|                                                                                                    |
|  5. delete_files_on_hpc()   ─── rm -f critical status files via SSH ─────────>  Scratch Cleaned    |
|                                                                                                    |
+----------------------------------------------------------------------------------------------------+

```

#### Communication Primitives and Protocol

All remote interaction with the HPC facility is mediated through the helper routine `run_hpc_command(command, config)`.

* **Authentication and Environment Setup**: The controller opens an SSH session to `hpc.host` (for example, `login.hpc.ugent.be`) using the designated VSC account credentials (`hpc.user`). Each remote invocation executes inside a login shell that explicitly swaps into the required cluster software architecture (such as `module swap cluster/shinx`) before changing into the remote scratch base directory (`hpc.base_dir`).

* **Timeout Protection**: Every remote SSH interaction is wrapped with a 30-second timeout to prevent network latency spikes or login-node unresponsiveness from hanging the local workflow orchestrator indefinitely.

---

#### The Five-Stage Remote Execution Cycle

Whenever a coupled simulation run commences under HPC mode, the orchestrator replaces `run_simulation.py` with five discrete remote execution steps:

##### 1. Pre-Execution Synchronisation (`rsync_to_hpc`)

Before job submission, the local controller uploads the active working folders defined in `config.yaml` (`CFD_1` and `CFD_2`) to the remote scratch filesystem:

```bash
rsync -avz --delete /cfdfile1/.../CFD_2 vscXXXXX@login.hpc.ugent.be:/scratch/gent/.../

```

The `--delete` flag guarantees that remote working directories match the local state, removing old residual files from earlier attempts while preserving network bandwidth through rsync delta-compression.

##### 2. Job Submission and Parsing (`submit_hpc_job`)

The controller submits the unified PBS submission script (`hpc.script`) via SSH using `qsub`:

```bash
ssh vscXXXXX@login.hpc.ugent.be "bash -l -c 'module swap cluster/shinx && cd ... && qsub job_shinx.pbs'"

```

The command output is parsed via regular expressions (`re.search(r"(\d+)", output)`) to extract the integer PBS Job ID. This ID is stored in memory as `current_job_id` so the orchestrator can track or cancel the job if interrupted.

##### 3. Job Monitoring and Filesystem Latency Verification (`monitor_job_completion`)

The orchestrator monitors job progression using a polling loop:

* **Queue Inspection**: Every 5 minutes, the script queries `qstat -a` filtered by the job ID. If the status returned is `Q` (Queued), `H` (Held), or `R` (Running), the orchestrator resets its strike counter and initiates another 5-minute countdown.

* **Shared Storage Latency Mitigation**: On large clusters, parallel filesystems (such as Lustre or GPFS) frequently suffer from metadata caching delays where files written by compute nodes do not appear immediately to login nodes. Once a job is no longer listed in `qstat`, the orchestrator does not assume immediate failure. It tests whether all entries in `sync.critical_files` (such as `case_results.pickle` and `case_restart_ts100.pickle`) are physically present on disk using remote file tests:

```bash
[ -f case_results.pickle ] && [ -f case_restart_ts100.pickle ] && echo 'ALL_FOUND'

```

* **Strike-Counter Guard**: If the job has disappeared from the queue but the critical files are not yet visible, the orchestrator registers a strike. It allows up to 6 consecutive strikes (30 minutes total) to provide sufficient time for remote filesystem buffers to flush. If files remain missing after 6 strikes, a `RuntimeError` is raised to trigger the crash handling routine.

##### 4. Post-Execution Synchronisation (`rsync_from_hpc`)

Once all critical files are confirmed present, the controller synchronises the results back to the local shared storage mount:

* Both solver directories (`CFD_1` and `CFD_2`) are updated locally using `rsync -avz --delete`.

* The critical pickle files (`case_results.pickle`, `case_restart*.pickle`) are copied into the simulation root directory to enable downstream processing and archiving.

##### 5. Remote Scratch Sanitisation (`delete_files_on_hpc`)

After data is safely mirrored locally, the controller removes old `critical_files` from the HPC scratch directory using `rm -f`. This step prevents future runs from falsely detecting old result files if a subsequent run crashes before writing its own data.

---

#### Remote Failure Handling and Cancellation

If a remote job crashes or if the user halts the local orchestrator using a keyboard interrupt (`Ctrl+C`), the orchestrator catches the exception inside `handle_crash()`.

If `hpc.enabled` is active and `current_job_id` is set, `cancel_hpc_job(job_id, config)` executes `qdel <job_id>` over SSH immediately. This stops remote compute jobs from continuing to consume compute allocations in the background after the local orchestrator has aborted.

---

## 2. Comprehensive Lifecycle and Execution Phases

The lifecycle progresses through three primary phases, governed by the central orchestrator (`orchestrator.py`) and indexed within `state.json`.

```
====================================================================================================
RUN TIMELINE & STATE PROGRESSION
====================================================================================================
state.json: {"current_run": 0}  --> Phase A: Initial Case Setup & Run 0
state.json: {"current_run": 1}  --> Phase B: Remesh Loop (Run 1 -> run_1 created, Run 2 computed)
...
state.json: {"current_run": 28} --> Phase B: Remesh Loop (Run 28 -> run_28 created, Run 29 computed)
state.json: {"current_run": 60} --> Phase C: Finalisation (Data collected, run_60 created)
====================================================================================================

```

### Phase A: Setup and Run 0 Initialisation

When `state['current_run'] == 0`, the orchestrator initialises the domain from scratch:

1. **Case Directory Creation (`src/setup_case.py`)**: Runs on `mac_setup`. It cleans existing `CFD_1` and `CFD_2` directories, copies the solid configuration from `setup_files/solid` into `CFD_1`, and populates `CFD_2` with liquid case files (`case.jou`, `eicosane`, `setup_fluent.sh`, `background.msh`, `component.msh`).
2. **Fluent Initialisation Batch (`setup_fluent.sh` calling `case.jou`)**:
* Assembles the overset grid system by reading Cartesian background mesh `background.msh` and appending component mesh `component.msh`.
* Defines physical zones, assigning `fluid.bg` and `fluid.com` to the fluid domain, `outer` to the overset boundary type, and creating the overset interface between background and component.
* Sets material properties and initialises the flow on the static overset grid
* Activates dynamic meshing on `fluid.com` using radial smoothing controls, establishes surface and volume monitors, exports initial boundary coordinates to `node_coords_inner`, and saves `pc_fluid.cas.h5`.


3. **Solid Interface Initialisation (`src/solid_interface.py`)**: Reads `CFD_2/node_coords_inner`, strips duplicate entries, formats $(x, y)$ coordinates, and exports `CFD_1/interface_nodes.csv`.
4. **Parameter Initialisation**: `prepare_parameters(config, is_restart=False)` sets `case_file` to `pc_fluid.cas.h5` and `restart_rb` to `0`.
5. **Execution**: Dispatches `src/run_simulation.py` to solve the initial chunk. Upon successful completion, the state updates to `state['current_run'] = 1`.

### Phase B: Cyclic Remeshing and Solution Loop

For every step from `start_run` to `total_runs - 1`, the orchestrator executes a sequence of data preservation, domain deconstruction, meshing, reconstruction, and solving:

```
[Start Run i] 
      │
      ├──> 1. collect_data.py (Aggregate log histories into post_processing/)
      │
      ├──> 2. remesh_pipeline.py (Archive to run_{i-1}, run step_1, mesh, run step_2)
      │
      ├──> 3. solid_interface.py (Update CFD_1/interface_nodes.csv)
      │
      ├──> 4. prepare_parameters() (Set parameters.json to restart mode)
      │
      └──> 5. run_simulation.py (Solve a fixed time interval via CoCoNuT)
            │
            ├── Success: state['current_run'] = i + 1  --> Proceed to next run
            └── Failure: handle_crash()                --> Manual recovery mode

```

1. **Continuous Data Harvesting (`src/collect_data.py`)**: Runs on `mac_process`. It reads transient report files generated in the completed run (`CFD_2/report-file.out` and `CFD_2/rigid-body-report-file.out`) and appends them to cumulative files in `post_processing/`.
2. **Remeshing Pipeline (`src/remesh_pipeline.py`)**: Dispatched to `mac_remesh`. Handles folder archiving, case stripping, GAMBIT mesh generation, and field interpolation.
3. **Interface Update (`src/solid_interface.py`)**: Regenerates `CFD_1/interface_nodes.csv` using the boundary output of the newly interpolated mesh.
4. **Parameter Update**: `prepare_parameters(config, is_restart=True)` sets the liquid solver case to `case_restart.cas.h5` and flags rigid body motion for restart mode (`restart_rb = 50` or `100`).
5. **Coupled Solution (`src/run_simulation.py`)**: Executes CoCoNuT on `mac_setup`, advancing the solution over the prescribed time-step block.
6. **State Persistence**: On clean termination, `state['current_run']` increments by 1, and the loop advances.

### Phase C: Finalisation and Archiving

When `state['current_run'] == total_runs`, the `while` loop terminates. The final completed run has not yet been processed by the start-of-loop hooks, so the orchestrator invokes a dedicated closing sequence:

1. **Final Data Harvesting**: Invokes `collect_data.py` on `mac_process` to append the final run time steps to the master logs in `post_processing/`.
2. **Terminal Archiving (`src/archive_run.py`)**: Evaluates `state.json`, constructs the directory `run_{total_runs}`, and performs a non-destructive copy of `CFD_1`, `CFD_2`, and final pickled solution states (`case_results.pickle`, `case_restart*.pickle`).

---

## 3. Deep Dive: Remeshing Subsystem Architecture

The remeshing subsystem operates via two distinct execution paths: an automated path driven by GAMBIT script generation, and a manual path using ANSYS DesignModeler.

```
====================================================================================================
FLUENT RESTRUCTURING & INTERPOLATION PIPELINE
====================================================================================================
Active Case (CFD_2)                Fluent step_1.jou                     remesh/case_restart.cas
[fluid.bg + fluid.com] ──> [Export inner coords]              ──> [Naked component geometry]
                           [Unload udf_thermal]                             │
                           [Delete overset interfaces]                      │ (Meshing Engine)
                           [Delete fluid.bg cell zone]                      ▼
                           [Slit wall: internal -> interior]      remesh/remesh.msh
                                                                            │
New Active Case (CFD_2)            Fluent step_2.jou                        │
[Interpolated Fields]  <── [Export new inner coords]          <─────────────┘
                           [Re-create overset interfaces]
                           [Read interpolated sol_dat.ip]
                           [Replace zone: remesh.msh]
====================================================================================================

```

### The Fluent Two-Step Journal Sequence

#### Step 1 Deconstruction (`step_1.jou`)

The component mesh cannot be remeshed while connected to the overset background domain in Fluent. `step_1.jou` strips the background elements away:

1. Loads the latest case and data file: `case_restart.cas.h5`.
2. Exports pre-remesh inner wall node coordinates to `node_coords_inner_before`.
3. Deletes dynamic mesh zones associated with `inner`.
4. Unloads the compiled C UDF library (`udf_thermal`) to release shared library memory locks.
5. Deletes all overset interface pairings.
6. Saves `case_restart.cas.h5` with full solution field data preserved.
7. Executes `mesh modify-zones delete-cell-zone fluid.bg () b` to delete the Cartesian background cell zone entirely.
8. Changes boundary condition types: switches `internal` to `wall`, slitting the interface into coupled walls (`internal` and `internal-shadow`).
9. Disables CFF file formatting (`file cff-files? n`) and writes the stripped component mesh to legacy format: `remesh/case_restart.cas`.

#### Step 2 Reconstruction and Interpolation (`step_2.jou`)

Once a new mesh file (`remesh/remesh.msh`) exists, `step_2.jou` reassembles the overset domain:

1. Re-opens the saved `case_restart.cas.h5` (which still contains both `fluid.bg` and the old `fluid.com`).
2. Interpolates current field solution variables (pressure, velocities, temperature) to file: `file interpolate write-data sol_dat.ip y y`.
3. Deletes active dynamic mesh specifications on `outer` and `fluid.com`.
4. Executes zone replacement: `mesh modify-zones replace-zone "remesh/remesh.msh" fluid.com fluid.com y`. Fluent matches boundary names (`inner`, `outer`, `internal`), deletes obsolete cells, inserts the new mesh topology, and automatically transfers boundary conditions.
5. Restores dynamic mesh zones: recreates `fluid.com` as a radial deforming zone and `outer` as an unspecified deforming boundary.
6. Restores boundary zone `outer` to type `overset`.
7. Re-establishes the overset interface pairing between cell zones `fluid.bg` and `fluid.com`.
8. Re-initialises flow and reads back field data from `sol_dat.ip`.
9. Restores volume and surface report definitions (`vmax`, `htr-int`).
10. Exports updated boundary node coordinates to `node_coords_inner`.
11. Writes updated solver files: `case_restart.cas.h5` and `case_restart.dat.h5`.

---

### Automated Meshing Implementation

The automated meshing engine relies on GAMBIT 2.4.6 running in batch mode, orchestrated by `setup_files/liquid/create_mesh/remesh.sh`:

```
+----------------------------------------------------------------------------------------------------+
|                                    AUTOMATED GAMBIT ENGINE FLOW                                    |
+----------------------------------------------------------------------------------------------------+
| 1. Obtain File Lock: /tmp/gambit_license.lock via flock 200 on cfdclu13                            |
| 2. Purge Stale Databases: rm -f identified_edges.* remesh.dbs remesh.jou                          |
| 3. Geometric Summarisation: gambit -inp identifyer.jou -> exports identified_edges.trn            |
| 4. Python Topological Analysis: generate_remesh.py reads .trn -> writes remesh.jou                 |
| 5. Mesh Generation: gambit -inp remesh.jou -> exports remesh.msh                                  |
| 6. Release File Lock: flock -u 200                                                                 |
+----------------------------------------------------------------------------------------------------+

```

In essence, the automated path in GAMBIT, strips away the existing, degenerated component mesh and created a new internal mesh within the same boundaries (solid interface and `outer` overset boundary).
The overset boundary cannot be regenerated in an automated way in GAMBIT and will eventually degenerate as well, necessitating the occasional manual remeshing intervention in Workbench.

#### Topological Edge Identification (`generate_remesh.py`)

Because the melting interface deforms non-uniformly, GAMBIT assigns dynamic edge names during file import. `generate_remesh.py` parses `identified_edges.trn` to extract edge names, lengths, boundary assignments, and node sizing. The values below are given as example to this specific simulation, but differ depending on the simulation:

* **Internal Edge Grading Reversal**: Internal edges connecting inner and outer boundaries require geometric expansion away from the melting wall. The script reads the first and last cell dimensions. If `first_size > last_size`, it appends the `reverse` flag to the GAMBIT meshing command (`ratio1 1.127 intervals 13 reverse`) to maintain appropriate boundary layer growth directions.
* **Length-Segregated Interface Sizing**: Inner melting edges evolve varying arc lengths. The script evaluates edge lengths against a threshold of $0.005\text{ m}$:
* Long edges ($L > 0.005\text{ m}$): Meshed with edge size $4.419\times 10^{-5}\text{ m}$.
* Short edges ($L \le 0.005\text{ m}$): Meshed with refined edge size $2.121\times 10^{-5}\text{ m}$.


* **Centroid-Based Edge Linking**: Outer overset boundary edges must maintain structured alignment with corresponding inner edges. The script computes Euclidean distances between edge centroids:

$$\text{dist} = \sqrt{(x_{\text{out}} - x_{\text{in}})^2 + (y_{\text{out}} - y_{\text{in}})^2 + (z_{\text{out}} - z_{\text{in}})^2}$$

Each outer edge links to its nearest inner counterpart (`edge link "outer_e" "inner_e"`), forcing uniform structured sub-mapping across quadrilaterals.

* **Physics Definition**: Modifies boundary entities, explicitly assigning `internal` edges to `INTERIOR`, and `outer` and `inner` edges to `WALL`.

---

### Manual Recovery Mode and DesignModeler Pipeline

When automated meshing fails due to severe overset boundary distortions, the orchestrator routes execution into manual recovery mode.

```
====================================================================================================
MANUAL INTERVENTION RECOVERY FLOW
====================================================================================================
Simulation Crash detected in Run i
  │
  ├──> 1. clean_local_fluent() (Terminates hanging background processes)
  │
  ├──> 2. User selects Option [2]
  │
  ├──> 3. export_designmodeler.py runs:
  │         Reads CFD_1/interface_nodes.csv
  │         Nearest-neighbour sorting of point cloud
  │         Adaptive RDP decimation until node count < 2000
  │         Writes manual_interventions/interface_nodes_designmodeler.txt
  │
  ├──> 4. User opens DesignModeler, loads 3D Curve, builds new mesh,
  │         and exports to manual_interventions/new_mesh.msh
  │
  ├──> 5. User presses ENTER
  │
  ├──> 6. remesh_pipeline.py --manual-mesh executes:
  │         Pulls clean restart files from run_{i-1}/CFD_2
  │         Injects new_mesh.msh into CFD_2/remesh/
  │         Executes step_1.jou and step_2.jou (Gambit bypassed)
  │
  ├──> 7. solid_interface.py updates CFD_1/interface_nodes.csv
  │
  ├──> 8. prepare_parameters() updates parameters.json
  │
  └──> 9. Orchestrator sets skip_remesh = True and re-executes Run i
====================================================================================================

```

#### Adaptive Curve Decimation (`src/export_designmodeler.py`)

ANSYS DesignModeler enforces a strict upper limit of 2000 points per imported 3D curve file. The simulation interface often contains 4000 to 6000 nodes. `src/export_designmodeler.py` transforms raw interface coordinates into compliant curves:

1. **Point Cloud Sorting**: Uses nearest-neighbour distance sorting to assemble an ordered traversal path along the interface perimeter, resolving random node ordering exported from unstructured boundary tables.
2. **Median Spacing Baseline**: Calculates the median Euclidean spacing between consecutive points:

$$\Delta s_k = \sqrt{(x_{k+1} - x_k)^2 + (y_{k+1} - y_k)^2}, \quad \bar{s} = \text{median}(\Delta s)$$

3. **Adaptive Ramer-Douglas-Peucker (RDP) Algorithm**: Line simplification applies a chordal perpendicular distance metric:

$$d_{\perp} = \frac{\vert{}(x_2 - x_1)(y_1 - y_0) - (y_2 - y_1)(x_1 - x_0)\vert{}}{\sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2}}$$

Points exceeding threshold $\epsilon = \text{TOLERANCE\_FACTOR} \times \bar{s}$ are retained. If the resulting point count exceeds 1999 nodes, the script executes an automated adaptation loop:

```python
while len(filtered_points) > DESIGNMODELER_NODE_LIMIT and len(filtered_points) > 2:
    adjustment_steps += 1
    tolerance_factor *= 1.20
    epsilon = tolerance_factor * median_spacing
    filtered_points = ramer_douglas_peucker(sorted_points, epsilon)

```

4. **DesignModeler 3D Curve Formatting**: Writes the filtered nodes using the format required by DesignModeler (`Group_ID Point_ID X Y Z`):
```text
1 1 -0.026000 0.052000 0
1 2 -0.025937 0.052000 0
...
1 0

```



#### Safe State Restoration in `remesh_pipeline.py`

When executed with the `--manual-mesh` argument, `remesh_pipeline.py` modifies its operation to protect simulation integrity:

* **Archive Protection**: Bypasses the archiving step in Step 3, preventing the corrupted or crashed active directory from overwriting earlier clean archives.
* **Restart Retrieval**: Instead of reading corrupted files from `CFD_2`, it sets `search_dir` to `run_{max_run_num}/CFD_2`, retrieving intact `case_timestep*.cas.h5` files from the last successful run.
* **UDF Library Preservation**: Restores compiled `udf_thermal` directories and source files from the clean archive.
* **Gambit Bypass**: Copies `new_mesh.msh` directly to `CFD_2/remesh/remesh.msh` and skips `remesh.sh`.

---

## 4. Error Handling and State Integrity

The framework incorporates multiple safeguards against process crashes, license conflicts, and loop desynchronisation.

```
+----------------------------------------------------------------------------------------------------+
|                                      CONCURRENCY & SAFETY LOCKS                                    |
+----------------------------------------------------------------------------------------------------+
| Concurrency Control   : flock 200 on /tmp/gambit_license.lock serialises GAMBIT execution          |
| Failure Detection     : Try/Except catches CalledProcessError, RuntimeError, KeyboardInterrupt     |
| Process Sanitisation  : clean_local_fluent executes cleanup scripts and kills apip-standalone      |
| Loop Desynchronisation: skip_remesh flag prevents automated overwriting after manual recovery      |
| State Atomicity       : state.json increments only after run_simulation.py exits cleanly           |
+----------------------------------------------------------------------------------------------------+

```

### 1. License Concurrency Control

Because GAMBIT operates under a single floating network license, concurrent execution attempts on `cfdclu13` result in immediate failure. `remesh.sh` wraps all GAMBIT operations inside an exclusive Linux file descriptor lock:

```bash
exec 200>"/tmp/gambit_license.lock"
flock 200

```

If a secondary simulation instance attempts to remesh simultaneously, the shell blocks execution at `flock 200` until the active instance exits and closes file descriptor 200.

### 2. Overwrite Protection via `skip_remesh`

In `orchestrator.py`, recovering from a crash using a manual mesh introduces a loop sequencing challenge. If the loop simply restarts, it would run `remesh_pipeline.py` in automated mode, immediately overwriting the manual mesh.

The orchestrator prevents this using a control flag:

```python
try:
    if not skip_remesh:
        execute_local_node("collect_data.py", mac_process, base_dir)
        execute_local_node("remesh_pipeline.py", mac_remesh, base_dir)
        execute_local_node("solid_interface.py", mac_process, base_dir)
        prepare_parameters(config, is_restart=True)
    skip_remesh = False
    execute_local_node("run_simulation.py", mac_setup, base_dir)
    state['current_run'] = i + 1
    save_state(state)
except Exception as e:
    recovered = handle_crash(i, mac_setup, mac_process, config, current_job_id, e)
    if recovered:
        skip_remesh = True

```

When manual recovery completes inside `handle_crash()`, it sets `skip_remesh = True`. The next iteration bypasses the automated remeshing pipeline, proceeding directly to `run_simulation.py` with the manual mesh in place.

### 3. Data Accumulation Mechanics (`collect_data.py`)

Because the coupled solver restarts repeatedly, local report logs (`CFD_2/report-file.out` and `CFD_2/rigid-body-report-file.out`) reset their row counts and relative timestamps at every run:

* **Overcoming Row 0 Overlap**: Fluent writes time step 0 upon restarting. `collect_data.py` discards row 0 (`df_local[df_local[0] > 0]`) to avoid duplicate timestamp entries.
* **Timeline Accumulation**: Extracts the final time step and flow time from the master log in `post_processing/`, adds them as offsets to the incoming local data, and appends the corrected entries.
* **Rigid Body Angular Unwrapping**: Rigid body rotation angle $\theta_z$ is offset by the final accumulated angle of the previous run, maintaining a continuous angular displacement history across remesh boundaries.

---

## 5. Complete File, Script, and Configuration Reference

```
HIRATA WORKFLOW REPOSITORY ROOT
├── config.yaml
├── orchestrator.py
├── parameters.json
├── state.json
├── manual_interventions/
│   ├── interface_nodes_designmodeler.txt
│   └── new_mesh.msh
├── post_processing/
│   ├── report-file.out
│   └── rigid-body-report-file.out
├── setup_files/
│   ├── solid/
│   │   └── solver_parameters.json
│   └── liquid/
│       ├── case.jou
│       ├── eicosane
│       ├── setup_fluent.sh
│       ├── background.msh
│       ├── component.msh
│       ├── step_1.jou
│       ├── step_2.jou
│       └── create_mesh/
│           ├── identifyer.jou
│           ├── generate_remesh.py
│           ├── remesh.sh
│           ├── create_mesh.sh
│           └── mesh_background.jou
└── src/
    ├── archive_run.py
    ├── collect_data.py
    ├── export_designmodeler.py
    ├── remesh_pipeline.py
    ├── run_simulation.py
    ├── setup_case.py
    ├── solid_interface.py
    └── utils.py

```

| File / Script Path | Host Node | Primary Responsibility | Input Artifacts | Output Artifacts |
| --- | --- | --- | --- | --- |
| `config.yaml` | All | Master workflow configuration | N/A | Configuration dictionary |
| `state.json` | All | Workflow progression counter | N/A | Persistent integer `current_run` |
| `parameters.json` | All | CoCoNuT master FSI configuration | `config.yaml` | Read by `run_simulation.py` |
| `orchestrator.py` | Controller | Master automation loop & crash handler | `config.yaml`, `state.json` | Drives entire simulation lifecycle |
| `src/utils.py` | Controller | SSH dispatch, cleanup, parameter update | Host configuration | Remote processes, modified JSON |
| `src/setup_case.py` | `local_setup` | Sets up clean directories for Run 0 | `setup_files/` | Fresh `CFD_1/`, `CFD_2/`, `pc_fluid.cas.h5` |
| `src/collect_data.py` | `local_process` | Aggregates discontinuous run output logs | `CFD_2/*report*.out` | Continuous `post_processing/*.out` |
| `src/remesh_pipeline.py` | `remesh` | Drives archiving, step 1, meshing, step 2 | `case_restart.cas.h5` | Archived `run_x/`, `case_restart.cas.h5` |
| `src/solid_interface.py` | `local_process` | Extracts boundary points for solid solver | `CFD_2/node_coords_inner` | `CFD_1/interface_nodes.csv` |
| `src/export_designmodeler.py` | `local_process` | Decimates interface curve for DesignModeler | `CFD_1/interface_nodes.csv` | `manual_interventions/*_designmodeler.txt` |
| `src/archive_run.py` | `local_process` | Finalises and archives the terminal run | Active `CFD_1/`, `CFD_2/` | Final directory `run_{total_runs}/` |
| `src/run_simulation.py` | `local_setup` | Executes coupled FSI solution chunk | `parameters.json` | Solved time steps, result pickles |
| `setup_files/liquid/step_1.jou` | `remesh` | Strips background mesh from Fluent case | `case_restart.cas.h5` | `remesh/case_restart.cas` |
| `setup_files/liquid/step_2.jou` | `remesh` | Injects new mesh and maps solution | `remesh/remesh.msh` | Updated `case_restart.cas.h5` |
| `setup_files/liquid/create_mesh/remesh.sh` | `remesh` | Coordinates GAMBIT batch meshing | `identifyer.jou`, `case_restart.cas` | `remesh.msh` |
| `setup_files/liquid/create_mesh/identifyer.jou` | `remesh` | Summarises edge topology in GAMBIT | `case_restart.cas` | `identified_edges.trn` |
| `setup_files/liquid/create_mesh/generate_remesh.py` | `remesh` | Generates custom GAMBIT meshing journal | `identified_edges.trn` | `remesh.jou` |


## References

<a id="1">[1]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “A partitioned interface-tracking method for close-contact melting of unconstrained saturated solids,” submitted to Int. J. Heat Mass Transfer., 2026.]

<a id="2">[2]</a>
[V. Van Riet, “The Development of a Partitioned Simulation Method for Close-Contact Melting of Unconstrained Solids,” Ghent University. Faculty of Engineering and Architecture, Ghent, Belgium, 2026.](http://hdl.handle.net/1854/LU-01M4T9NN21ND2PG4SZQW47D2XT)
