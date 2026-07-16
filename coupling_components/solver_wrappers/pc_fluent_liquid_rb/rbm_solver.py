import os
from os.path import join
import numpy as np
import pickle
from scipy.spatial import cKDTree


class RBMSolver:
    def __init__(self, rb_settings, pc_settings, contact_settings, delta_t, dimensions, dir_cfd, restart):
        self.delta_t = delta_t
        self.dimensions = dimensions
        self.dir_cfd = dir_cfd
        self.restart = restart

        # --- Rigid Body Motion Specific Settings ---
        self.restart_rb_only = rb_settings.get('restart_rb', 0)
        self.rb_tol = rb_settings.get('tolerance', 1E-6)
        self.rb_iter_min = rb_settings.get('iteration_min', 1)
        self.rb_iter_max = rb_settings.get('iteration_max', 10)

        if self.rb_iter_min > self.rb_iter_max:
            self.rb_iter_min = 1
            print('Info: rb_iter_min is set larger than rb_iter_max, rb_iter_min has been reset to 1.')

        self.rb_relax = rb_settings.get('relaxation', 1.0)  # Default: no relaxation
        self.fict_coeff = rb_settings.get('fict_coeff', 0.0)  # Default: no fictitious impedance
        self.multiplier = rb_settings.get('fict_multiplier', 1.0)  # Default: no multiplier for first timestep
        self.dyn_visc = rb_settings.get('liquid_dyn_visc', 0.0)  # Dynamic viscosity of liquid PCM
        self.fm_wall = rb_settings.get('fm_wall', None)
        self.rb_relax_method = rb_settings.get('relaxation_method', 'static')  # Options: 'static' or 'aitken'
        self.relax_first_ts = rb_settings.get('relax_first_ts', False)
        self.rot_update = rb_settings.get('rotational_update',
                                          'quaternion')  # Options: 'off', 'rot_mat' or 'quaternion'
        self.rb_predictor = rb_settings.get('predictor', 'constant')  # Options: 'constant' or 'linear'
        self.buoyancy = rb_settings.get('buoyancy', True)
        self.x_motion = rb_settings.get('x_motion', True)
        self.weight_ramp = rb_settings.get('weight_ramp', 0)  # Nr. of time steps over which full weight is added
        self.gravity = rb_settings.get('gravity', [0, -9.81, 0])

        if self.restart_rb_only != 0:
            self.weight_ramp = 0

        # --- Phase Change Specific Settings ---
        self.volume_change = pc_settings.get('volume_change', True)  # Account for volume change during melting
        self.liquid_density = pc_settings.get('liquid_density', None)
        self.solid_density = pc_settings.get('solid_density', 0.0)

        # Assume equal density if solid density is not explicitly provided
        if self.solid_density == 0.0 and self.liquid_density is not None:
            self.solid_density = self.liquid_density

        # --- Contact Model Specific Settings ---
        self.include_contact_force = contact_settings.get('contact_force', False)
        self.h_ul = contact_settings.get('upper_limit', 2e-4)  # [m] Upper limit for danger zone
        self.h_ll = contact_settings.get('lower_limit', 2e-4)  # [m] Lower limit for hard stop
        self.damping_ratio = contact_settings.get('damping_ratio', 1.0)
        self.k_mass = contact_settings.get('k_mass', 1.0)

        # --- Initialise Kinematics ---
        self.v_trans_prev = np.zeros(3)  # Previous timestep translational velocity
        self.v_trans = np.zeros(3)  # Current translational velocity
        self.omega_prev = np.zeros(3)  # Previous timestep rotational velocity
        self.omega = np.zeros(3)  # Current rotational velocity
        self.a_trans_prev = np.zeros(3)  # Previous timestep translational acceleration
        self.a_trans = np.zeros(3)  # Current translational acceleration
        self.a_rot_prev = np.zeros(3)  # Previous timestep rotational acceleration
        self.a_rot = np.zeros(3)  # Current rotational acceleration
        self.com = np.zeros(3)  # Current timestep center of mass
        self.com_prev = np.zeros(3)  # Previous timestep center of mass

        # Initialise rotational tracking matrices/quaternions if not restarting
        if not self.restart:
            if self.rot_update == 'rot_mat':
                self.orientation_prev = np.identity(3)
                self.orientation = np.identity(3)
            else:
                self.orientation_prev = np.array([1.0, 0.0, 0.0, 0.0])
                self.orientation = np.array([1.0, 0.0, 0.0, 0.0])

        # --- State Tracking Variables ---
        self.volume = 0.0
        self.M_sys = 0.0
        self.a_trans_prev_it = np.zeros(3)
        self.force_int = np.zeros(3)  # Integrated fluid forces
        self.moment_int = np.zeros(3)  # Integrated fluid moments
        self.force_pr_it = np.zeros(3)  # Relaxed force from previous iteration
        self.moment_pr_it = np.zeros(3)  # Relaxed moment from previous iteration
        self.h_min = 0.0
        self.contact_patches = []
        self.prev_patches = []
        self.new_patches = []
        self.gap_trees = {}
        self.wall_length = 0.0

        # --- Aitken Relaxation State Variables ---
        if self.rb_relax_method == 'aitken':
            self.force_res_prev = np.zeros(3)
            self.moment_res_prev = np.zeros(3)
            self.aitken_relax_factor_force = self.rb_relax
            self.aitken_relax_factor_moment = self.rb_relax

    @property
    def avg_v_trans(self):
        """Average velocity over the time step for 2nd order accurate position update."""
        return 0.5 * (self.v_trans_prev + self.v_trans)

    @property
    def avg_omega(self):
        """Average angular velocity over the time step."""
        return 0.5 * (self.omega_prev + self.omega)

    def initialise_walls(self, wall_coords):
        """Build KD-Trees for fast distance queries to boundary walls."""
        for wall_name, coords in wall_coords.items():
            wall_tree = cKDTree(coords)
            print(f"KD-Tree initialised with {len(coords)} nodes for wall {wall_name}.")
            self.gap_trees[wall_name] = wall_tree

            # Calculate Length of Correct Wall via Greedy Walk for fictitious mass
            if wall_name == self.fm_wall:
                curr_idx = np.argmin(coords[:, 0])  # Start at min X
                visited = {curr_idx}
                l_current = 0.0

                while len(visited) < len(coords):
                    # Query neighbors (k=10 buffer for visited nodes)
                    dists, idxs = wall_tree.query(coords[curr_idx], k=10)

                    # Find nearest unvisited neighbor
                    next_step = next(((d, i) for d, i in zip(dists, idxs) if i not in visited), None)

                    if next_step:
                        dist, idx = next_step
                        l_current += dist
                        visited.add(idx)
                        curr_idx = idx
                    else:
                        print(f"Warning: Discontinuity in wall {wall_name}")
                        break

                self.wall_length = l_current
                print(f"Wall Length: {self.wall_length:.6f} [m]")

    def initialise_solution_step(self, timestep):
        """Update historical variables at the start of a new timestep."""
        if self.rb_predictor == 'linear':
            v_trans_prev2 = self.v_trans_prev.copy()
            omega_prev2 = self.omega_prev.copy()

        # Shift current values to previous
        self.v_trans_prev = self.v_trans.copy()
        self.omega_prev = self.omega.copy()
        self.a_trans_prev = self.a_trans.copy()
        self.a_rot_prev = self.a_rot.copy()
        self.orientation_prev = self.orientation.copy()
        self.com_prev = self.com.copy()
        self.prev_patches = self.new_patches.copy()

        # Reset Aitken relaxation residuals for the new time step
        if self.rb_relax_method == 'aitken':
            self.force_res_prev = np.zeros(3)
            self.moment_res_prev = np.zeros(3)
            self.aitken_relax_factor_force = np.clip(self.aitken_relax_factor_force, 0.01, self.rb_relax)
            self.aitken_relax_factor_moment = np.clip(self.aitken_relax_factor_moment, 0.01, self.rb_relax)

        # Linear predictor step for v_trans & omega_z
        if self.rb_predictor == 'linear' and timestep > 1:
            self.v_trans = 2 * self.v_trans_prev - v_trans_prev2
            self.omega = 2 * self.omega_prev - omega_prev2

            # Make accelerations consistent with the predicted velocities
            self.a_trans = 2 * (self.v_trans - self.v_trans_prev) / self.delta_t - self.a_trans_prev
            self.a_rot = 2 * (self.omega - self.omega_prev) / self.delta_t - self.a_rot_prev

    def step(self, force_raw, moment_raw, volume, com, moi, r_interface, timestep, iteration, rb_iter):
        """Core physics solver step representing one internal rigid body iteration."""
        self.volume = volume
        self.com = com

        if timestep == 1 and iteration == 1 and rb_iter == 1:
            self.com_prev = self.com.copy()

        # --- Force & Moment Relaxation ---
        # Option to skip relaxation on the very first timestep and iteration to accelerate convergence
        if not self.relax_first_ts and timestep == 1 and iteration == 1:
            force_relaxed = force_raw.copy()
            moment_relaxed = moment_raw.copy()
        else:
            if self.rb_relax_method == 'static':
                # Simple under-relaxation
                force_relaxed = self.force_pr_it + self.rb_relax * (force_raw - self.force_pr_it)
                moment_relaxed = self.moment_pr_it + self.rb_relax * (moment_raw - self.moment_pr_it)
            elif self.rb_relax_method == 'aitken':
                # Aitken on the force and moment residuals
                self.aitken_relax_factor_force = self._calc_aitken(
                    force_raw, self.force_pr_it, self.force_res_prev, self.aitken_relax_factor_force, rb_iter
                )
                self.force_res_prev = force_raw - self.force_pr_it
                force_relaxed = self.force_pr_it + self.aitken_relax_factor_force * (force_raw - self.force_pr_it)

                self.aitken_relax_factor_moment = self._calc_aitken(
                    moment_raw, self.moment_pr_it, self.moment_res_prev, self.aitken_relax_factor_moment, rb_iter
                )
                self.moment_res_prev = moment_raw - self.moment_pr_it
                moment_relaxed = self.moment_pr_it + self.aitken_relax_factor_moment * (moment_raw - self.moment_pr_it)
            else:
                # Fallback: no relaxation
                force_relaxed = force_raw.copy()
                moment_relaxed = moment_raw.copy()

        # Store relaxed values as "previous" for the next iteration
        self.force_pr_it = force_relaxed.copy()
        self.moment_pr_it = moment_relaxed.copy()
        self.force_int = force_raw.copy()
        self.moment_int = moment_raw.copy()

        # Calculate fraction of weight that is accounted for (ramp-up to prevent initial instability)
        self.weight_factor = timestep / self.weight_ramp if timestep < self.weight_ramp else 1

        # --- Translational Update ---
        g = np.array(self.gravity)
        mass_solid = self.solid_density * self.volume

        # Fictitious Mass / Damping method to anticipate high added mass & viscous damping in fluid solver
        if self.gap_trees and r_interface is not None:
            h_used, self.h_min = self._calculate_h_min(r_interface)

            print("\n")
            print(f'Min. gap width = {self.h_min * 1000} mm')
            print(f'Number of active contact patches = {len(self.contact_patches)}')

            # Singularity protection for fictitious mass
            epsilon = 1e-4 * self.wall_length
            h_eff = max(h_used, epsilon) if h_used is not None else epsilon
            depth = 1.0

            fict_mass = self.liquid_density * (self.wall_length ** 3) * depth / h_eff
            fict_damping = depth * (self.wall_length ** 3) * (
                        self.solid_density / self.liquid_density) * self.dyn_visc / (h_eff ** 3)
        else:
            fict_mass = 0
            fict_damping = 0
            self.h_min = 0.0

        self.a_trans_prev_it = self.a_trans.copy()

        # Add buoyancy or direct gravity
        if self.buoyancy:
            F_net = force_relaxed + self.weight_factor * (self.solid_density - self.liquid_density) * g * self.volume
        else:
            F_net = force_relaxed + self.weight_factor * mass_solid * g

        # Include explicit contact forces if configured
        contact_F = np.zeros(3)
        contact_M = np.zeros(3)
        if self.include_contact_force and self.gap_trees:
            contact_F, contact_M = self._calc_contact_force_and_moment()

        F_net += contact_F

        # Apply Fictitious Mass: We add (M_sys * a_prev_iter) to the forces.
        # This dampens the acceleration update significantly without altering the final converged physics.
        self.M_sys = fict_mass + self.delta_t * fict_damping / 2
        M_sys_raw = self.M_sys
        self.M_sys *= self.fict_coeff
        if rb_iter == 1:
            self.M_sys *= self.multiplier

        F_tot = F_net + self.M_sys * self.a_trans_prev_it

        # Update velocity with Crank-Nicolson integration (2nd order accurate)
        self.a_trans = F_tot / (mass_solid + self.M_sys)
        self.v_trans = self.v_trans_prev + 0.5 * (self.a_trans_prev + self.a_trans) * self.delta_t

        if not self.x_motion:
            self.v_trans[0] = 0.0

        # --- Translational Wall Guard (Failsafe) ---
        if self.h_min < self.h_ll and self.contact_patches:
            deepest_contact = self.contact_patches[0][0]
            normal = deepest_contact['normal']
            # If the body is moving into the wall beyond the lower limit, reject that motion vector
            if np.dot(self.v_trans, normal) < 0:
                print("*** HARD STOP ACTIVATED ***")
                rejected = np.dot(self.v_trans, normal) * normal
                self.v_trans = self.v_trans - rejected
                self.a_trans = 2 * (self.v_trans - self.v_trans_prev) / self.delta_t - self.a_trans_prev

        # --- Rotational Update ---
        if self.rot_update == 'off':
            # Enforce zero rotation when turned off
            self.a_rot = np.zeros(3)
            self.omega = np.zeros(3)
            self.orientation = np.array([1.0, 0.0, 0.0, 0.0])  # Identity quaternion
        else:
            # NOTE: This implementation assumes a DIAGONAL moment of inertia tensor
            fict_moi = (self.M_sys / (self.solid_density * self.volume)) * moi if (
                                                                                              self.solid_density * self.volume) != 0 else np.zeros(
                3)
            moi_total = moi + fict_moi
            a_rot_prev_it = self.a_rot.copy()
            moment_tot = moment_relaxed + contact_M + fict_moi * a_rot_prev_it

            # Calculate rotational acceleration
            self.a_rot = np.divide(moment_tot, moi_total, out=np.zeros_like(moment_tot), where=moi_total != 0)

            # Update rotational velocity with Crank-Nicolson integration
            self.omega = self.omega_prev + 0.5 * (self.a_rot_prev + self.a_rot) * self.delta_t

            # Rotational wall guard
            if self.h_min < self.h_ll and self.contact_patches:
                deepest_contact = self.contact_patches[0][0]
                normal = deepest_contact['normal']
                pinch_point = deepest_contact['point']
                rot_vel = np.cross(self.omega, (pinch_point - self.com))
                if np.dot(rot_vel, normal) < 0:
                    self.omega = np.zeros_like(self.omega)
                    # Make acceleration consistent
                    self.a_rot = 2 * (self.omega - self.omega_prev) / self.delta_t - self.a_rot_prev

            if self.rot_update == 'quaternion':
                # Update orientation using quaternions for second-order accuracy
                avg_omega = 0.5 * (self.omega_prev + self.omega)
                delta_quat = quat_from_angular_velocity(avg_omega, self.delta_t)
                self.orientation = quat_multiply(delta_quat, self.orientation_prev)
                self.orientation = quat_normalize(self.orientation)
            elif self.rot_update == 'rot_mat':
                # Update rotation matrix (assumes only rotation along z-axis)
                avg_omega = 0.5 * (self.omega_prev + self.omega)
                theta = self.delta_t * avg_omega
                self.orientation = self.orientation_prev @ exp_map(theta)

                # Re-orthonormalize if numerical drift exceeds tolerance
                if np.linalg.norm(self.orientation.T @ self.orientation - np.identity(3)) > 1e-8:
                    U, _, Vt = np.linalg.svd(self.orientation)
                    self.orientation = U @ Vt

    def rotate_displacement(self, disp_step_melting):
        """Convert melting displacement from local (solid) frame to global (liquid) frame."""
        if self.rot_update == 'rot_mat':
            return np.dot(disp_step_melting, self.orientation_prev.T)
        else:
            return quat_rotate_vector_array(self.orientation_prev, disp_step_melting)

    def write_report(self, timestep):
        """Create or update the rigid-body report file each timestep."""
        tmp = "rigid-body-report-file.out"
        file_name = join(self.dir_cfd, tmp)

        # If first timestep: create file and write header
        if timestep == 1:
            with open(file_name, "w") as f:
                f.write("# CoCoNuT rigid body motion history\n")
                f.write("#\n")
                f.write("#  {:>10}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}\n"
                        .format("time", "CG_X", "CG_Y", "V_X", "V_Y", "THETA_Z", "F_X", "F_Y", "M_Z", "h_min",
                                "volume"))
                f.write("#  {:>10}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}\n"
                        .format("(s)", "(m)", "(m)", "(m/s)", "(m/s)", "(deg)", "(N)", "(N)", "(N*m)", "(m)", "(m^3)"))
                f.write("#\n")

        time = timestep * self.delta_t
        volume = self.volume
        cg_x = self.com[0]
        cg_y = self.com[1]
        v_x = self.v_trans[0]
        v_y = self.v_trans[1]

        if self.rot_update == 'rot_mat':
            # Extract 2D rotation part if orientation is 3x3
            R = self.orientation[:2, :2]
            theta_rad = np.arctan2(R[1, 0], R[0, 0])
        else:
            # Extract the w and z components from the orientation quaternion
            w = self.orientation[0]
            z = self.orientation[3]
            theta_rad = 2 * np.arctan2(z, w)

        # Convert to degrees for reporting
        theta_deg = np.degrees(theta_rad)
        force_x = self.force_int[0]
        force_y = self.force_int[1]
        moment_z = self.moment_int[2]
        h = self.h_min

        with open(file_name, "a") as f:
            f.write(f"{time:12.5e}  {cg_x:12.5e}  {cg_y:12.5e}  "
                    f"{v_x:12.5e}  {v_y:12.5e}  {theta_deg:12.5e}  "
                    f"{force_x:12.5e}  {force_y:12.5e}  {moment_z:12.5e}  {h:12.5e}  {volume:12.5e}\n")

    def save_restart_data(self, timestep, interface_rb):
        """Save rigid body motion state to a pickle file for restart capabilities."""
        state = {
            'v_trans_prev': self.v_trans_prev,
            'v_trans': self.v_trans,
            'omega_prev': self.omega_prev,
            'omega': self.omega,
            'a_trans_prev': self.a_trans_prev,
            'a_trans': self.a_trans,
            'a_rot_prev': self.a_rot_prev,
            'a_rot': self.a_rot,
            'orientation_prev': getattr(self, 'orientation_prev', np.array([1.0, 0.0, 0.0, 0.0])),
            'orientation': getattr(self, 'orientation', np.array([1.0, 0.0, 0.0, 0.0])),
            'com': self.com,
            'com_prev': self.com_prev,
            'prev_patches': self.prev_patches,
            'new_patches': self.new_patches,
            'interface_rb': interface_rb,
            'force_pr_it': self.force_pr_it,
            'moment_pr_it': self.moment_pr_it
        }

        tmp = f'restart_rb_timestep{timestep}.pickle'
        file_name = join(self.dir_cfd, tmp)

        with open(file_name, 'wb') as f:
            pickle.dump(state, f)

    def load_restart_data(self, timestep_start):
        """Load rigid body motion state from a pickle file."""
        if self.restart:
            tmp = f'restart_rb_timestep{timestep_start}.pickle'
        elif self.restart_rb_only != 0:
            tmp = f'restart_rb_timestep{self.restart_rb_only}.pickle'

        file_name = join(self.dir_cfd, tmp)

        if not os.path.exists(file_name):
            raise FileNotFoundError(f"Rigid body restart file not found: {file_name}")

        with open(file_name, 'rb') as f:
            state = pickle.load(f)

        self.v_trans_prev = state['v_trans_prev']
        self.v_trans = state['v_trans']
        self.omega_prev = state['omega_prev']
        self.omega = state['omega']
        self.a_trans_prev = state['a_trans_prev']
        self.a_trans = state['a_trans']
        self.a_rot_prev = state['a_rot_prev']
        self.a_rot = state['a_rot']
        self.com = state['com']
        self.com_prev = state['com_prev']
        self.prev_patches = state['prev_patches']
        self.new_patches = state['new_patches']
        self.force_pr_it = state['force_pr_it']
        self.moment_pr_it = state['moment_pr_it']

        if self.restart:
            self.orientation_prev = state['orientation_prev']
            self.orientation = state['orientation']

        print('Info: Rigid body restart data successfully loaded.')
        return state.get('interface_rb', None)

    def _calc_aitken(self, q_raw, q_prev, res_prev, relax_factor, rb_iter):
        """
        Generic Aitken Δ² relaxation factor update.
        Suitable for forces, moments, accelerations, or any vector residual.
        """
        # Startup Stabilization
        if rb_iter <= 2:
            return relax_factor

        res = q_raw - q_prev
        res_diff = res - res_prev
        denom = np.dot(res_diff, res_diff)

        raw_relax = relax_factor
        if denom > 1e-14:
            num = np.dot(res_prev, res_diff)
            raw_relax = - relax_factor * num / denom

        # Smoothing (Slew Rate Limiter): Limit change to avoid shocks
        max_change = 0.2
        change = raw_relax - relax_factor
        change = np.clip(change, -max_change, max_change)

        new_relax_factor = relax_factor + change

        # Safety Clamping
        return np.clip(new_relax_factor, 0.01, 0.99)

    def _calculate_h_min(self, r_interface):
        """
        Groups contact nodes into spatially distinct 'Patches' using KD-Trees and a Flood Fill approach.
        Returns a list of patches, where each patch is a list of candidate dictionaries.
        """
        h_used = None
        h_min = None
        candidates = []

        # 1. KDTree Query & Danger Zone Detection
        for wall_name, wall_tree in self.gap_trees.items():

            # Find nearest distance from any interface node to the wall
            dists, wall_ids = wall_tree.query(r_interface, k=1, workers=-1)
            h = np.min(dists)
            h_min = min(h_min, h) if h_min is not None else h

            # Save h_used for Fictitious Mass (specific to fm_wall)
            if wall_name == self.fm_wall:
                h_used = h

            # Identify nodes inside the Danger Zone
            danger_mask = dists < self.h_ul
            danger_indices = np.where(danger_mask)[0]

            if len(danger_indices) > 0:
                close_dists = dists[danger_indices]
                close_ids_wall = wall_ids[danger_indices]
                close_p_itf = r_interface[danger_indices]
                close_p_wall = wall_tree.data[close_ids_wall]

                # Calculate normals for these points
                diff_vecs = close_p_itf - close_p_wall
                norms = np.linalg.norm(diff_vecs, axis=1)
                norms[norms < 1e-12] = 1.0
                close_normals = diff_vecs / norms[:, None]

                # Handle 2D -> 3D conversion for normals/points if needed
                if self.dimensions == 2:
                    z_col = np.zeros((len(danger_indices), 1))
                    close_normals = np.hstack((close_normals, z_col))
                    close_p_itf = np.hstack((close_p_itf, z_col))

                for j in range(len(danger_indices)):
                    candidates.append({
                        'h': close_dists[j],
                        'normal': close_normals[j],
                        'point': close_p_itf[j],
                        'id': danger_indices[j]
                    })

        # 2. Clustering (True Chaining / Flood Fill)
        candidates.sort(key=lambda x: x['h'])  # Deepest nodes first

        self.contact_patches = []
        processed_indices = set()
        sep_tol = 2 * self.h_ul

        for cand in candidates:
            if cand['id'] in processed_indices:
                continue

            current_patch = [cand]
            processed_indices.add(cand['id'])
            search_queue = [cand]

            # Flood fill loop to link connected neighbors
            while len(search_queue) > 0:
                expansion_node = search_queue.pop(0)
                expansion_point = expansion_node['point']

                for potential_neighbor in candidates:
                    if potential_neighbor['id'] in processed_indices:
                        continue

                    dist = np.linalg.norm(potential_neighbor['point'] - expansion_point)

                    if dist < sep_tol:
                        processed_indices.add(potential_neighbor['id'])
                        current_patch.append(potential_neighbor)
                        search_queue.append(potential_neighbor)

            self.contact_patches.append(current_patch)

        if len(self.contact_patches) > 0:
            print(f"CONTACT: Found {len(self.contact_patches)} patch(es).")
            for idx, patch in enumerate(self.contact_patches):
                print(f"  Patch {idx}: {len(patch)} nodes")

        return h_used, h_min

    def _calc_contact_force_and_moment(self):
        """
        Calculates repulsive force using Normalized Weighted Average (NWA)
        to smooth transitions between individual contact nodes.
        """
        total_force = np.zeros(3)
        total_moment = np.zeros(3)

        # Define Stiffness (k) and Damping (c) for harmonic oscillator response
        mass = self.k_mass * self.solid_density * self.volume
        resolution_steps = 10
        k_stiff = mass * (np.pi / (resolution_steps * self.delta_t)) ** 2
        c_damp = self.damping_ratio * 2 * np.sqrt(mass * k_stiff)

        available_prev = list(self.prev_patches)
        self.new_patches = []
        match_tolerance = 2.0 * self.h_ul

        for patch in self.contact_patches:
            # Accumulate Weighted Averages using linear kernel
            w_sum = 0.0
            weighted_normal = np.zeros(3)
            weighted_point = np.zeros(3)
            weighted_h = 0.0

            for node in patch:
                h_i = node['h']
                w_i = max(0.0, self.h_ul - h_i) / self.h_ul
                w_sum += w_i
                weighted_normal += w_i * node['normal']
                weighted_point += w_i * node['point']
                weighted_h += w_i * h_i

            if w_sum <= 1e-12:
                continue

            # Normalize values
            n_eff = weighted_normal / w_sum
            n_eff = n_eff / np.linalg.norm(n_eff)  # Re-normalize to unit vector
            p_eff = weighted_point / w_sum
            h_eff = weighted_h / w_sum

            # Velocity Tracking: find closest previous patch
            v_eff = 0.0
            if available_prev:
                best_dist = float('inf')
                best_index = -1

                for i, prev in enumerate(available_prev):
                    dist = np.linalg.norm(p_eff - prev['p_eff'])
                    if dist < best_dist:
                        best_dist = dist
                        best_index = i

                if best_dist < match_tolerance and best_index != -1:
                    matched_patch = available_prev.pop(best_index)
                    v_eff = (h_eff - matched_patch['h_eff']) / self.delta_t

            self.new_patches.append({'p_eff': p_eff, 'h_eff': h_eff})

            # Calculate Forces on this Effective Contact
            penetration = max(0.0, self.h_ul - h_eff)
            f_spring_mag = k_stiff * penetration  # Spring Force (Linear)
            print(f'Spring force = {f_spring_mag} N')

            f_damp_mag = -c_damp * 2 * (penetration / (self.h_ul - self.h_ll)) * v_eff  # Hunt-Crossley Damping
            print(f'v_eff = {v_eff * 1e6} µm/s')
            print(f'Damper force = {f_damp_mag} N')

            # Combine forces, disallow attractive forces
            f_mag = max(0.0, f_spring_mag + f_damp_mag)
            f_vec = f_mag * n_eff

            total_force += f_vec
            r_vec = p_eff - self.com
            total_moment += np.cross(r_vec, f_vec)

        return total_force, total_moment


# --- Module Helper Functions (Quaternion Math) ---
def quat_multiply(q1, q2):
    """Multiplies two quaternions."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    return np.array([w, x, y, z])


def quat_normalize(q):
    """Normalizes a quaternion to unit length to prevent numerical drift."""
    norm = np.linalg.norm(q)
    if norm == 0:
        return np.array([1.0, 0.0, 0.0, 0.0])  # Return identity quaternion
    return q / norm


def quat_from_angular_velocity(omega_vec, dt):
    """Creates a rotation quaternion from an angular velocity vector and timestep."""
    angle = np.linalg.norm(omega_vec) * dt
    if angle < 1e-12:  # Avoid division by zero for very small rotations
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = omega_vec / (angle / dt)
    half_angle = angle / 2.0
    w = np.cos(half_angle)
    x, y, z = axis * np.sin(half_angle)
    return np.array([w, x, y, z])


def quat_rotate_vector_array(q, v_array):
    """Efficiently rotates an array of 3D vectors by a unit quaternion q."""
    q_w, q_vec = q[0], q[1:]
    # Fast, vectorized formula for quaternion rotation
    return v_array + 2 * np.cross(q_vec, np.cross(q_vec, v_array) + q_w * v_array)


def skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])


def exp_map(x):
    """Exponential map to convert rotation vectors to rotation matrices."""
    if np.linalg.norm(x) < 1e-12:
        return np.identity(3)
    else:
        return np.identity(3) + (np.sin(np.linalg.norm(x)) / np.linalg.norm(x)) * skew(x) + (1 / 2) * (
                (np.sin(np.linalg.norm(x)) ** 2) / ((np.linalg.norm(x) / 2) ** 2)) * (skew(x) @ skew(x))