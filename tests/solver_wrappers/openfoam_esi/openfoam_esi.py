from coconut.tools import create_instance, get_solver_env
import coconut.coupling_components.solver_wrappers.openfoam_esi.openfoam_esi_io as of_esi_io

import unittest
import numpy as np
import os
import multiprocessing
import shutil
import re
import json
import glob
from subprocess import check_call, DEVNULL
from contextlib import contextmanager

# Tolerances
LOW_TOL = 1e-3
MEDIUM_TOL = 1e-6
HIGH_TOL = 1e-9
VERY_HIGH_TOL = 1e-12


class TestSolverWrapperOpenFoamESI(unittest.TestCase):
    version = None # Openfoam version without dot

    # ==================== overloaded Methods ====================
    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = os.path.join(dir_name, f"test_v{cls.version}/tube3d/parameters.json")
        cls.working_dir = os.path.join(dir_name, f"test_v{cls.version}/tube3d/CFD")

        shutil.rmtree(cls.working_dir, ignore_errors=True)
        src = os.path.join(dir_name, f"test_v{cls.version}/tube3d/setup_openfoam_esi")
        shutil.copytree(src, cls.working_dir)

    def setUp(self):
        with open(self.file_name, "r") as parameter_file:
            self.parameters = json.load(parameter_file)

        # set working directory
        self.parameters["settings"]["working_directory"] = os.path.relpath(self.working_dir)
        settings = self.parameters["settings"]

        # set model parts
        self.mp_name_in = settings["interface_input"][0]["model_part"]
        self.mp_name_out = settings["interface_output"][0]["model_part"]

        self.folder_path = os.path.join(os.getcwd(), settings["working_directory"])
        self.delta_t = settings["delta_t"]
        self.t_prec = settings["time_precision"]
        self.max_cores = min(4, multiprocessing.cpu_count())

        solver_name = self.parameters["type"].replace("solver_wrappers.","")
        self.env = get_solver_env(solver_name, self.folder_path)
        self._clean_case()
        self._setup_case()

    def tearDown(self):
        self._clean_case()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.working_dir)

    # ==================== static methods ====================
    @staticmethod
    def calc_displacement(node_pos):
        '''Calculate displacement from a prescribed radial displacement'''
        x,y,z = node_pos[:,0], node_pos[:,1], node_pos[:,2]

        dx = np.zeros(x.shape)
        dr = 0.0001 * np.sin(2.0 * np.pi / 0.05 * x)
        theta = np.arctan2(z, y)
        dy = dr * np.cos(theta)
        dz = dr * np.sin(theta)
        return np.column_stack((dx, dy, dz))

    @staticmethod
    def rm_time_dirs(solver):
        '''Remove all time directories with pattern 0.*'''
        time_dirs = glob.glob(os.path.join(solver.working_directory, "0.*"))
        for time_dir in time_dirs:
            shutil.rmtree(time_dir)

    # ==================== helper functions ====================
    def _clean_case(self):
            check_call("sh " + os.path.join(self.folder_path, "Allclean"), shell=True, env=self.env)

    def _setup_case(self):
        check_call("sh " + os.path.join(self.folder_path, "Allprep"), shell=True, env=self.env)

    def _set_cores(self, cores):
        settings = self.parameters["settings"]
        if cores == 1:
            settings["parallel"] = False
        else:
            settings["parallel"] = True
            decomposepardit_path = os.path.join(self.folder_path, "system", "decomposeParDict")
            with open(decomposepardit_path, "r") as f:
                old_dict= f.read()
            new_dict = re.sub(r"numberOfSubdomains[\s\n]"+of_esi_io.int_pattern,
                              f"numberOfSubdomains {cores}",
                              old_dict)
            with open(decomposepardit_path, "w") as f:
                f.write(new_dict)

    def _set_tolerances(self, p_atol, p_rtol, U_atol, U_rtol):
        solution_file_name = os.path.join(self.folder_path, 'system', 'fvSolution')
        new_tol = f'residualControl\n' \
                  f'    {{\n' \
                  f'    p\n' \
                  f'        {{\n' \
                  f'            tolerance   {p_atol:g};\n' \
                  f'            relTol      {p_rtol:g};\n' \
                  f'        }}\n' \
                  f'    U\n' \
                  f'        {{\n' \
                  f'            tolerance   {U_atol:g};\n' \
                  f'            relTol      {U_rtol:g};\n' \
                  f'        }}\n' \
                  f'    }}'
        with open(solution_file_name, 'r') as f:
            old_dict = f.read()
        new_dict = re.sub(of_esi_io.get_nested_dict(old_dict, 'residualControl'), new_tol, old_dict)
        with open(solution_file_name, 'w') as f:
            f.write(new_dict)

    @contextmanager
    def run_solver(self, solver):
        solver.initialize()
        yield
        solver.finalize()

    @contextmanager
    def run_solution_step(self, solver):
        solver.initialize_solution_step()
        yield
        solver.finalize_solution_step()
        solver.output_solution_step()

    def iterate_solver(self, solver, interface_input=None, with_output=False):
        out = solver.solve_solution_step(interface_input)
        if with_output:
            solver.finalize_solution_step()
            solver.output_solution_step()
        return out

    # ==================== unittests ====================
    # test if order of nodes vary with different decomposition
    def test_model_part_nodes_with_different_cores(self):
        node_pos = {"x0":[], "y0":[], "z0":[], "id":[]}

        n_cores = [1, self.max_cores]
        for cores in n_cores:
            self._set_cores(cores)
            solver = create_instance(self.parameters)
            with self.run_solver(solver):
                pass

            # store model part data to nodes dict
            model_part = solver.model.get_model_part(self.mp_name_in)
            for key in node_pos.keys():
                node_pos[key].append(getattr(model_part, key))

            # remove solution dirs
            self.rm_time_dirs(solver)

        # compare ModelParts of both solvers
        for key in node_pos.keys():
            np.testing.assert_array_almost_equal(node_pos[key][0], node_pos[key][1])

    # test node displacement
    def test_displacement_on_nodes(self):
        n_cores = [1, self.max_cores]
        for cores in n_cores:
            self._setup_case()
            self._set_cores(cores)

            solver = create_instance(self.parameters)
            with self.run_solver(solver):
                interface_input = solver.get_interface_input()

                # set displacements
                model_part = interface_input.get_model_part(self.mp_name_in)
                node_pos0 = np.column_stack((
                    getattr(model_part, "x0"),
                    getattr(model_part, "y0"),
                    getattr(model_part, "z0")
                ))
                displacement = self.calc_displacement(node_pos0)
                interface_input.set_variable_data(self.mp_name_in, "displacement", displacement)

                # calculate reference node position
                node_pos_ref = node_pos0 + displacement

                # do one iteration
                solver.initialize_solution_step()
                self.iterate_solver(solver, interface_input, with_output=True)

            if cores > 1:
                reconstruct_cmd = f"cd {self.folder_path} && reconstructPar -latestTime -no-fields"
                check_call(reconstruct_cmd, shell=True, stdout=DEVNULL, env=self.env)

            time_folder_name = f"{self.delta_t:.{solver.time_precision}f}"
            _, node_pos = of_esi_io.get_boundary_points(solver.working_directory,
                                                           time_folder_name,
                                                           "mantle")

            np.testing.assert_allclose(node_pos, node_pos_ref, rtol=VERY_HIGH_TOL)

            self.rm_time_dirs(solver)

    # test pressure and traction by varying partition
    def test_stress_on_nodes_parallel(self):
        interface_output_list = []
        n_cores = [1, self.max_cores]
        for cores in n_cores:
            self._setup_case()
            self._set_cores(cores)

            solver = create_instance(self.parameters)
            with self.run_solver(solver):
                interface_input = solver.get_interface_input()

                # set displacement
                model_part = interface_input.get_model_part(self.mp_name_in)
                node_pos0 = np.column_stack((
                                getattr(model_part, "x0"),
                                getattr(model_part, "y0"),
                                getattr(model_part, "z0")
                            ))
                displacement = self.calc_displacement(node_pos0)
                interface_input.set_variable_data(self.mp_name_in, "displacement", displacement)

                # do one solver iteration
                solver.initialize_solution_step()
                out = self.iterate_solver(solver, interface_input, with_output=True)
                interface_output_list.append(out.get_interface_data())

            self.rm_time_dirs(solver)

        output_ref = interface_output_list[0]
        max_value = np.max(np.abs(output_ref))
        for output in interface_output_list[1:]:
            np.testing.assert_allclose(output/max_value, output_ref/max_value,
                                       atol=LOW_TOL, rtol=0)

    # test pressure and traction by varying node position
    def test_stress_on_nodes(self):
        self._set_cores(1)

        # define displacement scaling
        displacement_scaling = [1.0, 0.0, 1.0]
        # storage for pressure and traction
        stress = {"pressure":[None]*len(displacement_scaling), "traction":[None]*len(displacement_scaling)}

        # run solver
        solver = create_instance(self.parameters)
        with self.run_solver(solver):
            # get interface input
            interface_input = solver.get_interface_input()

            # set displacement
            model_part = interface_input.get_model_part(self.mp_name_in)
            node_pos0 = np.column_stack((
                            getattr(model_part, "x0"),
                            getattr(model_part, "y0"),
                            getattr(model_part, "z0")
                        ))
            displacement = self.calc_displacement(node_pos0)

            # do solver iteration three times
            solver.initialize_solution_step()
            for i, scale_i in enumerate(displacement_scaling):
                interface_input.set_variable_data(self.mp_name_in, "displacement", scale_i*displacement)
                interface_output = self.iterate_solver(solver, interface_input, with_output=False)
                stress["pressure"][i] = (interface_output.get_variable_data(self.mp_name_out,
                                                                     "pressure"))
                stress["traction"][i] = (interface_output.get_variable_data(self.mp_name_out,
                                                                     "traction"))

            solver.finalize_solution_step()
            solver.output_solution_step()

        # calc stress amplitudes
        avg_pressure = 0.5 * (stress["pressure"][0]+stress["pressure"][2])
        avg_traction = 0.5 * (stress["traction"][0]+stress["traction"][2])
        pressure_amplitude = 0.5 * (np.max(avg_pressure) - np.min(avg_pressure))
        traction_amplitude = 0.5 * (np.max(avg_traction) - np.min(avg_traction))

        # check stress tensors on same positions
        np.testing.assert_allclose(stress["pressure"][0]/pressure_amplitude, stress["pressure"][2]/pressure_amplitude,
                                   atol=LOW_TOL, rtol=0)
        np.testing.assert_allclose(stress["traction"][0]/traction_amplitude, stress["traction"][2]/traction_amplitude,
                                           atol=LOW_TOL, rtol=0)

        # check stress tensors on different positions
        p01 = np.linalg.norm(stress["pressure"][0][:] - stress["pressure"][1][:])
        p02 = np.linalg.norm(stress["pressure"][0][:] - stress["pressure"][2][:])
        t01 = np.linalg.norm(stress["traction"][0][:,:] - stress["traction"][1][:,:])
        t02 = np.linalg.norm(stress["traction"][0][:,:] - stress["traction"][2][:,:])
        self.assertTrue(p02/p01<LOW_TOL)
        self.assertTrue(t02/t01<LOW_TOL)

    # test coupling results in case of restart
    def test_restart(self):
        self._set_cores(4)
        self.parameters["settings"]["cores"] = 4
        self.parameters["settings"]["save_restart"] = 2
        n_time_steps = 3

        # pressure and traction
        results_1 = {"input": None, "output":None, "node_pos":None, "pressure":None, "traction":None}
        results_2 = {"input": None, "output":None, "node_pos":None, "pressure":None, "traction":None}

        # run solver from the start
        solver = create_instance(self.parameters)
        with self.run_solver(solver):
            interface_input = solver.get_interface_input()

            # set displacements
            model_part = interface_input.get_model_part(self.mp_name_in)
            node_pos0 = np.column_stack((
                getattr(model_part, "x0"),
                getattr(model_part, "y0"),
                getattr(model_part, "z0")
            ))
            displacement = self.calc_displacement(node_pos0)
            interface_input.set_variable_data(self.mp_name_in, "displacement", displacement)

            # iterate solver
            for i in range(n_time_steps):
                solver.initialize_solution_step()
                interface_input.set_variable_data(self.mp_name_in, "displacement", i*displacement)
                self.iterate_solver(solver, interface_input, with_output=True)

            results_1["input"] = solver.get_interface_input()
            results_1["output"] = solver.get_interface_output()

            check_call(f'reconstructPar -latestTime -no-fields', shell=True, cwd=self.folder_path, stdout=DEVNULL,
                        env=self.env)
            _, results_1["node_pos"] = of_esi_io.get_boundary_points(solver.working_directory,
                                                    f'{n_time_steps * self.delta_t:.{solver.time_precision}f}', 'mantle')

        # get data without restart
        interface_output = solver.get_interface_output()
        results_1["pressure"] = interface_output.get_variable_data(self.mp_name_out, "pressure")
        results_1["traction"] = interface_output.get_variable_data(self.mp_name_out, "traction")

        # run solver in restart mode
        t_restart = 2
        self.parameters["settings"]["timestep_start"] = t_restart
        solver = create_instance(self.parameters)
        with self.run_solver(solver):
            interface_input = solver.get_interface_input()

            for i in range(t_restart, n_time_steps):
                solver.initialize_solution_step()
                interface_input.set_variable_data(self.mp_name_in, "displacement", i*displacement)
                self.iterate_solver(solver, interface_input, with_output=True)

            results_2["input"] = solver.get_interface_input()
            results_2["output"] = solver.get_interface_output()

            check_call(f'reconstructPar -latestTime -no-fields', shell=True, cwd=self.folder_path, stdout=DEVNULL,
                                    env=self.env)
            _, results_2["node_pos"] = of_esi_io.get_boundary_points(solver.working_directory,
                                                    f'{n_time_steps * self.delta_t:.{solver.time_precision}f}', 'mantle')

        # get data from solver with restart
        interface_output = solver.get_interface_output()
        results_2["pressure"] = interface_output.get_variable_data(self.mp_name_out, "pressure")
        results_2["traction"] = interface_output.get_variable_data(self.mp_name_out, "traction")

         # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(results_1["input"].has_same_model_parts(results_2["input"]))
        self.assertTrue(results_1["output"].has_same_model_parts(results_2["output"]))

        # check if deformed coordinates are equal
        np.testing.assert_allclose(results_1["node_pos"], results_2["node_pos"], rtol=VERY_HIGH_TOL)

        # check if pressure and traction are equal
        np.testing.assert_allclose(results_1["pressure"], results_2["pressure"], rtol=HIGH_TOL)
        np.testing.assert_allclose(results_1["traction"], results_2["traction"], rtol=HIGH_TOL)

    # test coupling convergence considering different applied displacements
    def test_check_coupling(self):
        # apply different displacement between 1st and 2nd iter and same displacement between 2nd and 3rd iteration
        displacement_scaling = [1.0, 2.0, 2.0]

        # change tolerances in order to make coupling convergence work
        self._set_tolerances(1e-3, 0.0, 1e-3, 0.0)

        n_cores = [1]#, self.max_cores]
        for cores in n_cores:
            self._set_cores(cores)
            self.parameters["settings"]["cores"] = cores

            solver = create_instance(self.parameters)
            solver.check_coupling_convergence = True
            with self.run_solver(solver):
                interface_input = solver.get_interface_input()

                # set displacements
                model_part = interface_input.get_model_part(self.mp_name_in)
                node_pos0 = np.column_stack((
                    getattr(model_part, "x0"),
                    getattr(model_part, "y0"),
                    getattr(model_part, "z0")
                ))
                displacement = self.calc_displacement(node_pos0)

                with self.run_solution_step(solver):
                    for i, scale_i in enumerate(displacement_scaling):
                        interface_input.set_variable_data(self.mp_name_in, "displacement", scale_i*displacement)
                        solver.solve_solution_step(interface_input)

                        if i==2:
                            self.assertTrue(solver.coupling_convergence)
                        else:
                            self.assertFalse(solver.coupling_convergence)

            self.rm_time_dirs(solver)
