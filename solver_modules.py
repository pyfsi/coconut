"""
This file contains the module load commands for the solvers, for different machines in form of python-dict. This
is used for setting up the environment variables required by the solvers. During the simulation, the solver wrapper
uses the environment variables in all the spawned processes, used for executing the solver-specific terminal commands
in the various stages of the simulation.
"""

machine_name = 'ugent_cluster_SL6.3'

solver_load_cmd_dict = {
    'ugent_cluster_SL6.3': {
        'fluent.v2019R1': 'ml ANSYS_CFD/2019R1',
        'fluent.v2019R2': 'ml ANSYS_CFD/2019R2',
        'fluent.v2019R3': 'ml ANSYS_CFD/2019R3',
        'fluent.v2020R1': 'ml ANSYS_CFD/2020R1',
        'abaqus.v614': 'ml intel/2018a && ml ABAQUS/6.14',
        'kratos.structural_mechanics_application.v60': 'ml Kratos/6.0-foss-2018a-Python-3.6.4',
        'openfoam.v41': 'ml OpenFOAM/4.1 && source $FOAM_BASH'
    },
    'ugent_cluster_CO7': {
        'fluent.v2019R1': 'ml ANSYS_CFD/2019R1',
        'fluent.v2019R3': 'ml ANSYS_CFD/2019R3',
        'fluent.v2020R1': 'ml ANSYS_CFD/2020R1',
        'abaqus.v614': 'ml intel/2018a && ml ABAQUS/6.14',
        'kratos.structural_mechanics_application.v60': 'ml Kratos/6.0-foss-2018a-Python-3.6.4',
        'openfoam.v41': 'ml OpenFOAM/4.1-foss-2016b && source $FOAM_BASH',
        'openfoam.v8': 'ml OpenFOAM/8-foss-2020b && source $FOAM_BASH'
    },
    'ugent_hpc': {
        'fluent.v2019R1': 'ml ANSYS_CFD/2019R1',
        'fluent.v2019R2': 'ml ANSYS_CFD/2019R2',
        'fluent.v2019R3': 'ml FLUENT/2019R3',
        'abaqus.v614': 'ml intel/2018a && ml ABAQUS/6.14.1-linux-x86_64 && unset SLURM_GTIDS',
        'kratos.structural_mechanics_application.v60': 'ml Kratos/6.0-foss-2018a-Python-3.6.4',
    }
}


def get_solver_cmd(solver_name):
    return solver_load_cmd_dict[machine_name][solver_name]
