"""
This file contains the module load commands for the solvers, for different machines in form of python-dict. This
is used for setting up the environment variables required by the solvers. During the simulation, the solver wrapper
uses the environment variables in all the spawned processes, used for executing the solver-specific terminal commands
in the various stages of the simulation.
"""

machine_name = 'ugent_cluster_RHEL8'

solver_load_cmd_dict = {
    'ugent_cluster_RHEL8': {
        'fluent.v2023R1': 'ml ANSYS_CFD/2023R1',
        'cht_fluent.v2023R1': 'ml ANSYS_CFD/2023R1',
        'pc_fluent.v2023R1': 'ml ANSYS_CFD/2023R1',
        'fluent.v2024R1': 'ml ANSYS_CFD/2024R1',
        'cht_fluent.v2024R1': 'ml ANSYS_CFD/2024R1',
        'pc_fluent.v2024R1': 'ml ANSYS_CFD/2024R1',
        'fluent.v2024R2': 'ml ANSYS_CFD/2024R2',
        'kratos_structure.v94': 'ml Kratos/9.4.5-Anaconda3-2023.09-Python-3.11',
        'openfoam.v10': 'ml OpenFOAM/10-foss-2023a && source $FOAM_BASH',
        'openfoam.v11': 'ml OpenFOAM/11-foss-2023a && source $FOAM_BASH',
        'abaqus.v2023': 'ml intel/2022b && ml ABAQUS/2023 '
                        '&& export LM_LICENSE_FILE=@ir03lic1.ugent.be:27000@ea11serv03.private.ugent.be',
        'abaqus.v2024': 'ml intel/2022b && ml ABAQUS/2024 '
                        '&& export LM_LICENSE_FILE=@ir03lic1.ugent.be:27000@ea11serv03.private.ugent.be',
    },
    'ugent_hpc': {
        'fluent.v2023R1': 'ml FLUENT/2023R1 ; ml iimpi/2023a ; unset SLURM_GTIDS '
                          '&& export ANSYSLI_SERVERS=2325@ir03lic1.ugent.be '
                          '&& export ANSYSLMD_LICENSE_FILE=1055@ir03lic1.ugent.be',
        'fluent.v2024R1': 'ml FLUENT/2024R1 ; unset SLURM_GTIDS '
                          '&& export ANSYSLI_SERVERS=2325@ir03lic1.ugent.be '
                          '&& export ANSYSLMD_LICENSE_FILE=1055@ir03lic1.ugent.be',
        'fluent.v2024R2': 'ml FLUENT/2024R2 ; unset SLURM_GTIDS '
                          '&& export ANSYSLI_SERVERS=2325@ir03lic1.ugent.be '
                          '&& export ANSYSLMD_LICENSE_FILE=1055@ir03lic1.ugent.be',
        'abaqus.v2023': 'ml intel/2022b ; ml ABAQUS/2023 '
                        '; export LM_LICENSE_FILE=@ir03lic1.ugent.be:27000@ea11serv03.private.ugent.be',
        'abaqus.v2024': 'ml intel/2022b ; ml ABAQUS/2024-hotfix-2405 '
                        '; export LM_LICENSE_FILE=@ir03lic1.ugent.be:27000@ea11serv03.private.ugent.be',
        'kratos_structure.v94': 'export OMP_PROC_BIND=TRUE',  # needs to be FALSE in job script, but TRUE for Kratos
        'openfoam.v10': 'ml OpenFOAM/10-foss-2023a ; source $FOAM_BASH',
        'openfoam.v11': 'ml OpenFOAM/11-foss-2023a ; source $FOAM_BASH'
    },
    'hortense': {
        'fluent.v2023R1': 'ml FLUENT/2023R1 ; ml intel/2021a '
                          '; export ANSYSLI_SERVERS=2325@ir03lic1.ugent.be '
                          '&& export ANSYSLMD_LICENSE_FILE=1055@ir03lic1.ugent.be',
        'kratos_structure.v94': 'export OMP_PROC_BIND=TRUE',  # needs to be FALSE in job script, but TRUE for Kratos
    }
}


def get_solver_cmd(solver_name):
    return solver_load_cmd_dict[machine_name][solver_name]
