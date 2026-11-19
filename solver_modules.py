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
        'fluent.v2024R1': 'ml ANSYS_CFD/2024R1',
        'fluent.v2024R2': 'ml ANSYS_CFD/2024R2',
        'fluent.v2025R2': 'ml ANSYS_CFD/2025R2',
        'kratos_structure.v94': 'ml Kratos/9.4.5-Anaconda3-2023.09-Python-3.11',
        'openfoam_org.v10': 'ml OpenFOAM/10-foss-2023a && source $FOAM_BASH',
        'openfoam_org.v11': 'ml OpenFOAM/11-foss-2023a && source $FOAM_BASH',
        'openfoam_esi.v2312': 'ml OpenFOAM/v2312-foss-2023a && source $FOAM_BASH',
        'abaqus.v2023': 'ml intel/2022b && ml ABAQUS/2023 && [ -f ~/abaqus_v6.env ] || cat > ~/abaqus_v6.env << EOF\n'
                        'import os\nlicense_server_type=DSLS\ndsls_license_config=os.path.join(os.path.expanduser("~"),'
                        '"DSLicSrv.txt")\ndel os\nEOF\n cat > ~/DSLicSrv.txt << EOF\ncfdfile1.ugent.be:4085\n'
                        'ea11serv03.private.ugent.be:4085\nEOF\n echo ""',
        'abaqus.v2024': 'ml intel/2022b && ml ABAQUS/2024 && [ -f ~/abaqus_v6.env ] || cat > ~/abaqus_v6.env << EOF\n'
                        'import os\nlicense_server_type=DSLS\ndsls_license_config=os.path.join(os.path.expanduser("~"),'
                        '"DSLicSrv.txt")\ndel os\nEOF\n cat > ~/DSLicSrv.txt << EOF\ncfdfile1.ugent.be:4085\n'
                        'ea11serv03.private.ugent.be:4085\nEOF\n echo ""',
        'abaqus_cse.v2025': 'ml ABAQUS/2025 && [ -f ~/abaqus_v6.env ] || cat > ~/abaqus_v6.env << EOF\n'
                        'import os\nlicense_server_type=DSLS\ndsls_license_config=os.path.join(os.path.expanduser("~"),'
                        '"DSLicSrv.txt")\ndel os\nEOF\n echo ea11serv03.private.ugent.be:4085 > ~/DSLicSrv.txt',
        'abaqus_cse.v2026': 'ml ABAQUS/2026 && [ -f ~/abaqus_v6.env ] || cat > ~/abaqus_v6.env << EOF\n'
                        'import os\nlicense_server_type=DSLS\ndsls_license_config=os.path.join(os.path.expanduser("~"),'
                        '"DSLicSrv.txt")\ndel os\nEOF\n echo ea11serv03.private.ugent.be:4085 > ~/DSLicSrv.txt'
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
        'abaqus.v2024': 'ml intel/2024a; ml ABAQUS/2024-hotfix-2405; unset SLURM_GTIDS; [ -f ~/abaqus_v6.env ] || '
                        'cat > ~/abaqus_v6.env << EOF\nimport os\nlicense_server_type=DSLS\ndsls_license_config='
                        'os.path.join(os.path.expanduser("~"), "DSLicSrv.txt")\ndel os\nEOF\n'
                        'echo ea11serv03.private.ugent.be:4085 > ~/DSLicSrv.txt',
        'abaqus_cse.v2025': 'ml intel/2024a; ml ABAQUS/2025; unset SLURM_GTIDS; [ -f ~/abaqus_v6.env ] || '
                        'cat > ~/abaqus_v6.env << EOF\nimport os\nlicense_server_type=DSLS\ndsls_license_config='
                        'os.path.join(os.path.expanduser("~"), "DSLicSrv.txt")\ndel os\nEOF\n'
                        'echo ea11serv03.private.ugent.be:4085 > ~/DSLicSrv.txt',
        'kratos_structure.v94': 'export OMP_PROC_BIND=TRUE',  # needs to be FALSE in job script, but TRUE for Kratos
        'openfoam_org.v10': 'ml OpenFOAM/10-foss-2023a; source $FOAM_BASH',
        'openfoam_org.v11': 'ml OpenFOAM/11-foss-2023a; source $FOAM_BASH',
        'openfoam_esi.v2312': 'ml OpenFOAM/v2312-foss-2023a && source $FOAM_BASH',
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
