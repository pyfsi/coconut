#!/bin/bash

# 1. Load Modules
ml -ANSYS_CFD
ml GAMBIT/2.4.6

# 2. Run Identifier Step
# This creates identified_edges.dbs and identified_edges.trn
gambit -inp identifyer.jou > setup_gambit_id.log 2>&1

# 3. Run Python to Generate remesh.jou
# Ensure python is available (usually standard, or load a module)
python3 generate_remesh.py

# 4. Run the Generated Remesh Script
gambit -inp remesh.jou > setup_gambit_remesh.log 2>&1

# 5. Cleanup
ml -GAMBIT