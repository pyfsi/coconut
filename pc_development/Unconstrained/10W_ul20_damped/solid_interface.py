import pandas as pd

# Define filenames
input_filename = 'CFD_2/node_coords_inner'
output_filename = 'CFD_1/interface_nodes.csv'

# 1. Read the input file
# skipinitialspace=True handles the spaces after commas in your file
df = pd.read_csv(input_filename, sep=r'\s+', skipinitialspace=True)

# 2. Extract specific columns
# Your file structure: nodenumber (col 0), x (col 1), y (col 2), y (col 3), x (col 4)
# We select columns at index 1 and 2
coords_df = df.iloc[:, [1, 2]]

# Set the proper column names for the output
coords_df.columns = ['x-coordinate', 'y-coordinate']

# 3. Remove duplicate nodes
unique_df = coords_df.drop_duplicates()

# 4. Save to CSV
# float_format='%.18e' ensures high precision scientific notation matches your example
unique_df.to_csv(output_filename, index=False, float_format='%.18e')

print(f"Successfully created {output_filename} with {len(unique_df)} unique nodes.")