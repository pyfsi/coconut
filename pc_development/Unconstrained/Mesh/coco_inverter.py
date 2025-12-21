import numpy as np

# 1. Load all data (x, y, and id)
data = np.loadtxt('nodes_timestep0_thread11.dat', skiprows=1)

# 2. Find indices of unique IDs (looking only at column 2)
# return_index=True gives us the index of the first occurrence of each unique ID
_, indices = np.unique(data[:, 2], return_index=True)

# 3. Sort indices to preserve the original file order
indices.sort()

# 4. Filter the data by these indices, then slice off the ID column (keep only cols 0 and 1)
final_data = data[indices][:, :2]

# 5. Save to CSV
np.savetxt('interface_nodes.csv', final_data, delimiter=',',
           header='x-coordinate,y-coordinate', comments='', fmt='%.18e')