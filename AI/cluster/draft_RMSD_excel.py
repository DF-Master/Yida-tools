import pandas as pd
import numpy as np
from scipy.spatial.transform import Rotation as R

# Load the data from the Excel file
data = pd.read_excel(
    "C:\\Users\\jiang\\OneDrive\\Research\\tc\\Data\\2024\\04\\justify结果_hxn_重算3-2-bk.xlsx"
)

# Use iloc to access the columns by their integer locations
# Excel columns:      A  B  C  D ---> Z  AA AB AC ... AZ BA BB, BC, BD ...
# Python positions:   0  1  2  3 ---> 25 26 27 28 ... 51 52 53, 54, 55 ...
ref_points = data.iloc[:, [54, 55, 56]].values
mobile_points = data.iloc[:, [51, 52, 53]].values

# Compute the centroids of each point set
ref_centroid = np.mean(ref_points, axis=0)
mobile_centroid = np.mean(mobile_points, axis=0)

# Translate the points to bring the centroids to the origin
ref_points -= ref_centroid
mobile_points -= mobile_centroid

#Determine the rotation matrix that aligns the two point sets using Kabsch algorithm
rotation_matrix = R.from_rotvec(ref_points).as_matrix()

# Rotate the mobile points to align them with the reference points
mobile_points = np.dot(mobile_points, rotation_matrix)

# Compute the root-mean-square deviation (RMSD) between the aligned points
rmsd = np.sqrt(np.mean((ref_points - mobile_points)**2))

print(f"RMSD: {rmsd}")