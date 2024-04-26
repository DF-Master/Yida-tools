import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# Read the data
file_path = r'C:\Users\jiang\OneDrive\Research\tc\Data\2024\04\justify结果_hxn_重算3-2.xlsx'
df_xyz = pd.read_excel(
    file_path,
    usecols='AZ:BB',
)[:-4]
df_XL = pd.read_excel(file_path, usecols='S:AR')[:-4]
# Remove the title row
df_xyz.columns = ['x', 'y', 'z']

# Check if there are NaN values in the dataframe
print(df_XL.isna().any())

# Counts the number of NaNs in each column
print(df_XL.isna().sum())

# Find rows with NaN values
nan_rows = df_XL[df_XL.isna().any(axis=1)]

# Display the rows with NaN values
print(nan_rows)

# Using sklearn kmeans to do clustering
kmeans = KMeans(n_clusters=10)
kmeans.fit(df_XL)

# Predicting clusters on our data
labels = kmeans.predict(df_XL)

# Adding labels to the original dataframe
df_XL['labels'] = labels
df_output = pd.concat([df_xyz, df_XL], axis=1)
with pd.ExcelWriter(file_path, engine='openpyxl', mode='a') as writer:
    df_output.to_excel(writer, sheet_name='XL_labels_10clusters_allXL_2')

# 3D Plotting
fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection="3d")

ax.scatter3D(
    df_xyz['x'],
    df_xyz['y'],
    df_xyz['z'],
    c=df_XL['labels'],
    cmap='viridis',
    s=5,
)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()