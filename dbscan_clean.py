import numpy as np
from sklearn.cluster import DBSCAN

def load_point_cloud(filename):
    """Load a point cloud from an ASCII XYZ file."""
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y, z = map(float, line.strip().split(","))
            points.append([x, y, z])
    return np.array(points)

def save_point_cloud(filename, points):
    """Save a point cloud to an ASCII XYZ file."""
    with open(filename, 'w') as file:
        for point in points:
            file.write(f"{point[0]} {point[1]} {point[2]}\n")

def despike_point_cloud(points, eps, min_samples):
    """Despike a point cloud using DBSCAN."""
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(points)
    labels = clustering.labels_
    filtered_points = points[labels != -1]
    return filtered_points

# Example usage
input_file = 'input.xyz'
output_file = 'output.xyz'
eps = 0.1  # DBSCAN epsilon parameter
min_samples = 5  # DBSCAN minimum number of points

# Load the point cloud
point_cloud = load_point_cloud(input_file)

# Despike the point cloud
filtered_points = despike_point_cloud(point_cloud, eps, min_samples)

# Save the filtered point cloud
save_point_cloud(output_file, filtered_points)
