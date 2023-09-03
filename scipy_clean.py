import numpy as np
from scipy import stats

def load_point_cloud(filename):
    """Load a point cloud from an ASCII XYZ file."""
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y, z = map(float, line.strip().split())
            points.append([x, y, z])
    return np.array(points)

def save_point_cloud(filename, points):
    """Save a point cloud to an ASCII XYZ file."""
    with open(filename, 'w') as file:
        for point in points:
            file.write(f"{point[0]} {point[1]} {point[2]}\n")

def despike_point_cloud(points, threshold=3.0):
    """Despike a point cloud using the Z-score method."""
    z_scores = np.abs(stats.zscore(points[:, 2]))  # Compute Z-scores for the Z-coordinate
    filtered_points = points[z_scores < threshold]
    return filtered_points

# Example usage
input_file = 'input.xyz'
output_file = 'output.xyz'
threshold = 3.0

# Load the point cloud
point_cloud = load_point_cloud(input_file)

# Despike the point cloud
filtered_points = despike_point_cloud(point_cloud, threshold)

# Save the filtered point cloud
save_point_cloud(output_file, filtered_points)


