import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Simulate depth measurements
num_measurements = 100
depth_values = np.concatenate([
    np.random.normal(10, 1, int(num_measurements * 0.8)),  # Majority of data
    np.random.normal(50, 10, int(num_measurements * 0.2))  # Outliers
])

# Calculate Z-scores
z_scores = np.abs(stats.zscore(depth_values))

# Set a Z-score threshold for outlier detection
z_score_threshold = 2.5

# Identify outliers
outliers = np.where(z_scores > z_score_threshold)[0]

# Remove outliers
cleaned_depth_values = np.delete(depth_values, outliers)

# Plot the original and cleaned data
plt.figure(figsize=(10, 6))
plt.scatter(range(num_measurements), depth_values, color='b', label='Original Data')
plt.scatter(outliers, depth_values[outliers], color='r', label='Outliers')
plt.plot(cleaned_depth_values, 'g', label='Cleaned Data')
plt.xlabel('Measurement Index')
plt.ylabel('Depth')
plt.title('Z-score Outlier Cleaning')
plt.legend()
plt.grid(True)
plt.show()
