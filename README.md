# ggoutlier.py
ggoutlier uses machine learning...

##What is an Outlier?
In bathymetry, an outlier typically refers to an isolated or anomalous depth measurement or feature on a seafloor depth map or chart. These outliers can be depths that are significantly different from the surrounding seafloor topography. Outliers might be caused by various factors such as errors in data collection, equipment malfunction, or unique geological features like wrecks, obstructions, seamounts or underwater volcanoes that stand out from the surrounding seabed. The scale of an outlier can be considerable. Identifying and understanding  outliers in bathymetric data is important for accurate navigation, scientific research, and ocean exploration. Separating a real feature from noise is a complex issue. The final decision comes down to the skill and experience of the Surveyor In Charge. GGOutlier efficiently analyse and highlight outliers for validation.

##Installation Requirements
* python 3.10 or newer
* pip install reportlab
* pip install pyproj
* pip install pyshp
* pip install scikit-learn

##GGOutlier Principles
GGOutlier is a tool developed by Guardian Geomatics to Quality Control processed a multibeam bathymetry surface, and validate a processed depth surface against a standard such as those published by IHO SP44 or HIPP. The principle is similar in methodology to a traditional review by a surveyor-in-charge (SIC) process in which the SIC would review the depth surface by identifying outliers relative to its nearest neighours, determine if the outlier is significant and if so flag it
forinvestigation.
GGOutlier primary purpose is to positively, rigorously identify each and every depth which is considered an outlier relative to the required total vertical uncertainty at that depth. This provides the SIC and client with full confidence that the quality of the depth surface meets the required specification and any remaining outliers are known, have been investigated and are considered features rather than noise.
Inputs are very simple. A depth surface (a floating point TIF file) and a IHO SP44 specification such as 'order1a' 'specialorder'.
The tif file of depths are converted to point clouds. Machine learning is used to identify the signal to noise ratio in the data. Filters are autommatically adapted to identify the 'most noisy' depths in the file. These candidates need to be evaluated against surrounding depths. The nearest neighbours are used to compute a 'regional depth surface' using a median value of the surrouding depths. The difference in candidate-regional depth (which we call DeltaZ) is then assessed against the TVU for
that depth and either flagged as an outlier or accepted as within specification. The remaining depths are called 'outliers'.
Inliers are points which do meet the required specification for allowable total vertical uncertainty.
Outliers are points which do NOT meet the required specification for allowable Total Vertical Uncertainty (TVU).
Outliers are saved to a point cloud file and a shape file. The shape file contains the processed depth, the Regional Depth, the AllowableTVU, the DeltaZ (difference between Regional Depth and
processed depth) and a field for Review/Approval by SIC.
GGOutlier does NOT modify the input file in any way. It is a read-only process.
GGOutlier will generate a QC report (this document) in order to enable rapid assessment of results.