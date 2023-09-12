import math
import rasterio
from rasterio.transform import from_origin
from rasterio.transform import Affine
import numpy as np

import geodetic

###############################################################################
def saveastif(outfilename, geo, pcd, resolution=1, fill=False):
	'''given a numpy array of point cloud, make a floating point geotif file using rasterio'''
	'''the numpy point clouds define the bounding box'''

	if len(pcd)==0:	
		return
		
	NODATA = -999
	xmin = pcd.min(axis=0)[0]
	ymin = pcd.min(axis=0)[1]
	zmin = pcd.min(axis=0)[2]
	
	xmax = pcd.max(axis=0)[0]
	ymax = pcd.max(axis=0)[1]
	zmax = pcd.max(axis=0)[2]

	xres 	= resolution
	yres 	= resolution
	width 	= math.ceil((xmax - xmin) / resolution)
	height 	= math.ceil((ymax - ymin) / resolution)

	transform = Affine.translation(xmin - xres / 2, ymin - yres / 2) * Affine.scale(xres, yres)
	
	print("Creating tif file... %s" % (outfilename))
	transform = from_origin(xmin, ymax, xres, yres)

	# save to file...
	src= rasterio.open(
			outfilename,
			mode="w",
			driver="GTiff",
			height=height,
			width=width,
			count=1,
			dtype='float32',
			crs=geo.projection.srs,
			transform=transform,
			nodata=NODATA,
	) 
	# populate the numpy array with the values....
	arr = np.full((height, width), fill_value=NODATA, dtype=float)
	
	from numpy import ma
	arr = ma.masked_values(arr, NODATA)

	for row in pcd:
		px = math.floor((row[0] - xmin) / xres)
		py = math.floor(height - (row[1] - ymin) / yres) - 1 #lord knows why -1
		# py, px = src.index(row[0], row[1])
		arr[py, px] = row[2]
		
	#we might want to fill in the gaps. useful sometimes...
	if fill:
		from rasterio.fill import fillnodata
		arr = fillnodata(arr, mask=None, max_search_distance=xres*2, smoothing_iterations=0)

	src.write(arr, 1)
	src.close()
	print("Creating tif file Complete.")

	return outfilename
###############################################################################
