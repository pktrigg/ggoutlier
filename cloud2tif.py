import math
import rasterio
from rasterio.transform import from_origin
from rasterio.transform import Affine
import numpy as np

import geodetic
import logging

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
	
	log("Creating tif file... %s" % (outfilename))
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
	log("Creating tif file Complete.")

	return outfilename

###############################################################################
def pcd2meantif(outfilename, geo, pcd, resolution=1, fill=False):
    # Current (inefficient) code to quantize into XY 'bins' and take mean Z values in each bin

	if len(pcd)==0:	
		return
		
	pcd[:, 0:2] = np.round(pcd[:, 0:2]/float(resolution))*float(resolution) # Round XY values to nearest resolution value

	NODATA = -999
	xmin = pcd.min(axis=0)[0]
	ymin = pcd.min(axis=0)[1]
	
	xmax = pcd.max(axis=0)[0]
	ymax = pcd.max(axis=0)[1]

	xres 	= resolution
	yres 	= resolution
	width 	= math.ceil((xmax - xmin) / resolution)
	height 	= math.ceil((ymax - ymin) / resolution)
	mean_height = np.zeros((height, width))

	# Loop over each x-y bin and calculate mean z value 
	x_val = xmin
	for x in range(width):
		y_val = ymax
		for y in range(height):
			height_vals = pcd[(pcd[:,0] == float(x_val)) & (pcd[:,1] == float(y_val)), 2]
			if height_vals.size != 0:
				mean_height[y,x] = np.mean(height_vals)
			y_val -= resolution
		x_val += resolution

	# return mean_height
	arr = mean_height
	arr[mean_height == 0] = NODATA
	
	log("Creating tif file... %s" % (outfilename))
	transform = from_origin(xmin-(xres/2), ymax + (yres/2), xres, yres)

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
	#we might want to fill in the gaps. useful sometimes...
	if fill:
		from rasterio.fill import fillnodata
		arr = fillnodata(arr, mask=None, max_search_distance=xres*2, smoothing_iterations=0)

	src.write(arr, 1)
	src.close()
	log("Creating tif file Complete.")

	return outfilename

###############################################################################
def point2raster(outfilename, geo, pcd, resolution=1, bintype="mean", fill=False):
	'''given a numpy array of point cloud, make a floating point geotif file using rasterio'''
	'''the numpy point clouds define the bounding box'''
	# https://stackoverflow.com/questions/54842690/how-to-efficiently-convert-large-numpy-array-of-point-cloud-data-to-downsampled

	NODATA = -999
	
	if len(pcd)==0:
		return

	#take the point cloud array and transpose the xyz,xyz,xyz into xxx,yyy so we can bin them efficienctly without looping thru the data
	xy = pcd.T[:2]
	#bin the xy data into buckets.  at present this is only integer based so 1m resolution is minimum
	xy = ((xy + resolution / 2) // resolution).astype(int)
	# xy = ((xy - resolution / 2) // resolution).astype(int)
	#compute the range of the data 
	mn, mx = xy.min(axis=1), xy.max(axis=1)
	#compute the size of the data
	sz = mx + 1 - mn

	if bintype == 'mean':
		#Converts a tuple of index arrays into an array of flat indices, applying boundary modes to the multi-index.
		#RETURNS An array of indices into the flattened version of an array of dimensions dims.
		flatidx = np.ravel_multi_index(xy-mn[:, None], dims=sz)
		#compute the mean of each bin as efficiently as possible
		histo = np.bincount(flatidx, pcd[:, 2], sz.prod()) / np.maximum(1, np.bincount(flatidx, None, sz.prod()))
		arr = histo.reshape(sz).T
		arr = np.flip(arr, axis = 0)

	if bintype == 'count':
		#Converts a tuple of index arrays into an array of flat indices, applying boundary modes to the multi-index.
		#RETURNS An array of indices into the flattened version of an array of dimensions dims.
		flatidx = np.ravel_multi_index(xy-mn[:, None], dims=sz)
		#we can compute the count rapidly as well...
		histo = np.maximum(0, np.bincount(flatidx, None, sz.prod()))
		arr = histo.reshape(sz).T
		arr = np.flip(arr, axis = 0)

	if bintype == 'median':
		#calculate the medians...
		#https://stackoverflow.com/questions/10305964/quantile-median-2d-binning-in-python
		# Median is a bit harder
		flatidx = np.ravel_multi_index(xy-mn[:, None], dims=sz)
		order = flatidx.argsort()
		bin = flatidx[order]
		w = pcd[:, 2][order]
		edges = (bin[1:] != bin[:-1]).nonzero()[0] + 1
		# Median 
		median = [np.median(i) for i in np.split(w, edges)]
		#construct BINSxBINS matrix with median values
		binvals=np.unique(bin)
		medvals=np.zeros([sz.prod()])
		medvals[binvals]=median
		medvals=medvals.reshape(sz)
		arr = np.asarray(medvals).reshape(sz).T
		arr = np.flip(arr, axis = 0)

	if bintype == 'stddev':
		#https://stackoverflow.com/questions/10305964/quantile-median-2d-binning-in-python
		# Median is a bit harder
		flatidx = np.ravel_multi_index(xy-mn[:, None], dims=sz)
		order = flatidx.argsort()
		bin = flatidx[order]
		w = pcd[:, 2][order]
		edges = (bin[1:] != bin[:-1]).nonzero()[0] + 1
		# Standard Deviation
		stddev = [np.std(i) for i in np.split(w, edges)]
		#construct BINSxBINS matrix with median values
		binvals=np.unique(bin)
		sdvals=np.zeros([sz.prod()])
		sdvals[binvals]=stddev
		sdvals=sdvals.reshape(sz)
		arr = np.asarray(sdvals).reshape(sz).T
		arr = np.flip(arr, axis = 0)

	# clear out the empty nodes and set to NODATA value
	arr[arr == 0] = NODATA

	xmin = mn[0]
	ymin = mn[1]
	xmax = mx[0]
	ymax = mx[1]
	xres 	= resolution
	yres 	= resolution
	
	width 	= math.ceil((xmax - xmin) / resolution)
	height 	= math.ceil((ymax - ymin) / resolution)

	log("Creating tif file... %s" % (outfilename))
	transform = from_origin(xmin-(xres/2), ymax + (yres/2), xres, yres)

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
	#we might want to fill in the gaps. useful sometimes...
	if fill:
		from rasterio.fill import fillnodata
		arr = fillnodata(arr, mask=None, max_search_distance=xres*2, smoothing_iterations=0)

	src.write(arr, 1)
	src.close()
	log("Creating tif file Complete.")

	return outfilename
###############################################################################

###############################################################################
def	log(msg, error = False, printmsg=True):
		if printmsg:
			print (msg)
		if error == False:
			logging.info(msg)
		else:
			logging.error(msg)
