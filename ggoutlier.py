#name:		  	GGOutlier
#created:		August 2023
#by:			paul.kennedy@guardiangeomatics.com
#description:   python module to read a point cloud or raster file of depths, identify outliers, write out points to a shap or laz file for further QC

#requirements
# python 3.10 (must be 3.10)
# pip install open3d
# pip install reportlab
# pip install pyproj
# pip install pyshp
# pip install scikit-learn
#done##########################################
# load a laz file
# compute the cloud2raster MEAN
# compute the cloud2raster MEDIAN
# compute the cloud2raster COUNT
# compute the cloud2raster STDDEV
# to clean outliers to IHO standard we need to do the following
# make a class to contain the IHO standard so we can easily compute the TVU for a given depth
# we need to compute a rough depth surface so we can compare the point cloud to the approximate surface.  
# This could be a surface of larger bin size so we get more soundings or it could be the median depth
# we then get the a subset of z percent of the most noisy depths
# we then compute the difference from the gridded surface
# if the depth difference is more than the permitted TVU we have an outlier.
# we need to flag the outlier
# allow user to specify threshold by iho order
# improve shop file with surface point depth att
# improve shop file with surface depth att
# improve shop file with SIC approval attribute
# load a tif file
# write results to laz
# allow user to specify threshold by percentage
# allow user to specify threshold by delta Z
#improve how we handle edge of data issues. most of outliers are edge of survey area.
#reduce marker size in pdf report
#improve metrics in pdf report
#improve pdf report to show the regional surface
#move to specific git repo

#todo##########################################
#load point cloud file from ascii fileutils

import os.path
from argparse import ArgumentParser
from datetime import datetime, timedelta
import math
import numpy as np
import open3d as o3d
import sys
import time
import glob
import rasterio
import multiprocessing as mp
import shapefile
import logging
import gc

# locals
import fileutils
import geodetic
import cloud2tif
import lashelper
import ggmbesstandard
import pdfdocument
###########################################################################
def main():

	iho = ggmbesstandard.sp44()
	msg = str(iho.getordernames())

	parser = ArgumentParser(description='Read a floating point TIF file of depths and find all outliers exceeding a user specified threshold.')
	parser.add_argument('-epsg', 	action='store', 		default="0",		dest='epsg', 			help='Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326')
	parser.add_argument('-i', 		action='store',			default="", 		dest='inputfile', 		help='Input filename/folder to process.')
	parser.add_argument('-odir', 	action='store', 		default="",			dest='odir', 			help='Specify a relative output folder e.g. -odir GIS')
	parser.add_argument('-n', 		action='store', 		default="20",		dest='numpoints', 		help='ADVANCED:Specify the number of nearest neighbours points to use.  More points means more data will be rejected. ADVANCED ONLY [Default:20]')
	parser.add_argument('-p', 		action='store', 		default="5",		dest='outlierpercentage',help='ADVANCED:Specify the approximate percentage of data to remove.  the engine will analyse the data and learn what filter settings are appropriate for your waterdepth and data quality. This is the most important (and only) parameter to consider spherical radius to find the nearest neightbours. Maximum is 5% [Default:5]')
	parser.add_argument('-z', 		action='store', 		default="10",		dest='zscale',			help='ADVANCED:Specify the ZScale to accentuate the depth difference ove the horizontal distance between points. Think of this as how you exxagerate the vertical scale in a swath editor to more easily spot the outliers. [Default:10]')
	parser.add_argument('-smooth', 	action='store', 		default="5",		dest='smooth',			help='ADVANCED:Specify the MEDIAN filter kernel width for computation of the regional surface so nearest neghbours can be calculated. [Default:5]')
	parser.add_argument('-standard',action='store', 		default="order1a",	dest='standard',		help='Specify the IHO SP44 survey order so we can set the filters to match the required specification. Select from :' + msg + ' [Default:order1a]' )
	parser.add_argument('-debug', 	action='store_true', 	default=False,		dest='debug',			help='DEbug to write LAZ files and other supproting file.s  takes some additional time!,e.g. -debug [deafult:false]')
	
	matches = []
	args = parser.parse_args()

	if os.path.isfile(args.inputfile):
		matches.append(args.inputfile)

	if len (args.inputfile) == 0:
		# no file is specified, so look for a .pos file in terh current folder.
		inputfolder = os.getcwd()
		matches = fileutils.findFiles2(False, inputfolder, "*.tif")

	if os.path.isdir(args.inputfile):
		matches = fileutils.findFiles2(False, args.inputfile, "*.tif")

	if len(matches) == 0:
		log("oops, no files found to process, quitting")
		exit(0)

	#make an output folder
	if len(args.odir) == 0:
		args.odir = str("GGOutlier_%s" % (time.strftime("%Y%m%d-%H%M%S")))
	odir = os.path.join(os.path.dirname(matches[0]), args.odir)
	makedirs(odir)

	logfilename = os.path.join(odir,"GGOutlier_log.txt").replace('\\','/')
	logging.basicConfig(filename = logfilename, level=logging.INFO)
	log("Configuration: %s" % (str(args)))
	log("Output Folder: %s" % (odir))
	log("GGOutlier Version: 1.01")
	log("GGOutlier started at: %s" % (datetime.now()))
	log("Username: %s" %(os.getlogin()))
	log("Computer: %s" %(os.environ['COMPUTERNAME']))
	log("Number of CPUs %d" %(mp.cpu_count()))	

	m = fileutils.MemoryStatusEx()
	log('You have %0.2f GiB of RAM installed' % (m.totalPhys / (1024.)**3))


	args.outlierpercentage = min(5.0, float(args.outlierpercentage))
	start_time = time.time() # time the process
	for file in matches:
		process(file, args)
		log("QC Duration:%.3fs" % (time.time() - start_time))
		pdfdocument.GGOutlierreport(logfilename, odir)
		log("Report Complete")

############################################################
def process(filename, args):
	'''we will try to find outliers using machine learning to determine the typical noise level in the file and adapt filter accordingly so the user-requested noise level is obatined.'''

	log("Processing file: %s" % (filename))

	log("QC to Survey Standard: %s" % (args.standard))
	iho = ggmbesstandard.sp44()
	standard = iho.loadstandard(args.standard)
	log("Survey_Standard: %s" %(standard.details()))

	#load the python proj projection object library if the user has requested it
	geo = geodetic.geodesy(args.epsg)
	log("EPSGCode for geodetic conversions: %s" % (args.epsg))

	pingcounter = 0
	beamcountarray = 0
	ZSCALE = float(args.zscale) # we might prefer 5 for this as this is how we like to 'look' for spikes in our data.  this value exaggerates the Z values thereby placing more emphasis on the Z than then X,Y

	#load the tif file...	
	with rasterio.open(filename) as src:
		band1 = src.read(1)
		z = band1.flatten()
		print('Band1 has shape', band1.shape)
		height = band1.shape[0]
		width = band1.shape[1]
		cols, rows = np.meshgrid(np.arange(width), np.arange(height))
		xs, ys = rasterio.transform.xy(src.transform, rows, cols)
		x = np.array(xs).flatten()
		y = np.array(ys).flatten()
		print('lons shape', x.shape)
		# src._crs.wkt
		xyz = np.stack((x,y,z), axis=1)
		NODATA = src.nodatavals[0]
		SRCRESOLUTION = src.res[0]
	
		# Delete the array
		del x
		del y
		del z
		del band1
		del xs
		del ys
		# Force garbage collection
		gc.collect()
		src.close()

	log("Loading Point Cloud...")
	pcd = o3d.geometry.PointCloud()
	#remove the NODATA values
	xyz = xyz[np.all(xyz != NODATA, axis=1)]
	#scale up the Z data so we accentuate the Z axis compared to the X,Y axis.  this is important as we use 
	xyz[:,2] *= ZSCALE
	pcd.points = o3d.utility.Vector3dVector(xyz)
	xyz[:,2] /= ZSCALE
	log("Depths loaded for quality control: %s" % (f'{len(pcd.points):,}'))

	if args.debug:
		#RAW report on RAW POINTS
		outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.splitext(os.path.basename(filename))[0] + "_RawPoints.txt")
		log ("Creating raw laz file of input raw points: %s " % (outfile))
		np.savetxt(outfile, xyz, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
		fname = lashelper.txt2las(outfile, epsg=args.epsg)
		fileutils.deletefile(outfile)

	del xyz
	# Force garbage collection
	gc.collect()

	# Populate the 'counter' field automatically so we can track the points accepted/rejected status
	beamcountarray = np.arange(0, len(pcd.points))  # This will populate 'counter' with values 1, 2, 3

	#PERCENTAGE PASSMARK
	log("Learning about your signal to noise levels...")
	start_time = time.time() # time the process
	low = 0
	high = 100
	TARGET = float(args.outlierpercentage)
	NUMPOINTS = max(int(args.numpoints),1)
	#MACHINE LEARNING TO DETERMINE CORRECT LEVEL OF FILTER...
	pcd, inlier_cloud, outlier_cloud, inlierindex = findoutlier(pcd, low, high, TARGET, NUMPOINTS)

	log ("Points accepted by machine learning: %d" % (len(inlier_cloud.points)))
	log ("Points tagged for further evaluation: %d" % (len(outlier_cloud.points)))
	# inliers = np.asarray(inlier_cloud.points)
	outliers = np.asarray(outlier_cloud.points)
	# inliers[:,2] /= ZSCALE
	outliers[:,2] /= ZSCALE
	#   log("QC Duration: %.3fs" % (time.time() - start_time))
	# log("QC Duration: %.3fseconds" % (time.time() - start_time)) # log the processing time. It is handy to keep an eye on processing performance.
	########

	#we need 1 list of ALL beams which are either accepted or rejected.
	beamqualityresult = np.isin(beamcountarray, inlierindex)

	#RAW save as a tif file so we can easily view in GIS...
	# outfilename = os.path.join(outfile + "_RegionalDepth.tif")
	#we need to make a regional grid which uses the nearest neighbours, so make this 3 times larger than the source grid.  this means 1 pixel each side of the current point
	# regionalfilename = lashelper.lasgridsubcircle( fname, outfilename, resolution= SRCRESOLUTION/2, epsg=args.epsg, subcircle=SRCRESOLUTION)
	regionalfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_RegionalDepth.tif")
	fname = cloud2tif.smoothtif(filename, regionalfilename, smooth=int(args.smooth))
	log ("Created REGIONAL TIF file for IHO validation: %s " % (fname))

	if args.debug:
		outfilename = os.path.join(outfile + "_RegionalDepth.laz")
		outfilename = lashelper.demzip2(regionalfilename, outfilename, nodata=-9999)
		log ("Created REGIONAL LAZ file of input raw points: %s " % (filename))

	# we can now double check the outliers to see how far they are away from the resulting inlier raster file of mean depths.  
	# if they are close then we can re-accept them
	rio = rasterio.open(regionalfilename)
	outlieridx = 0
	confirmedoutliers = []
	confirmedinliers = []
	for idx, validity in enumerate(beamqualityresult):
		if validity == True:
			continue
		else:
			#the beam quality is false which means we need to get the next value from the outliers list
			pt = outliers[outlieridx]
			#the outlier IDX is the sequential number used in conjunction with the beam quality results.  Its not great but thats how the open3d cleaning works
			outlieridx = outlieridx + 1
			depth = pt[2]
			# EDGE CLEAN UP: check if at edge of data by building a mask for regional and see if anything is an empty cell...			
			isedge=False
			radius = 2 #the median filter is set to 5 so we use 2 each side of central point
			mask = []
			for x in range(-radius, radius):
				for y in range(-radius, radius):
					mask.append([pt[0] + (x * SRCRESOLUTION), pt[1] + (y * SRCRESOLUTION)])
			for val in rio.sample(mask): 
				if val == rio.nodatavals[0]:
					griddepth = next(rio.sample([(pt[0], pt[1])]))[0]
					tvu = standard.gettvuat(griddepth)
					deltaz = abs(griddepth-depth)
					pt = np.append(pt, [deltaz, tvu, griddepth])
					confirmedinliers.append(pt)
					#re-accept the point as it is actually in specification
					inlierindex.append(idx)
					isedge=True
					break

			if isedge==True:
				continue
			#check the depth against the regional surface and the IHO standarON
			griddepth = next(rio.sample([(pt[0], pt[1])]))[0]
			tvu = standard.gettvuat(griddepth)
			deltaz = abs(griddepth-depth)
			if griddepth == rio.nodatavals[0]:
				#skip point if there is no raster surface value
				pt = np.append(pt, [deltaz, tvu, griddepth])
				confirmedinliers.append(pt)
				#re-accept the point as it is actually in specification
				inlierindex.append(idx)
				continue
			if deltaz >= tvu:
				pt = np.append(pt, [deltaz, tvu, griddepth])
				confirmedoutliers.append(pt)
			else:
				pt = np.append(pt, [deltaz, tvu, griddepth])
				confirmedinliers.append(pt)
				#re-accept the point as it is actually in specification
				inlierindex.append(idx)

	rio.close()
	
	#write the outliers to a point SHAPE file
	shpfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_OutlierPoints" + ".shp")
	w = shapefile.Writer(shpfilename)
	w.field('DeltaZ', 'N',8,3) # 8 byte floats, 3 decimal places
	w.field('AllowedTVU', 'N',8,3)
	w.field('Depth', 'N',8,3)
	w.field('GridDepth', 'N',8,3)
	w.field('SICApproved','C','254')
	log ("Writing outliers to: %s" % (shpfilename))

	for idx, pt in enumerate(confirmedoutliers):
		w.pointz(pt[0], pt[1], pt[2])
		w.record(pt[3], pt[4], pt[2], pt[5])
	w.close()

	# q: what is a median filter and how does it work?
	# a: https://www.youtube.com/watch?v=VvQ2mo3EXHk
	#we need 1 list of ALL beams which are either accepted or rejected.  This is now a revised list folowing Standard validation
	inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
	outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)
	percentage 		= (100 * (len(outlier_cloud.points) / len(pcd.points)))
	log ("Points accepted: %s" % (f'{len(inlier_cloud.points):,}'))
	log ("Points outside specification: %s" % (f'{len(outlier_cloud.points):,}'))
	log ("Percentage outside specification: %.4f" % (percentage))
	inliers = np.asarray(inlier_cloud.points)
	outliers = np.asarray(outlier_cloud.points)
	inliers[:,2] /= ZSCALE
	outliers[:,2] /= ZSCALE

	if args.debug:
		#INLIERS reporting...
		outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_InlierPoints" + ".txt")
		np.savetxt(outfile, inliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
		outfilename = os.path.join(outfile + "_Depth.tif")
		# cloud2tif.saveastif(outfilename, geo, inliers, fill=False)
		inlierraster = cloud2tif.point2raster(outfilename, geo, inliers, resolution = 1, bintype='mean', fill=False)
		#write the outliers to a point cloud laz file
		fname = lashelper.txt2las(outfile, epsg=args.epsg)
		lashelper.lasgrid4( fname, outfilename, resolution=1, epsg=args.epsg)
		# fileutils.deletefile(outfile)
		log ("Created LAZ file of inliers: %s " % (fname))

	#OUTLIERS reporting...
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_OutlierPoints" + ".txt")
	np.savetxt(outfile, outliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	log ("Created TXT file of outliers: %s " % (outfile))
	#write the outliers to a point cloud laz file
	fname = lashelper.txt2las(outfile, epsg=args.epsg)
	log ("Created LAZ file of outliers: %s " % (fname))

	# msg = "GGOutlier complete. outlier_cloud.points
	log("QC complete at: %s" % (datetime.now()))
	return shpfilename

##################################################################################
def findoutlier(pcd, low, high, TARGET=1.0, NUMPOINTS=3):
	'''clean outliers using binary chop to control how many points we reject'''
	'''use spherical radius to identify outliers and clusters'''
	'''binary chop will aim for target percentage of data deleted rather than a fixed filter level'''
	'''this way the filter adapts to the data quality'''
	'''TARGET is the percentage of the input points we are looking to reject'''
	'''NUMPOINTS is the number of nearest neighbours within the spherical radius which is the threshold we use to consider a point an outlier.'''
	'''If a point has no friends, then he is an outlier'''
	'''if a point has moew the NUMPOINTS in the spherical radius then he is an inlier, ie good'''

	# Force garbage collection
	gc.collect()

	#outlier removal by radius
	# http://www.open3d.org/docs/latest/tutorial/geometry/pointcloud_outlier_removal.html?highlight=outlier
	# http://www.open3d.org/docs/latest/tutorial/Advanced/pointcloud_outlier_removal.html
	
	#cl: The pointcloud as it was fed in to the model (for some reason, it seems a bit pointless to return this).
	#ind: The index of the points which are NOT outliers
	currentfilter = (high+low)/2

	cl, inlierindex = pcd.remove_statistical_outlier(nb_neighbors=NUMPOINTS,	std_ratio=currentfilter)
	# cl, inlierindex = pcd.remove_radius_outlier(nb_points = NUMPOINTS, radius = currentfilter)

	inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
	outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)
	percentage 		= (100 * (len(outlier_cloud.points) / len(pcd.points)))
	log ("Current filter StdDev %.2f" % (currentfilter))
	log ("Percentage Rejection %.2f" % (percentage))

	decimals = len(str(TARGET).split(".")[1])
	percentage = round(percentage, decimals)
	if percentage < TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		log ("Filter level increasing to reject a few more points...")
		del inlier_cloud
		del outlier_cloud
		del cl 				
		del inlierindex
		pcd, inlier_cloud, outlier_cloud, inlierindex = findoutlier(pcd, low, currentfilter, TARGET, NUMPOINTS)
	elif percentage > TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		log ("Filter level decreasing to reject a few less points...")
		del inlier_cloud
		del outlier_cloud
		del cl 				
		del inlierindex
		pcd, inlier_cloud, outlier_cloud, inlierindex = findoutlier(pcd, currentfilter, high, TARGET, NUMPOINTS)

	return (pcd, inlier_cloud, outlier_cloud, inlierindex)

###############################################################################
def update_progress(job_title, progress):
	'''progress value should be a value between 0 and 1'''
	length = 20 # modify this to change the length
	block = int(round(length*progress))
	msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
	if progress >= 1: msg += " DONE\r\n"
	sys.stdout.write(msg)
	sys.stdout.flush()

###############################################################################
def	makedirs(odir):
	if not os.path.isdir(odir):
		os.makedirs(odir, exist_ok=True)

###############################################################################
def	log(msg, error = False, printmsg=True):
		if printmsg:
			print (msg)
		if error == False:
			logging.info(msg)
		else:
			logging.error(msg)

###############################################################################
if __name__ == "__main__":
		main()

	########################v#######################################################
	# log("Statistical outlier removal")
	# voxel_down_pcd = pcd.voxel_down_sample(voxel_size=0.001)
	# voxel_down_pcd = pcd
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=3.0) # 1.51
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=10, std_ratio=3.0) # 1.89
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=10, std_ratio=1.0) 
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=2.0) # 3.54%
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=1.0) # 9.56

	# obb = pcd.get_oriented_bounding_box()
	# obb.color = (0,0,0)
	# display_inlier_outlier(voxel_down_pcd, ind)

	# o3d.visualization.draw_geometries([pcd, obb])

	# pc = open3d.io.read_point_cloud(outfile, format='xyz')
	# print (pcd)
	# eps = 0.1  # DBSCAN epsilon parameter
	# min_samples = 1  # DBSCAN minimum number of points
	# despike_point_cloud(xyz, eps, min_samples)

	# eps = 0.1  # DBSCAN epsilon parameter
	# min_samples = 3  # DBSCAN minimum number of points
	# despike_point_cloud(xyz, eps, min_samples)

	# eps = 0.1  # DBSCAN epsilon parameter
	# min_samples = 10  # DBSCAN minimum number of points
	# despike_point_cloud(xyz, eps, min_samples)


	# eps = 0.01  # DBSCAN epsilon parameter
	# min_samples = 3  # DBSCAN minimum number of points
	# despike_point_cloud(xyz, eps, min_samples)

	# eps = 0.05  # DBSCAN epsilon parameter
	# min_samples = 3  # DBSCAN minimum number of points
	# despike_point_cloud(xyz, eps, min_samples)

	# print ("DBSCAN...")
	# xrange = max(xyz[:,0]) - min(xyz[:,0])
	# yrange = max(xyz[:,1]) - min(xyz[:,1])
	# maxrange = max(xrange, yrange)
	# mediandepth = statistics.median(xyz[:, 2])
	# print ("WaterDepth %.2f" % (mediandepth))
	# eps = mediandepth * 0.05 # 1% waterdepth  bigger number rejects fewer points
	# # eps = 0.1  # DBSCAN epsilon parameter
	# min_samples = 5  # DBSCAN minimum number of points
	# rejected = despike_point_cloud(xyz, eps, min_samples)
	# print ("DBSCAN Complete")
	# print ("Percentage rejected %.2f" % (len(rejected)/ len(xyz) * 100))	
	# fig = plt.figure(figsize=(10, 6))
	# ax = fig.add_subplot(111, projection='3d')
	# # create light source object.
	# # ls = LightSource(azdeg=0, altdeg=65)
	
	# # shade data, creating an rgb array.
	# # rgb = ls.shade(z, plt.cm.RdYlBu)
	
	# zrange = max(xyz[:,2]) - min(xyz[:,2])
	# xyzdisplay = xyz[::2]
	# ax.scatter(xyzdisplay[:, 0], xyzdisplay[:, 1], xyzdisplay[:, 2], color = 'lightgrey', s=5)
	# ax.scatter(rejected[:, 0], rejected[:, 1], rejected[:, 2], color = 'red', s=50)
	# ax.set_xlim3d(min(xyz[:,0]), min(xyz[:,0]) + maxrange)
	# zscale = 5
	# ax.set_zlim3d(min(xyz[:,1]), (min(xyz[:,1]) + maxrange) * 5)
	# ax.set_zlim3d(min(xyz[:,2]), (min(xyz[:,2]) + maxrange) * 5)

	# plt.show()

	# cloud2tif.point2raster(filename + "_median.tif", geo, xyz, resolution=1, bintype='median', fill=False)
	# cloud2tif.point2raster(filename + "_count.tif", geo, xyz, resolution=1, bintype='count', fill=False)
	# cloud2tif.point2raster(filename + "_STD.tif", geo, xyz, resolution=1, bintype='stddev', fill=False)

	# Code to generate example dataset
	# xyz = np.random.uniform(0.0, 1000.0, size=(1000,3))
	# xyz = np.asarray([[0,5,10],[0,6,10],[0,7,10],[1,5,11],[1,6,11],[1,7,11],[2,5,12],[2,6,12],[2,7,12], [0,5,999], [0,5,50], [0,5,51], [0,5,52]])

	#SINGLE PASS...
	# start_time = time.time() # time the process
	# currentfilter = 20 #bigger removes fewer points
	# NUMPOINTS = int(args.numpoints)
	# cl, inlierindex = pcd.remove_statistical_outlier(nb_neighbors=NUMPOINTS,	std_ratio=currentfilter)
	# inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
	# outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)
	# cl, inlierindex = pcd.remove_radius_outlier(nb_points = NUMPOINTS, radius = currentfilter)
