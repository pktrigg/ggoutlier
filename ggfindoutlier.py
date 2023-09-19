#name:		  	ggfindoutlier
#created:		August 2023
#by:			paul.kennedy@guardiangeomatics.com
#description:   python module to read a point cloud or raster file of depths, identify outliers, write out points to a shap or laz file for further QC

#done##########################################
# load a laz file
# compute the cloud2raster MEAN
# compute the cloud2raster MEDIAN
# compute the cloud2raster COUNT
# compute the cloud2raster STDDEV

#todo##########################################
#load point cloud file from ascii fileutils
# load a tif file
# write results to laz
# allow user to specify threshol by iho order
# allow user to specify threshold by percentage
# allow user to specify threshold by delta Z

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

import kmall
import fileutils
import geodetic
import multiprocesshelper 
import cloud2tif
import lashelper
import ggmbesstandard

###########################################################################
def main():

	iho = ggmbesstandard.sp44()
	msg = str(iho.getorders())

	parser = ArgumentParser(description='Read a KMALL file.')
	parser.add_argument('-epsg', 	action='store', 		default="0",	dest='epsg', 			help='Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326')
	parser.add_argument('-i', 		action='store',			default="", 	dest='inputfile', 		help='Input filename/folder to process.')
	parser.add_argument('-cpu', 	action='store', 		default='0', 	dest='cpu', 			help='number of cpu processes to use in parallel. [Default: 0, all cpu]')
	parser.add_argument('-odir', 	action='store', 		default="",	dest='odir', 			help='Specify a relative output folder e.g. -odir GIS')
	parser.add_argument('-n', 		action='store', 		default="1",	dest='numpoints', 		help='Specify the number of nearest neighbours points to use.  More points means more data will be rejected. ADVANCED ONLY [Default:1]')
	parser.add_argument('-p', 		action='store', 		default="0.1",	dest='outlierpercentage',help='Specify the approximate percentage of data to remove.  the engine will analyse the data and learn what filter settings are appropriate for your waterdepth and data quality. This is the most important (and only) parameter to consider spherical radius to find the nearest neightbours. [Default:0.1]')
	parser.add_argument('-z', 		action='store', 		default="1.0",	dest='zscale',			help='Specify the ZScale to accentuate the depth difference ove the horizontal distance between points. Think of this as how you exxagerate the vertical scale in a swath editor to more easily spot the outliers. [Default:1.0]')
	parser.add_argument('-sp44', 	action='store', 		default="order1a",	dest='sp44',		help='Specify the IHO SP44 survey order so we can set the filters to match the required specification. Select from :' + msg + ' [Default:order1a]' )
	parser.add_argument('-debug', 	action='store', 		default="-1",	dest='debug', 			help='Specify the number of pings to process.  good only for debugging. [Default:-1]')
	
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
		args.odir = str("NearestNeighbours_%d_OutliersPercent_%.2f_zscale_%.2f" % (int(args.numpoints), float(args.outlierpercentage), float(args.zscale)))
	odir = os.path.join(os.path.dirname(matches[0]), args.odir)
	makedirs(odir)

	logging.basicConfig(filename = os.path.join(odir,"kmallclean_log.txt"), level=logging.INFO)
	log("configuration: %s" % (str(args)))
	log("Output Folder: %s" % (odir))

	results = []
	if args.cpu == '1':
		for file in matches:
			findoutlier(file, args)
	else:
		multiprocesshelper.log("Files to Import: %d" %(len(matches)))		
		cpu = multiprocesshelper.getcpucount(args.cpu)
		log("Processing with %d CPU's" % (cpu))

		pool = mp.Pool(cpu)
		multiprocesshelper.g_procprogress.setmaximum(len(matches))
		poolresults = [pool.apply_async(findoutlier, (file, args), callback=multiprocesshelper.mpresult) for file in matches]
		pool.close()
		pool.join()
		# for idx, result in enumerate (poolresults):
		# 	results.append([file, result._value])
		# 	print (result._value)

############################################################
def findoutlier(filename, args):
	'''we will try to find outliers using machine learning to determine the typical noise level in the file and adapt filter accordingly so the user-requested noise level is obatined.'''

	log("Processing file: %s" % (filename))

	maxpings = int(args.debug)
	if maxpings == -1:
		maxpings = 999999999

	pingcounter = 0
	beamcountarray = 0
	ZSCALE = float(args.zscale) # we might prefer 5 for this as this is how we like to 'look' for spikes in our data.  this value exaggerates the Z values thereby placing more emphasis on the Z than then X,Y
	
	# if os.path.splitext(filename)[1].lower() == ".tif":
		# rio = rasterio.open(filename)
		# arr = rio.read(1)  # read all raster values

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
	# if args.epsg == '0':
	# 	approxlongitude, approxlatitude = r.getapproximatepositon()
	# 	args.epsg = geodetic.epsgfromlonglat (approxlongitude, approxlatitude)

	#load the python proj projection object library if the user has requested it
	geo = geodetic.geodesy(args.epsg)
	log("EPSGCode for geodetic conversions: %s" % (args.epsg))
	
	log("Loading Point Cloud...")
	#get the record count so we can show a progress bar

	pcd = o3d.geometry.PointCloud()

	#remove the NODATA values
	xyz = xyz[np.all(xyz != NODATA, axis=1)]
	xyz[:,2] *= ZSCALE
	pcd.points = o3d.utility.Vector3dVector(xyz)
	xyz[:,2] /= ZSCALE

	# Code to generate example dataset
	# xyz = np.random.uniform(0.0, 1000.0, size=(1000,3))
	# xyz = np.asarray([[0,5,10],[0,6,10],[0,7,10],[1,5,11],[1,6,11],[1,7,11],[2,5,12],[2,6,12],[2,7,12], [0,5,999], [0,5,50], [0,5,51], [0,5,52]])

	# cloud2tif.point2raster(filename + "_mean.tif", geo, xyz, resolution=1, bintype='mean', fill=False)
	# cloud2tif.point2raster(filename + "_median.tif", geo, xyz, resolution=1, bintype='median', fill=False)
	# cloud2tif.point2raster(filename + "_count.tif", geo, xyz, resolution=1, bintype='count', fill=False)
	cloud2tif.point2raster(filename + "_STD.tif", geo, xyz, resolution=1, bintype='stddev', fill=False)

	log("Depths loaded for cleaning: %s" % (f'{len(pcd.points):,}'))

	# Populate the 'counter' field automatically
	beamcountarray = np.arange(0, len(pcd.points))  # This will populate 'counter' with values 1, 2, 3

	#PERCENTAGE PASS
	log("Understanding your data noise levels...")
	start_time = time.time() # time the process
	low = 0
	high = 100
	TARGET = float(args.outlierpercentage)
	NUMPOINTS = max(int(args.numpoints),1)
	pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier2(pcd, low, high, TARGET, NUMPOINTS)

	#SINGLE PASS
	# start_time = time.time() # time the process
	# currentfilter = 20 #bigger removes fewer points
	# NUMPOINTS = int(args.numpoints)
	# cl, inlierindex = pcd.remove_statistical_outlier(nb_neighbors=NUMPOINTS,	std_ratio=currentfilter)
	# inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
	# outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)

	# cl, inlierindex = pcd.remove_radius_outlier(nb_points = NUMPOINTS, radius = currentfilter)

	log ("Points accepted: %.2f" % (len(inlier_cloud.points)))
	log ("Points rejected: %.2f" % (len(outlier_cloud.points)))
	inliers = np.asarray(inlier_cloud.points)
	outliers = np.asarray(outlier_cloud.points)
	inliers[:,2] /= ZSCALE
	outliers[:,2] /= ZSCALE
	log("Clean Duration: %.3f seconds" % (time.time() - start_time)) # log the processing time. It is handy to keep an eye on processing performance.
	########

	#we need 1 list of ALL beams which are either accepted or rejected.
	beamqualityresult = np.isin(beamcountarray, inlierindex)

	#report on RAW POINTS
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_R.txt")
	xyz[:,2] /= ZSCALE
	np.savetxt(outfile, xyz, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	fname = lashelper.txt2las(outfile, epsg=args.epsg)
	#save as a tif file...
	outfilename = os.path.join(outfile + "_Depth.tif")
	lashelper.lasgrid4( fname, outfilename, resolution=1, epsg=args.epsg)
	
	fileutils.deletefile(outfile)
	log ("Created LAZ file of input raw points: %s " % (fname))
	outfilename = os.path.join(outfile + "_Raw_Depth.tif")
	# raw = np.asarray(pcd.points)
	# raw[:,2] /= ZSCALE
	# cloud2tif.saveastif(outfilename, geo, raw, fill=False)
	# outfilename = os.path.join(outfile + "_R_NEW.tif")
	# cloud2tif.pcd2meantif2(outfilename, geo, raw, fill=False)

	#report on INLIERS
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_Inlier" + ".txt")
	np.savetxt(outfile, inliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	outfilename = os.path.join(outfile + "_Depth.tif")
	# cloud2tif.saveastif(outfilename, geo, inliers, fill=False)
	# inlierraster = cloud2tif.pcd2meantif(outfilename, geo, inliers, fill=False)
	#write the outliers to a point cloud laz file
	fname = lashelper.txt2las(outfile, epsg=args.epsg)
	lashelper.lasgrid4( fname, outfilename, resolution=1, epsg=args.epsg)
	# fileutils.deletefile(outfile)
	log ("Created LAZ file of inliers: %s " % (fname))

	#report on OUTLIERS
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_Outlier" + ".txt")
	np.savetxt(outfile, outliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	#write the outliers to a point cloud laz file
	fname = lashelper.txt2las(outfile, epsg=args.epsg)
	log ("Created LAZ file of outliers: %s " % (fname))

	# we can now double check the outliers to see how far they are away from the resulting inlier raster file of mean depths.  
	# if they are close then we can re-accept them
	# outlieridx = 0
	# rio = rasterio.open(inlierraster)
	# for idx, validity in enumerate(beamqualityresult):
	# 	if validity == True:
	# 		continue
	# 	else:
	# 		pt = outliers[outlieridx]
	# 		outlieridx = outlieridx + 1
	# 		# griddepth = rio.sample([(pt[0], pt[1])])
	# 		griddepth = next(rio.sample([(pt[0], pt[1])]))[0]
	# 		log ((griddepth-pt[2]))
	# 		if abs(pt[2] - griddepth) > abs((griddepth*0.05)):
	# 			log ("confirmed")
	# 		else:
	# 			log ("re-accept")


	#write the outliers to a point SHAPE file
	# outfilename = os.path.join(outfile + ".shp")
	# w = shapefile.Writer(outfilename)
	# # for point in outlier_cloud.poin:
	# w.multipoint(outliers.tolist())
	# w.field('name', 'C')
	# w.record('outlier')
	# w.close()

	#now lets write out a NEW KMALL file with the beams modified...
	#create an output file....
	outfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename))
	outfilename = fileutils.addFileNameAppendage(outfilename, "_CLEANED")
	outfileptr = open(outfilename, 'wb')

	log("Writing NEW KMALL file %s" % (outfilename))
	pingcounter = 0
	beamcounter = 0

	log("Cleaning complete at: %s" % (datetime.now()))
	return outfilename

###############################################################################
# def setbeamquality(datagram, beamcounter, inlierindex):
# 	'''apply the cleaning results to the ping of data'''
# 	test_set=set(inlierindex)
# 	for idx, beam in enumerate(datagram.beams):
# 		if not beamcounter+idx in test_set: 
# 		# if not beamcounter+idx in inlierindex: 
# 			beam.detectionType = 2
# 		else:
# 			pk=1
# 	return

###############################################################################
def validateoutliers(inlierraster, outlier_cloud):

	pcd = np.asarray(outlier_cloud.points)
	for row in pcd:
		# py, px = inlierraster.index(row[0], row[1])
		v = inlierraster.sample(row[0], row[1])

	# py, px = inlierraster.index(row[0], row[1])

	return outlier_cloud

##################################################################################
def cleanoutlier2(pcd, low, high, TARGET=1.0, NUMPOINTS=3):
	'''clean outliers using binary chop to control how many points we reject'''
	'''use spherical radius to identify outliers and clusters'''
	'''binary chop will aim for target percentage of data deleted rather than a fixed filter level'''
	'''this way the filter adapts to the data quality'''
	'''TARGET is the percentage of the input points we are looking to reject'''
	'''NUMPOINTS is the number of nearest neighbours within the spherical radius which is the threshold we use to consider a point an outlier.'''
	'''If a point has no friends, then he is an outlier'''
	'''if a point has moew the NUMPOINTS in the spherical radius then he is an inlier, ie good'''

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
	log ("Current filter Nearest Neighbours %d" % (NUMPOINTS))
	log ("Current filter StdDEv %.2f" % (currentfilter))
	log ("Percentage rejection %.2f" % (percentage))

	decimals = len(str(TARGET).split(".")[1])
	percentage = round(percentage, decimals)
	if percentage < TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		log ("Filter level increasing to reject a few more points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier2(pcd, low, currentfilter, TARGET, NUMPOINTS)
	elif percentage > TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		log ("Filter level decreasing to reject a few less points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier2(pcd, currentfilter, high, TARGET, NUMPOINTS)

	return (pcd, inlier_cloud, outlier_cloud, inlierindex)

##################################################################################
def cleanoutlier(pcd, low, high, TARGET=1.0, NUMPOINTS=3):
	'''clean outliers using binary chop to control how many points we reject'''
	'''use spherical radius to identify outliers and clusters'''
	'''binary chop will aim for target percentage of data deleted rather than a fixed filter level'''
	'''this way the filter adapts to the data quality'''
	'''TARGET is the percentage of the input points we are looking to reject'''
	'''NUMPOINTS is the number of nearest neighbours within the spherical radius which is the threshold we use to consider a point an outlier.'''
	'''If a point has no friends, then he is an outlier'''
	'''if a point has moew the NUMPOINTS in the spherical radius then he is an inlier, ie good'''

	#outlier removal by radius
	# http://www.open3d.org/docs/latest/tutorial/geometry/pointcloud_outlier_removal.html?highlight=outlier
	# http://www.open3d.org/docs/latest/tutorial/Advanced/pointcloud_outlier_removal.html
	
	#cl: The pointcloud as it was fed in to the model (for some reason, it seems a bit pointless to return this).
	#ind: The index of the points which are NOT outliers
	currentfilter = (high+low)/2
	cl, inlierindex = pcd.remove_radius_outlier(nb_points = NUMPOINTS, radius = currentfilter)

	inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
	outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)
	percentage 		= (100 * (len(outlier_cloud.points) / len(pcd.points)))
	log ("Current filter radius %.2f" % (currentfilter))
	log ("Percentage rejection %.2f" % (percentage))

	decimals = len(str(TARGET).split(".")[1])
	percentage = round(percentage, decimals)
	if percentage < TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		log ("Filter level increasing to reject a few more points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier(pcd, low, currentfilter, TARGET, NUMPOINTS)
	elif percentage > TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		log ("Filter level decreasing to reject a few less points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier(pcd, currentfilter, high, TARGET, NUMPOINTS)

	return (pcd, inlier_cloud, outlier_cloud, inlierindex)

###############################################################################
def clipper(datagram, clip):
	'''using the datagram, reject if the take off angle is outside the clip limit'''

	for beam in datagram.beams:
		if abs(beam.beamAngleReRx_deg) > clip:
			beam.detectionType = 2 # reject the beam
			beam.detectionType = 0 # no valid detect
			# beam.rejectionInfo1 = ??

###############################################################################
def display_inlier_outlier(cloud, ind):
	inlier_cloud = cloud.select_by_index(ind)
	outlier_cloud = cloud.select_by_index(ind, invert=True)
	log (inlier_cloud)
	log (outlier_cloud)
	log ("Percentage rejection %.2f" % (100 * (len(outlier_cloud.points) / len(inlier_cloud.points))))
	log("Showing outliers (red) and inliers (gray): ")
	outlier_cloud.paint_uniform_color([1, 0, 0])
	inlier_cloud.paint_uniform_color([0.8, 0.8, 0.8])

	# hull = inlier_cloud.compute_convex_hull()
	# hull_ls = o3d.geometry.LineSet.create_from_triangle_mesh(hull)
	# hull_ls.paint_uniform_color((1, 0, 0))
	# hull_ls = o3d.geometry.LineSet.create_from_triangle_mesh(hull.to to_legacy())
	# hull.paint_uniform_color((1, 0, 0))

	o3d.visualization.draw_geometries([inlier_cloud, outlier_cloud])
										# zoom=0.3412,
										# front=[0.4257, -0.2125, -0.8795],
										# lookat=[2.6172, 2.0475, 1.532],
										# up=[-0.0694, -0.9768, 0.2024])

###############################################################################
###############################################################################
def computebathypointcloud(datagram, geo):
	'''using the MRZ datagram, efficiently compute a numpy array of the point clouds  '''

	for beam in datagram.beams:
		beam.east, beam.north = geo.convertToGrid((beam.deltaLongitude_deg + datagram.longitude), (beam.deltaLatitude_deg + datagram.latitude))
		beam.depth = (beam.z_reRefPoint_m + datagram.txTransducerDepth_m) * -1.0 #invert depths so we have negative depths.
		# beam.depth = beam.z_reRefPoint_m - datagram.z_waterLevelReRefPoint_m
		# beam.id			= datagram.pingCnt

	npeast = np.fromiter((beam.east for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	npnorth = np.fromiter((beam.north for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	npdepth = np.fromiter((beam.depth for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	npq = np.fromiter((beam.rejectionInfo1 for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	# npid = np.fromiter((beam.id for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)

	# we can now comput absolute positions from the relative positions
	# npLatitude_deg = npdeltaLatitude_deg + datagram.latitude_deg	
	# npLongitude_deg = npdeltaLongitude_deg + datagram.longitude_deg
	return (npeast, npnorth, npdepth, npq)

###############################################################################
# def despike_point_cloud(points, eps, min_samples):
# 	"""Despike a point cloud using DBSCAN."""
# 	clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(points)
# 	labels = clustering.labels_
# 	filtered_points = points[labels != -1]
# 	rejected_points = points[labels == -1]
    
# 	log("EPS: %f MinSample: %f Rejected: %d Survivors: %d InputCount %d" % (eps,  min_samples, len(rejected_points), len(filtered_points), len(points)))
# 	return rejected_points


# ###############################################################################
# def findFiles2(recursive, filespec, filter):
# 	'''tool to find files based on user request.  This can be a single file, a folder start point for recursive search or a wild card'''
# 	matches = []
# 	if recursive:
# 		matches = glob(os.path.join(filespec, "**", filter), recursive = False)
# 	else:
# 		matches = glob(os.path.join(filespec, filter), recursive = False)
	
# 	mclean = []
# 	for m in matches:
# 		mclean.append(m.replace('\\','/'))
		
# 	# if len(mclean) == 0:
# 	# 	log ("Nothing found to convert, quitting")
# 		# exit()
# 	return mclean

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
		# exit()

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

