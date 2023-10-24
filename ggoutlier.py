#name:		  	GGOutlier
#created:		August 2023
#by:			paul.kennedy@guardiangeomatics.com
#description:   python module to read a point cloud or raster file of depths, identify outliers, write out points to a shap or laz file for further QC

#requirements
# python 3.10 or newer
# pip install reportlab
# pip install pyproj
# pip install pyshp
# pip install scikit-learn
# pip install rasterio
# pip install numpy #must be 1.26 or better
# pip install matplotlib

#done##########################################
# added support for multiband and singleband tif files.  if multiband, we will look for a band called 'depth' and use that.  if not found we will use the first band.
# trap if zero outliers found
# trap if tif file has no WKT georeferencing
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
# improve how we handle edge of data issues. most of outliers are edge of survey area.
# reduce marker size in pdf report
# improve metrics in pdf report
# improve pdf report to show the regional surface
# move to specific git repository
# spelling errors in report
# write out a prj file alongside the shp file to ESRI standards
# try to improve performance by replacing open3d with numpy
# tidy up the pdf report spelling
# clean memory as we go and use 32 bit floats instead of 64bit.  this helps it a lot!
# epsg code is no longer required.
# need to implement tiling as some tif files are too large.
# replace lastools with pylasfile
# fix word wrapping
# code clean up

# todo ##########################################

import os.path
from argparse import ArgumentParser
from datetime import datetime, timedelta
import math
import numpy as np
import sys
import time
import glob
import rasterio
import multiprocessing as mp
import shapefile
import logging

# locals
import fileutils
import geodetic
import cloud2tif
import ggmbesstandard
import pdfdocument
import pylasfile

###########################################################################
def main():

	iho = ggmbesstandard.sp44()
	msg = str(iho.getordernames())

	parser = ArgumentParser(description='Analyse a floating point TIF file of depths and find all outliers exceeding a user specified threshold.')
	parser.add_argument('-i', 		action='store',			default="", 		dest='inputfile', 		help='(required) Input filename/folder to process.')
	parser.add_argument('-epsg', 	action='store', 		default="",		dest='epsg', 			help='(optional) Specify an output EPSG code for transforming from WGS84 to East,North. If the TIF file is georeferenced this is not required to be specified. e.g. -epsg 32751')
	parser.add_argument('-odir', 	action='store', 		default="",			dest='odir', 			help='(optional) Specify a relative output folder. normally leave this empty.  a folder will be created for you. e.g. -odir GIS')
	parser.add_argument('-near', 	action='store', 		default="5",		dest='near',			help='(optional) ADVANCED:Specify the MEDIAN filter kernel width for computation of the regional surface so nearest neighbours can be calculated. [Default:5]')
	parser.add_argument('-standard',action='store', 		default="order1a",	dest='standard',		help='(optional) Specify the IHO SP44 survey order so we can set the filters to match the required specification. Select from :' + ''.join(msg) + ' [Default:order1a]' )
	parser.add_argument('-unc',		action='store', 		default="",			dest='uncertaintyfilename',		help='(optional) Specify the Uncertainty TIF filename, which is used with the allowable TVU to compute the TVU barometer  [Default:<nothing>]' )
	parser.add_argument('-verbose', 	action='store_true', 	default=False,		dest='verbose',			help='verbose to write LAZ files and other supproting file.s  takes some additional time!,e.g. -verbose [Default:false]')

	matches = []
	args = parser.parse_args()

	if os.path.isfile(args.inputfile):
		matches.append(args.inputfile)

	if len(args.inputfile) == 0:
		# no file is specified, so look for a .pos file in terh current folder.
		inputfolder = os.getcwd()
		matches = fileutils.findFiles2(False, inputfolder, "*.tif")

	if os.path.isdir(args.inputfile):
		matches = fileutils.findFiles2(False, args.inputfile, "*.tif")

	if len(matches) == 0:
		print("oops, no files found to process, quitting")
		exit(0)

	#get the WKT from the TIF file and add it to the args so we can make available....
	WKT = cloud2tif.getWKT(matches[0])
	parser.add_argument('-wkt',		action='store', 		default=WKT,			dest='wkt')
	args = parser.parse_args()
	#load the python proj projection object library if the user has requested it
	if len(args.epsg) == 0:
		args.epsg = str(geodetic.wkt2epsg(wkt=WKT))
	#add the parameters to the args object
	geo = geodetic.geodesy(EPSGCode=args.epsg, wkt=WKT)

	# args.outlierpercentage = min(5.0, float(args.outlierpercentage))
	start_time = time.time() # time the process
	for file in matches:
		#make an output folder
		if len(args.odir) == 0:
			args.odir = os.path.join(os.path.dirname(file), str("GGOutlier_%s" % (time.strftime("%Y%m%d-%H%M%S"))))
		if not os.path.isdir(args.odir):
			args.odir = os.path.join(os.path.dirname(file), args.odir)
		makedirs(args.odir)

		logfilename = os.path.join(args.odir,"GGOutlier_log.txt").replace('\\','/')
		logging.basicConfig(filename = logfilename, level=logging.INFO)
		log("Output Folder: %s" % (args.odir))
		log("EPSGCode for geodetic conversions: %s" % (args.epsg))
		
		log("Configuration: %s" % (str(args)))
		log("GGOutlier Version: 3.01")
		log("GGOutlier started at: %s" % (datetime.now()))
		log("Username: %s" %(os.getlogin()))
		log("Computer: %s" %(os.environ['COMPUTERNAME']))
		log("Number of CPUs %d" %(mp.cpu_count()))	
		log("QC to Survey Standard: %s" % (args.standard))
		iho = ggmbesstandard.sp44()
		standard = iho.loadstandard(args.standard)
		log("Survey_Standard: %s" %(standard.details()))

		#skip the file if its the same as the uncertianty file name
		if os.path.basename(file) == os.path.basename(args.uncertaintyfilename):
			continue

		# args.odir = str("GGOutlier_%s" % (time.strftime("%Y%m%d-%H%M%S")))
		# odir = os.path.join(os.path.dirname(matches[0]), args.odir)
		# makedirs(odir)
		# process(file, args)

		# args.odir = str("GGOutlier_%s_V2" % (time.strftime("%Y%m%d-%H%M%S")))
		# odir = os.path.join(os.path.dirname(matches[0]), args.odir)
		# makedirs(odir)
		process2(file, args)
		log("QC Duration:%.3fs" % (time.time() - start_time))
		log("Creating Report...")
		pdfdocument.GGOutlierreport(logfilename, args.odir)
		log("Report Complete")
		log("QC complete at: %s" % (datetime.now()))	
############################################################
def cleanup(args):
	'''clean up the odir folder '''
	#remove the tif files
	matches = fileutils.findFiles2(False, args.odir, "*.tif")
	for file in matches:
		os.remove(file)

############################################################
def process2(filename, args):
	'''we will try to find outliers using fast array matrix processing.'''

	# if os.path.exists(args.uncertaintyfilename):
	# 	outfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.splitext(os.path.basename(filename))[0] + "_TVU_Barometer.tif")
	# 	standard.computeTVUBarometer(allowabletvufilename, args.uncertaintyfilename, outfilename)
	# 	log ("Created TVU Barometer TIF file for validation: %s " % (outfilename))

	#####################################
	#####################################
	#VERSION 2 of ENGINE to find OUTLIERS
	# tile the file so we are good for memory

	# status, filename = fileutils.copyfile(filename, os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename)))
					
	#get size and resolution...
	log("Processing file: %s" % (filename))
	pixels, SRCRESOLUTION = cloud2tif.getsize(filename)

	#before we tile we need to check if a multiband raster and if so, we need to split it into single band rasters
	#we need to do this because the tileraster function only works on single band rasters
	bandnames = cloud2tif.getbandnames(filename)
	if len(bandnames) == 1:
		bandfilename = filename
	else:
		log("Bands in tif file: %s" % (str(bandnames)))
		for band in bandnames:
			if 'depth' in band.lower():
				outfilename = os.path.join(args.odir, os.path.splitext(os.path.basename(filename))[0] + "_Depth.tif")
				bandfilename = cloud2tif.multibeand2singleband(filename, outfilename, band)
				break
	
	if not os.path.exists(bandfilename):
		log("Oops, no file found to process, please ensure you have provided a single band file or multi-band file with 'depth' and a band name, quitting")
		return ""
	
	log("Processing band file: %s" % (bandfilename))
	#its a single band raster so we can just tile it...
	tilefolder = cloud2tif.tileraster(bandfilename, args.odir, tilewidth = 4096, tileheight = 4096, tileoverlap= 0)
	matches = fileutils.findFiles2(False, tilefolder, "*.tif")

	#now we can process each tile in turn...
	confirmedoutliers = []
	confirmedinliers = []
	ptout = []
	originalfilename = filename
	for tileidx, filename in enumerate(matches):
		depthfilename = filename
		# log("Processing Tile: %s" % (filename))
		#we need to make a regional grid which uses the nearest neighbours, so make this 3 times larger than the source grid.  this means 1 pixel each side of the current point
		regionalfilename = os.path.join(os.path.dirname(filename), os.path.splitext(os.path.basename(filename))[0] + "_RegionalDepth.tif")
		makedirs(os.path.dirname(regionalfilename))
		fname = cloud2tif.smoothtif(filename, regionalfilename, near=int(args.near))
		# log ("Tiled Created REGIONAL TIF file for IHO validation: %s " % (fname))

		# log("QC to Survey Standard: %s" % (args.standard))
		iho = ggmbesstandard.sp44()
		standard = iho.loadstandard(args.standard)
		# log("Survey_Standard: %s" %(standard.details()))
		outfilename = os.path.join(os.path.dirname(filename), os.path.splitext(os.path.basename(filename))[0] + "_TVU_Allowable.tif")
		log ("Computing TVU surface: %s " % (outfilename))
		allowabletvufilename = standard.computeTVUSurface(filename, outfilename)

		deltazfilename = os.path.join(os.path.dirname(filename), os.path.splitext(os.path.basename(filename))[0] + "_DeltaZ.tif")
		log ("Computing DeltaZ surface: %s " % (deltazfilename))
		standard.computeDeltaZ(regionalfilename, depthfilename, deltazfilename)

		# log ("Created DeltaZ TIF file for validation of ALL depths: %s " % (deltazfilename))

		outliersfilename = os.path.join(os.path.dirname(filename),  os.path.splitext(os.path.basename(filename))[0] + "_Outliers.tif")
		log ("Pass1: Finding Outliers...")
		outliersfilename, xydz = standard.findoutliers(allowabletvufilename, deltazfilename, outliersfilename)

		# we can now double check the outliers to see how far they are away from the resulting inlier raster file of mean depths.  
		# if they are close then we can re-accept them
		# log ("Validating candidates against TVU standard")
		rio = rasterio.open(regionalfilename)
		depthio = rasterio.open(depthfilename)

		log ("Pass2: double checking outliers...")
		# EDGE CLEAN UP: check if at edge of data by building a mask for regional and see if anything is an empty cell...			
		radius = 2 #the mask is set to 5 so we use 2 each side of central point.  this means we will skip outlier detection 2 pixels from the edge of the raster
		mask = []
		for x in range(-radius, radius):
			for y in range(-radius, radius):
				mask.append([(x * SRCRESOLUTION), (y * SRCRESOLUTION)])
		npmask = np.asarray(mask)

		for idx, pt in enumerate(xydz):
			npmaskrealworld = np.array(npmask)
			npmaskrealworld[:,0] += pt[0]
			npmaskrealworld[:,1] += pt[1]

			if rio.nodatavals[0] in list(rio.sample(npmaskrealworld)):
			# if rio.nodatavals[0] in list(rio.sample(mask)):
				#this is an edge case so drop it
				confirmedinliers.append(pt)
			else:
				#this is a real outlier so keep it
				confirmedoutliers.append(pt)
				griddepth = next(rio.sample([(pt[0], pt[1])]))[0]
				depth = next(depthio.sample([(pt[0], pt[1])]))[0]
				tvu = standard.gettvuat(griddepth)
				deltaz = abs(pt[2])
				ptout.append([pt[0], pt[1], depth, deltaz, tvu, griddepth])
		
		rio.close()
		depthio.close()
		#tell the user something is happening
		update_progress("Analysing Tiles for Outliers...", tileidx / len(matches))

	log ("Points checked: %s" % (f'{pixels:,}'))
	# log ("Points checked: %s" % (f'{len(xyz):,}'))
	log ("Points outside specification: %s" % (f'{len(ptout):,}'))
	if pixels > 0:
		log ("Percentage outside specification: %.7f" % (100 * (len(ptout)/ pixels)))

	#write the outliers to a point SHAPE file
	shpfilename = os.path.join(args.odir, os.path.splitext(os.path.basename(originalfilename))[0] + "_OutlierPoints" + ".shp")
	# shpfilename = os.path.join(os.path.dirname(originalfilename), args.odir, os.path.basename(originalfilename) + "_OutlierPoints" + ".shp")
	w = shapefile.Writer(shpfilename)
	w.field('DeltaZ', 'N',8,3) # 8 byte floats, 3 decimal places
	w.field('AllowedTVU', 'N',8,3)
	w.field('Depth', 'N',8,3)
	w.field('GridDepth', 'N',8,3)
	w.field('SICApproved','C','254')
	log ("Writing outliers to: %s" % (shpfilename))

	for idx, pt in enumerate(ptout):
		w.pointz(pt[0], pt[1], pt[2])
		w.record(pt[3], pt[4], pt[2], pt[5])
	w.close()
	#write out the prj so it opens nicely in GIS
	cloud2tif.createprj(shpfilename.replace(".shp",".prj"), args.wkt)
	# cloud2tif.createprj(shpfilename.replace(".shp",".prj"), args.epsg, args.wkt)

	#OUTLIERS reporting...
	outfile = os.path.join(args.odir, os.path.splitext(os.path.basename(originalfilename))[0] + "_OutlierPoints" + ".txt")
	# outfile = os.path.join(os.path.dirname(originalfilename), args.odir, os.path.basename(originalfilename) + "_OutlierPoints" + ".txt")
	np.savetxt(outfile, ptout, fmt='%.4f', delimiter=',', newline='\n')
	log ("Created TXT file of outliers: %s " % (outfile))
	#write the outliers to a point cloud laz file
	# fname = lashelper.txt2las(outfile, epsg=args.epsg)

	if len(ptout) > 0 :
		#write to the las file using pylasfile...
		outfile = os.path.join(args.odir, os.path.splitext(os.path.basename(originalfilename))[0] + "_OutlierPoints" + ".las")
		# outfile = os.path.join(os.path.dirname(originalfilename), args.odir, os.path.basename(originalfilename) + "_OutlierPoints" + ".las")
		writer = pylasfile.laswriter(filename=outfile, lasformat=1.4)
		pointsourceID = 1
		writer.hdr.FileSourceID = pointsourceID
		# write out a WGS variable length record so users know the coordinate reference system
		writer.writeVLR_WKT(args.wkt)
		writer.writeVLR_WGS84()
		writer.hdr.PointDataRecordFormat = 1
		a = np.array(ptout)
		# columns = zip(*ptout) #transpose rows to columns
		writer.writepointlist(a[:,0],a[:,1],a[:,2])
		writer.close()	
		log ("Created LAS file of outliers: %s " % (outfile))
	else:
		log ("No outliers found, no las file created.")

	#clean up the tiles...
	cleanup(args)

	log("Creating a regional file for QC purposes")
	try:
		regionalfilename = os.path.join(args.odir, os.path.splitext(os.path.basename(originalfilename))[0] + "_RegionalDepth.tif")
		makedirs(os.path.dirname(regionalfilename))
		fname = cloud2tif.smoothtif(originalfilename, regionalfilename, near=int(args.near))
		log("Creating a regional file for QC purposes: %s" % (fname))
	except:			
		log("Error while creating regional file. Maybe memory is an issue?")

	if args.verbose:
		#load the tif file...	
		with rasterio.open(originalfilename) as src:
			band1 = src.read(1)
			z = band1.flatten() 
			height = band1.shape[0]
			width = band1.shape[1]
			cols, rows = np.meshgrid(np.arange(width), np.arange(height))
			xs, ys = rasterio.transform.xy(src.transform, rows, cols)
			x = np.array(xs).flatten()
			y = np.array(ys).flatten()
			xyz = np.stack((x,y,z), axis=1)
			NODATA = src.nodatavals[0]
			xyz = xyz[np.all(xyz != NODATA, axis=1)]

		src.close()
		log ("Creating LAS file of RAW points...")
		#write to the las file using pylasfile...
		outfile = os.path.join(args.odir, os.path.splitext(os.path.basename(originalfilename))[0] + "_RAWPoints" + ".las")
		writer = pylasfile.laswriter(filename=outfile, lasformat=1.4)
		pointsourceID = 1
		writer.hdr.FileSourceID = pointsourceID
		# write out a WGS variable length record so users know the coordinate reference system
		writer.writeVLR_WKT(args.wkt)
		writer.writeVLR_WGS84()
		writer.hdr.PointDataRecordFormat = 1
		a = np.array(xyz)
		# columns = zip(*ptout) #transpose rows to columns
		writer.writepointlist(a[:,0],a[:,1],a[:,2])
		writer.close()	
		log ("Created LAS file of RAW points: %s " % (outfile))

	log("QC complete at: %s" % (datetime.now()))
	return shpfilename

# ##################################################################################
# def findoutlier(pcd, low, high, TARGET=1.0, NUMPOINTS=3):
# 	'''clean outliers using binary chop to control how many points we reject'''
# 	'''use spherical radius to identify outliers and clusters'''
# 	'''binary chop will aim for target percentage of data deleted rather than a fixed filter level'''
# 	'''this way the filter adapts to the data quality'''
# 	'''TARGET is the percentage of the input points we are looking to reject'''
# 	'''NUMPOINTS is the number of nearest neighbours within the spherical radius which is the threshold we use to consider a point an outlier.'''
# 	'''If a point has no friends, then he is an outlier'''
# 	'''if a point has moew the NUMPOINTS in the spherical radius then he is an inlier, ie good'''

# 	# Force garbage collection
# 	gc.collect()

# 	#outlier removal by radius
# 	# http://www.open3d.org/docs/latest/tutorial/geometry/pointcloud_outlier_removal.html?highlight=outlier
# 	# http://www.open3d.org/docs/latest/tutorial/Advanced/pointcloud_outlier_removal.html
	
# 	#cl: The pointcloud as it was fed in to the model (for some reason, it seems a bit pointless to return this).
# 	#ind: The index of the points which are NOT outliers
# 	currentfilter = (high+low)/2

# 	cl, inlierindex = pcd.remove_statistical_outlier(nb_neighbors=NUMPOINTS,	std_ratio=currentfilter)
# 	# cl, inlierindex = pcd.remove_radius_outlier(nb_points = NUMPOINTS, radius = currentfilter)

# 	inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
# 	outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)
# 	percentage 		= (100 * (len(outlier_cloud.points) / len(pcd.points)))
# 	log ("Current filter StdDev %.2f" % (currentfilter))
# 	log ("Percentage Rejection %.2f" % (percentage))

# 	decimals = len(str(TARGET).split(".")[1])
# 	percentage = round(percentage, decimals)
# 	if percentage < TARGET:
# 		#we have rejected too few, so run again setting the low to the pervious value
# 		log ("Filter level decreasing to reject a few more points...")
# 		del inlier_cloud
# 		del outlier_cloud
# 		del cl 				
# 		del inlierindex
# 		pcd, inlier_cloud, outlier_cloud, inlierindex = findoutlier(pcd, low, currentfilter, TARGET, NUMPOINTS)
# 	elif percentage > TARGET:
# 		#we have rejected too few, so run again setting the low to the pervious value
# 		log ("Filter level increasing to reject a few less points...")
# 		del inlier_cloud
# 		del outlier_cloud
# 		del cl 				
# 		del inlierindex
# 		pcd, inlier_cloud, outlier_cloud, inlierindex = findoutlier(pcd, currentfilter, high, TARGET, NUMPOINTS)

# 	return (pcd, inlier_cloud, outlier_cloud, inlierindex)

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



############################################################
# def process(filename, args):
# 	'''we will try to find outliers using machine learning to determine the typical noise level in the file and adapt filter accordingly so the user-requested noise level is obatined.'''

# 	depthfilename = filename
# 	log("Processing file: %s" % (filename))

# 	log("QC to Survey Standard: %s" % (args.standard))
# 	iho = ggmbesstandard.sp44()
# 	standard = iho.loadstandard(args.standard)
# 	log("Survey_Standard: %s" %(standard.details()))
# 	outfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.splitext(os.path.basename(filename))[0] + "_TVU_Allowable.tif")
# 	allowabletvufilename = standard.computeTVUSurface(filename, outfilename)
# 	if os.path.exists(args.uncertaintyfilename):
# 		outfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.splitext(os.path.basename(filename))[0] + "_TVU_Barometer.tif")
# 		standard.computeTVUBarometer(allowabletvufilename, args.uncertaintyfilename, outfilename)

# 		log("TVU Surface created: %s" % (outfilename))
# 	#load the python proj projection object library if the user has requested it
# 	geo = geodetic.geodesy(args.epsg)
# 	log("EPSGCode for geodetic conversions: %s" % (args.epsg))

# 	pingcounter = 0
# 	beamcountarray = 0
# 	ZSCALE = float(args.zscale) # we might prefer 5 for this as this is how we like to 'look' for spikes in our data.  this value exaggerates the Z values thereby placing more emphasis on the Z than then X,Y

# 	#load the tif file...	
# 	with rasterio.open(filename) as src:
# 		band1 = src.read(1)
# 		z = band1.flatten()
# 		print('Band1 has shape', band1.shape)
# 		height = band1.shape[0]
# 		width = band1.shape[1]
# 		cols, rows = np.meshgrid(np.arange(width), np.arange(height))
# 		xs, ys = rasterio.transform.xy(src.transform, rows, cols)
# 		x = np.array(xs).flatten()
# 		y = np.array(ys).flatten()
# 		print('lons shape', x.shape)
# 		# src._crs.wkt
# 		xyz = np.stack((x,y,z), axis=1)
# 		NODATA = src.nodatavals[0]
# 		SRCRESOLUTION = src.res[0]
	
# 		# Delete the array
# 		del x
# 		del y
# 		del z
# 		del band1
# 		del xs
# 		del ys
# 		# Force garbage collection
# 		gc.collect()
# 		src.close()

# 	log("Loading Point Cloud...")
# 	pcd = o3d.geometry.PointCloud()
# 	#remove the NODATA values
# 	xyz = xyz[np.all(xyz != NODATA, axis=1)]
# 	#scale up the Z data so we accentuate the Z axis compared to the X,Y axis.  this is important as we use 
# 	xyz[:,2] *= ZSCALE
# 	pcd.points = o3d.utility.Vector3dVector(xyz)
# 	xyz[:,2] /= ZSCALE
# 	log("Depths loaded for quality control: %s" % (f'{len(pcd.points):,}'))

# 	if args.verbose:
# 		#RAW report on RAW POINTS
# 		outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.splitext(os.path.basename(filename))[0] + "_RawPoints.txt")
# 		log ("Creating raw laz file of input raw points: %s " % (outfile))
# 		np.savetxt(outfile, xyz, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
# 		fname = lashelper.txt2las(outfile, epsg=args.epsg)
# 		fileutils.deletefile(outfile)

# 	del xyz
# 	# Force garbage collection
# 	gc.collect()

# 	# Populate the 'counter' field automatically so we can track the points accepted/rejected status
# 	beamcountarray = np.arange(0, len(pcd.points))  # This will populate 'counter' with values 1, 2, 3
# 	#PERCENTAGE PASSMARK
# 	log("Learning about your signal to noise levels...")
# 	start_time = time.time() # time the process
# 	low = 0
# 	high = 100
# 	TARGET = float(args.outlierpercentage)
# 	NUMPOINTS = max(int(args.numpoints),1)
# 	#MACHINE LEARNING TO DETERMINE CORRECT LEVEL OF FILTER...
# 	pcd, inlier_cloud, outlier_cloud, inlierindex = findoutlier(pcd, low, high, TARGET, NUMPOINTS)
# 	log ("Points accepted by machine learning: %d" % (len(inlier_cloud.points)))
# 	log ("Points tagged for further evaluation: %d" % (len(outlier_cloud.points)))
# 	# inliers = np.asarray(inlier_cloud.points)
# 	outliers = np.asarray(outlier_cloud.points)
# 	# inliers[:,2] /= ZSCALE
# 	outliers[:,2] /= ZSCALE
# 	########

# 	#we need 1 list of ALL beams which are either accepted or rejected.
# 	beamqualityresult = np.isin(beamcountarray, inlierindex)

# 	#RAW save as a tif file so we can easily view in GIS...
# 	#we need to make a regional grid which uses the nearest neighbours, so make this 3 times larger than the source grid.  this means 1 pixel each side of the current point
# 	regionalfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_RegionalDepth.tif")
# 	fname = cloud2tif.smoothtif(filename, regionalfilename, near=int(args.near))
# 	log ("Created REGIONAL TIF file for IHO validation: %s " % (fname))

# 	if args.verbose:
# 		outfilename = os.path.join(outfile + "_RegionalDepth.laz")
# 		outfilename = lashelper.demzip2(regionalfilename, outfilename, nodata=-9999)
# 		log ("Created REGIONAL LAZ file of input raw points: %s " % (filename))

# 	# we can now double check the outliers to see how far they are away from the resulting inlier raster file of mean depths.  
# 	# if they are close then we can re-accept them
# 	log ("Validating candidates against TVU standard")
# 	rio = rasterio.open(regionalfilename)
# 	outlieridx = 0
# 	confirmedoutliers = []
# 	confirmedinliers = []
# 	for idx, validity in enumerate(beamqualityresult):
# 		if validity == True:
# 			continue
# 		else:
# 			update_progress("Validating Candidates", idx/max(1,len(beamqualityresult)))
# 			#the beam quality is false which means we need to get the next value from the outliers list
# 			pt = outliers[outlieridx]
# 			#the outlier IDX is the sequential number used in conjunction with the beam quality results.  Its not great but thats how the open3d cleaning works
# 			outlieridx = outlieridx + 1
# 			depth = pt[2]
# 			# EDGE CLEAN UP: check if at edge of data by building a mask for regional and see if anything is an empty cell...			
# 			isedge=False
# 			radius = 2 #the median filter is set to 5 so we use 2 each side of central point
# 			mask = []
# 			for x in range(-radius, radius):
# 				for y in range(-radius, radius):
# 					mask.append([pt[0] + (x * SRCRESOLUTION), pt[1] + (y * SRCRESOLUTION)])
# 			for val in rio.sample(mask): 
# 				if val == rio.nodatavals[0]:
# 					griddepth = next(rio.sample([(pt[0], pt[1])]))[0]
# 					tvu = standard.gettvuat(griddepth)
# 					deltaz = abs(griddepth-depth)
# 					pt = np.append(pt, [deltaz, tvu, griddepth])
# 					confirmedinliers.append(pt)
# 					#re-accept the point as it is actually in specification
# 					inlierindex.append(idx)
# 					isedge=True
# 					break

# 			if isedge==True:
# 				continue
# 			#check the depth against the regional surface and the IHO standarON
# 			griddepth = next(rio.sample([(pt[0], pt[1])]))[0]
# 			tvu = standard.gettvuat(griddepth)
# 			deltaz = abs(griddepth-depth)
# 			if griddepth == rio.nodatavals[0]:
# 				#skip point if there is no raster surface value
# 				pt = np.append(pt, [deltaz, tvu, griddepth])
# 				confirmedinliers.append(pt)
# 				#re-accept the point as it is actually in specification
# 				inlierindex.append(idx)
# 				continue
# 			if deltaz >= tvu:
# 				pt = np.append(pt, [deltaz, tvu, griddepth])
# 				confirmedoutliers.append(pt)
# 			else:
# 				pt = np.append(pt, [deltaz, tvu, griddepth])
# 				confirmedinliers.append(pt)
# 				#re-accept the point as it is actually in specification
# 				inlierindex.append(idx)

# 	rio.close()
	
# 	#write the outliers to a point SHAPE file
# 	shpfilename = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_OutlierPoints" + ".shp")
# 	w = shapefile.Writer(shpfilename)
# 	w.field('DeltaZ', 'N',8,3) # 8 byte floats, 3 decimal places
# 	w.field('AllowedTVU', 'N',8,3)
# 	w.field('Depth', 'N',8,3)
# 	w.field('GridDepth', 'N',8,3)
# 	w.field('SICApproved','C','254')
# 	log ("Writing outliers to: %s" % (shpfilename))

# 	for idx, pt in enumerate(confirmedoutliers):
# 		w.pointz(pt[0], pt[1], pt[2])
# 		w.record(pt[3], pt[4], pt[2], pt[5])
# 	w.close()
# 	#write out the prj so it opens nicely in GIS
# 	cloud2tif.createprj(shpfilename.replace(".shp",".prj"), args.epsg)

# 	# q: what is a median filter and how does it work?
# 	# a: https://www.youtube.com/watch?v=VvQ2mo3EXHk
# 	#we need 1 list of ALL beams which are either accepted or rejected.  This is now a revised list folowing Standard validation
# 	inlier_cloud 	= pcd.select_by_index(inlierindex, invert = False)
# 	outlier_cloud 	= pcd.select_by_index(inlierindex, invert = True)
# 	percentage 		= (100 * (len(outlier_cloud.points) / len(pcd.points)))
# 	log ("Points accepted: %s" % (f'{len(inlier_cloud.points):,}'))
# 	log ("Points outside specification: %s" % (f'{len(outlier_cloud.points):,}'))
# 	log ("Percentage outside specification: %.7f" % (percentage))
# 	inliers = np.asarray(inlier_cloud.points)
# 	outliers = np.asarray(outlier_cloud.points)
# 	inliers[:,2] /= ZSCALE
# 	outliers[:,2] /= ZSCALE

# 	if args.verbose:
# 		#INLIERS reporting...
# 		outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_InlierPoints" + ".txt")
# 		np.savetxt(outfile, inliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
# 		outfilename = os.path.join(outfile + "_Depth.tif")
# 		# cloud2tif.saveastif(outfilename, geo, inliers, fill=False)
# 		inlierraster = cloud2tif.point2raster(outfilename, geo, inliers, resolution = 1, bintype='mean', fill=False)
# 		#write the outliers to a point cloud laz file
# 		fname = lashelper.txt2las(outfile, epsg=args.epsg)
# 		lashelper.lasgrid4( fname, outfilename, resolution=1, epsg=args.epsg)
# 		# fileutils.deletefile(outfile)
# 		log ("Created LAZ file of inliers: %s " % (fname))

# 	#OUTLIERS reporting...
# 	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_OutlierPoints" + ".txt")
# 	np.savetxt(outfile, outliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
# 	log ("Created TXT file of outliers: %s " % (outfile))
# 	#write the outliers to a point cloud laz file
# 	fname = lashelper.txt2las(outfile, epsg=args.epsg)
# 	log ("Created LAZ file of outliers: %s " % (fname))

# 	# msg = "GGOutlier complete. outlier_cloud.points
# 	log("QC complete at: %s" % (datetime.now()))
# 	return shpfilename


############################################################
############################################################
############################################################
############################################################
############################################################

	# with rasterio.open(deltazfilename) as deltazsrc:
		# pixels = deltazsrc.height * deltazsrc.width
		# band1 = deltazsrc.read(1)
	# 	z = band1.flatten()
	# 	# print('Band1 has shape', band1.shape)
		# height = band1.shape[0]
		# width = band1.shape[1]
	# 	cols, rows = np.meshgrid(np.arange(width), np.arange(height))
	# 	xs, ys = rasterio.transform.xy(deltazsrc.transform, rows, cols)
	# 	xs = np.float32(xs)
	# 	ys = np.float32(ys)
	# 	x = np.array(xs).flatten()
	# 	y = np.array(ys).flatten()
	# 	print('Array Size', x.shape)
	# 	xyz = np.stack((x,y,z), axis=1)
	# 	# xyz = np.stack((x,y,z), axis=1, dtype=np.float32)
	# 	NODATA = deltazsrc.nodatavals[0]
		# SRCRESOLUTION = deltazsrc.res[0]
	# 	#remove the NODATA values
	# 	xyz = xyz[np.all(xyz != NODATA, axis=1)]
	# deltazsrc.close()
	# #garbage collect
	# gc.collect()	
