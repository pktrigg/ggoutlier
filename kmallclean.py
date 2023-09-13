#name:		  	kmallclean
#created:		August 2023
#by:			paul.kennedy@guardiangeomatics.com
#description:   python module to read a Kongsberg KMALL file, create a point cloud, identify outliers, write out a NEW kmall file with flags set

#done##########################################
#reading of a kmall file to a point cloud
#pass pcd to open3d
#view pcd file
#find outliers
#save inliers, outliers to a file
#add option to clip on angle
#create tif file from inliers
#option to reject n percent of the pcd
#create tif file of raw data
#create a tif file of inliers
#create a tif file of outliers
#optionally fill the tif file to interpolate.  we need this for the revalidation
#rewrite rejected records to a new kmall file
#added percentage to args
#added numpoints to args
#added debug to args
#fixed descaling after cleaning to txt file export
#write outliers to a shape file point cloud so we can visualise them easily in GIS
#profile to improve performance
#scale the Z values so we accentuate the outlier noise from the horizontal noise
#improve the indexing to tif file so its faster
#code clean up
#point clouds save to laz files so we can import and view in qgis neat profile view
#trap max numpoints to be 1 or more
#ignore filter for outer beams at edge of swath

#todo##########################################
#test with 10X vertical
#test with 50x vertical
#validate each outlier against the results and re-approve if it is now acceptable

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

import kmall
import fileutils
import geodetic
import multiprocesshelper
import cloud2tif
import lashelper

###########################################################################
def main():

	parser = ArgumentParser(description='Read a KMALL file.')
	parser.add_argument('-epsg', 	action='store', 		default="0",	dest='epsg', 			help='Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326')
	parser.add_argument('-i', 		action='store',			default="", 	dest='inputfile', 		help='Input filename/folder to process.')
	parser.add_argument('-c', 		action='store', 		default="-1",	dest='clip', 			help='clip outer beams each side to this max angle. Set to -1 to disable [Default: -1]')
	parser.add_argument('-cpu', 	action='store', 		default='0', 	dest='cpu', 			help='number of cpu processes to use in parallel. [Default: 0, all cpu]')
	parser.add_argument('-odir', 	action='store', 		default="",	dest='odir', 			help='Specify a relative output folder e.g. -odir GIS')
	parser.add_argument('-n', 		action='store', 		default="3",	dest='numpoints', 		help='Specify the number of nearest neighbours points to use.  More points means more data will be rejected. ADVANCED ONLY [Default:5]')
	parser.add_argument('-p', 		action='store', 		default="1.0",	dest='outlierpercentage',help='Specify the approximate percentage of data to remove.  the engine will analyse the data and learn what filter settings are appropriate for your waterdepth and data quality. This is the most important (and only) parameter to consider spherical radius to find the nearest neightbours. [Default:1.0]')
	parser.add_argument('-z', 		action='store', 		default="1.0",	dest='zscale',			help='Specify the ZScale to accentuate the depth difference ove the horizontal distance between points. Think of this as how you exxagerate the vertical scale in a swath editor to more easily spot the outliers. [Default:1.0]')
	parser.add_argument('-debug', 	action='store', 		default="-1",	dest='debug', 			help='Specify the number of pings to process.  good only for debugging. [Default:-1]')
	
	matches = []
	args = parser.parse_args()
	# args.inputfile = "C:/sampledata/kmall/B_S2980_3005_20220220_084910.kmall"

	if os.path.isfile(args.inputfile):
		matches.append(args.inputfile)

	if len (args.inputfile) == 0:
		# no file is specified, so look for a .pos file in terh current folder.
		inputfolder = os.getcwd()
		matches = fileutils.findFiles2(False, inputfolder, "*.kmall")

	if os.path.isdir(args.inputfile):
		matches = fileutils.findFiles2(False, args.inputfile, "*.kmall")

	if len(args.odir) == 0:
		args.odir = str("NearestNeighbours_%d_OutliersPercent_%.2f" % (int(args.numpoints), float(args.outlierpercentage)))

	#make an output folder
	if (len(matches) > 0):
		odir = os.path.join(os.path.dirname(matches[0]), args.odir)
		print("Output Folder: %s" % (odir))
		makedirs(odir)

	results = []
	if args.cpu == '1':
		for file in matches:
			kmallcleaner(file, args)
	else:
		multiprocesshelper.log("Files to Import: %d" %(len(matches)))		
		cpu = multiprocesshelper.getcpucount(args.cpu)
		print("Processing with %d CPU's" % (cpu))

		pool = mp.Pool(cpu)
		multiprocesshelper.g_procprogress.setmaximum(len(matches))
		poolresults = [pool.apply_async(kmallcleaner, (file, args), callback=multiprocesshelper.mpresult) for file in matches]
		pool.close()
		pool.join()
		# for idx, result in enumerate (poolresults):
		# 	results.append([file, result._value])
		# 	print (result._value)

############################################################
def kmallcleaner(filename, args):
	'''we will try to auto clean beams by extracting the beam xyzF flag data and attempt to clean in scipy'''
	'''we then set the beam flags to reject files we think are outliers and write the kmall file to a new file'''
	
	maxpings = int(args.debug)
	if maxpings == -1:
		maxpings = 999999999

	pingcounter = 0
	clip = float(args.clip)
	beamcountarray = 0
	ZSCALE = float(args.zscale) # we might prefer 5 for this as this is how we like to 'look' for spikes in our data.  this value exaggerates the Z values thereby placing more emphasis on the Z than then X,Y
	
	print("Loading Point Cloud...")
	pointcloud = kmall.Cpointcloud()

	r = kmall.kmallreader(filename)

	if args.epsg == '0':
		approxlongitude, approxlatitude = r.getapproximatepositon()
		args.epsg = geodetic.epsgfromlonglat (approxlongitude, approxlatitude)

	#load the python proj projection object library if the user has requested it
	geo = geodetic.geodesy(args.epsg)
	print("EPSGCode for geodetic conversions: %s" % (args.epsg))
	
	#get the record count so we can show a progress bar
	recordcount, starttimestamp, enftimestamp = r.getRecordCount("MRZ")

	# demonstrate how to load the navigation records into a list.  this is really handy if we want to make a trackplot for coverage
	while r.moreData():
		# read a datagram.  If we support it, return the datagram type and aclass for that datagram
		# The user then needs to call the read() method for the class to undertake a fileread and binary decode.  This keeps the read super quick.
		typeofdatagram, datagram = r.readDatagram()
		if typeofdatagram == '#MRZ':
			datagram.read()
			x, y, z, q = computebathypointcloud(datagram, geo)
			pointcloud.add(x, y, z, q)
			update_progress("Extracting Point Cloud", pingcounter/recordcount)
			pingcounter = pingcounter + 1

		if pingcounter == maxpings:
			break

	print("")
	r.close()

	pcd = o3d.geometry.PointCloud()
	xyz = np.column_stack([pointcloud.xarr,pointcloud.yarr, pointcloud.zarr])
	xyz[:,2] *= ZSCALE
	pcd.points = o3d.utility.Vector3dVector(xyz)
	print("Depths loaded for cleaning: %s" % (f'{len(pcd.points):,}'))

	# Populate the 'counter' field automatically
	beamcountarray = np.arange(0, len(pcd.points))  # This will populate 'counter' with values 1, 2, 3

	#lets clean the data to a user specified threshold using the input data quality to control the filter.  this means the machine learns about the data...
	########
	print("Understanding your data noise levels...")
	start_time = time.time() # time the process
	low = 0
	high = 10
	TARGET = float(args.outlierpercentage)
	NUMPOINTS = max(int(args.numpoints),1)
	pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier(pcd, low, high, TARGET, NUMPOINTS)
	print ("Points accepted: %.2f" % (len(inlier_cloud.points)))
	print ("Points rejected: %.2f" % (len(outlier_cloud.points)))
	inliers = np.asarray(inlier_cloud.points)
	outliers = np.asarray(outlier_cloud.points)
	inliers[:,2] /= ZSCALE
	outliers[:,2] /= ZSCALE
	print("Clean Duration: %.3f seconds" % (time.time() - start_time)) # print the processing time. It is handy to keep an eye on processing performance.
	########

	#we need 1 list of ALL beams which are either accepted or rejected.
	beamqualityresult = np.isin(beamcountarray, inlierindex)

	#report on RAW POINTS
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_R.txt")
	xyz[:,2] /= ZSCALE
	np.savetxt(outfile, xyz, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	fname = lashelper.txt2las(outfile)
	fileutils.deletefile(outfile)
	print ("Created LAZ file of input raw points: %s " % (fname))
	# outfilename = os.path.join(outfile + "_R.tif")
	# raw = np.asarray(pcd.points)
	# raw[:,2] /= ZSCALE
	# cloud2tif.saveastif(outfilename, geo, raw, fill=False)

	#report on INLIERS
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_Inlier" + ".txt")
	np.savetxt(outfile, inliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	# outfilename = os.path.join(outfile + ".tif")
	# inlierraster = cloud2tif.saveastif(outfilename, geo, inliers, fill=True)
	#write the outliers to a point cloud laz file
	fname = lashelper.txt2las(outfile)
	fileutils.deletefile(outfile)
	print ("Created LAZ file of inliers: %s " % (fname))

	#report on OUTLIERS
	outfile = os.path.join(os.path.dirname(filename), args.odir, os.path.basename(filename) + "_Outlier" + ".txt")
	np.savetxt(outfile, outliers, fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	#write the outliers to a point cloud laz file
	fname = lashelper.txt2las(outfile)
	fileutils.deletefile(outfile)
	print ("Created LAZ file of outliers: %s " % (fname))

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

	print("Writing NEW KMALL file %s" % (outfilename))
	pingcounter = 0
	beamcounter = 0

	r = kmall.kmallreader(filename)
	while r.moreData():
		# read a datagram.  If we support it, return the datagram type and aclass for that datagram
		# The user then needs to call the read() method for the class to undertake a fileread and binary decode.  This keeps the read super quick.
		typeofdatagram, datagram = r.readDatagram()
		bbytes = datagram.loadbytes() # get a hold of the bytes for the ping so we can modify them and write to a new file.
		if typeofdatagram == '#MRZ':
			datagram.read()
			# clip the outer beams...
			if clip > 0:
				clipper(datagram, clip)

			update_progress("Writing cleaned data", pingcounter/recordcount)
			pingcounter = pingcounter + 1

			#write out the kmall datagrem with modified beam flags
			barray=bytearray(bbytes)
			for idx, beam in enumerate(datagram.beams):
				#apply the results of the cleaning process...

				if not beamqualityresult[beamcounter]:
					if (idx > 0) & (idx < len(datagram.beams)): # do not reject the outer edge of the swath.  these ar
						beam.detectionType = 2
				# beam flag offset is 3 bytes into the beam structure so we can now set that flag to whatever we want it to be
				barray [beam.beambyteoffset + 3] = beam.detectionType
				#this is the beam counter for ALL pings + current ping since start of file.  we use it to keep track of the cleaned points so this is super important.
				beamcounter += 1
			# now write out the modified byte array
			outfileptr.write(bytes(barray))

		else:
			outfileptr.write(bbytes)

		if pingcounter == maxpings:
			break
		# continue
	return

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
	print ("Current filter radius %.2f" % (currentfilter))
	print ("Percentage rejection %.2f" % (percentage))

	decimals = len(str(TARGET).split(".")[1])
	percentage = round(percentage, decimals)
	if percentage < TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		print ("Filter level increasing to reject a few more points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier(pcd, low, currentfilter, TARGET, NUMPOINTS)
	elif percentage > TARGET:
		#we have rejected too few, so run again setting the low to the pervious value
		print ("Filter level decreasing to reject a few less points...")
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
	print (inlier_cloud)
	print (outlier_cloud)
	print ("Percentage rejection %.2f" % (100 * (len(outlier_cloud.points) / len(inlier_cloud.points))))
	print("Showing outliers (red) and inliers (gray): ")
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
    
# 	print("EPS: %f MinSample: %f Rejected: %d Survivors: %d InputCount %d" % (eps,  min_samples, len(rejected_points), len(filtered_points), len(points)))
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
# 	# 	print ("Nothing found to convert, quitting")
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
if __name__ == "__main__":
		main()
		# exit()

	########################v#######################################################
	# print("Statistical outlier removal")
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

