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

#todo##########################################
#validate each outlier against the results and re-approve if it is now acceptable
#scale the Z values so we accentuate the outlier noise from the horizontal noise

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
from rasterio.transform import Affine

import kmall
import fileutils
import geodetic

###########################################################################
def main():

	parser = ArgumentParser(description='Read a KMALL file.')
	parser.add_argument('-epsg', 	action='store', 		default="0",	dest='epsg', 			help='Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326')
	parser.add_argument('-i', dest='inputFile', action='store', default="", help='Input filename.pos to process.')
	parser.add_argument('-c', 		action='store', 		default="-1",	dest='clip', 			help='clip outer beams each side to this max angle. Set to -1 to disable [Default: -1]')
	
	files = []
	args = parser.parse_args()
	# args.inputFile = "/Users/paulkennedy/Documents/development/sampledata/0822_20210330_091839.kmall"
	args.inputFile = "c:/sampledata/EM304_0002_20220406_122446.kmall"
	# args.inputFile = "c:/sampledata/EM2040_0822_20210330_091839.kmall"
	args.inputFile = "c:/sampledata/0494_20210530_165628.kmall"
	args.inputFile = "C:/sampledata/kmall/B_S2980_3005_20220220_084910.kmall"
	if len (args.inputFile) == 0:
		# no file is specified, so look for a .pos file in terh current folder.
		inputfolder = os.getcwd()
		files = findFiles2(False, inputfolder, "*.kmall")
	else:
		files.append(args.inputFile)

	for file in files:
		print ("processing file: %s" % (file))
		kmallcleaner(file, args)

############################################################
def kmallcleaner(filename, args):
	'''we will try to auto clean beams by extracting the beam xyzF flag data and attempt to clean in scipy'''
	'''we then set the beam flags to reject files we think are outliers and write the kmall file to a new file'''
	
	counter = 0
	clip = float(args.clip)
	beamcounter = 0

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
	recordcount, starttimestamp, enftimestamp = r.getRecordCount()

	# demonstrate how to load the navigation records into a list.  this is really handy if we want to make a trackplot for coverage
	start_time = time.time() # time the process
	print("Modifying Flags...")
	while r.moreData():
		# read a datagram.  If we support it, return the datagram type and aclass for that datagram
		# The user then needs to call the read() method for the class to undertake a fileread and binary decode.  This keeps the read super quick.
		typeofdatagram, datagram = r.readDatagram()
		counter = counter + 1
		if typeofdatagram == '#MRZ':
			datagram.read()
			x, y, z, q = computebathypointcloud(datagram, geo)
			pointcloud.add(x, y, z, q)
			update_progress("Extracting Point Cloud", counter/recordcount)

		# if counter == 1000:
		# 	break
		# continue

	print("")
	r.close()

	outfile = os.path.join(os.path.dirname(filename), os.path.basename(filename) + "_R.txt")
	xyz = np.column_stack([pointcloud.xarr,pointcloud.yarr, pointcloud.zarr])
	xyz[:,2] *= 10.0

	pcd = o3d.geometry.PointCloud()
	pcd.points = o3d.utility.Vector3dVector(xyz)

	# outfilename = os.path.join(outfile + "_R.tif")
	# saveastif(outfilename, geo, pcd)

	#lets clean the data to a user specified threshold using the input data quality to control the filter.  this means the machine learns about the data...
	########
	low = 0
	high = 10
	target = 1.5
	pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier1(pcd, low, high, target)
	########

	# outfile = os.path.join(os.path.dirname(filename), os.path.basename(filename) + "_C_Inlier" + ".txt")
	np.savetxt(outfile, (np.asarray(inlier_cloud.points)), fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	outfilename = os.path.join(outfile + ".tif")
	inlierraster = saveastif(outfilename, geo, inlier_cloud, fill=True)

	#we can now revalidate the outliers and re-accept if they fit the surface
	# outlier_cloud = validateoutliers(inlierraster, outlier_cloud)

	outfile = os.path.join(os.path.dirname(filename), os.path.basename(filename) + "_C_Outlier" + ".txt")
	# np.savetxt(outfile, (np.asarray(outlier_cloud.points)), fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	outfilename = os.path.join(outfile + ".tif")
	# saveastif(outfilename, geo, outlier_cloud, fill=False)

	#now lets write out a NEW KMALL file with the beams modified...
	#create an output file....
	outfilename = fileutils.addFileNameAppendage(filename, "_CLEANED")
	outfileptr = open(outfilename, 'wb')

	print("Writing NEW KMALL file %s" % (outfilename))
	counter = 0
	r = kmall.kmallreader(filename)
	while r.moreData():
		# read a datagram.  If we support it, return the datagram type and aclass for that datagram
		# The user then needs to call the read() method for the class to undertake a fileread and binary decode.  This keeps the read super quick.
		typeofdatagram, datagram = r.readDatagram()
		bbytes = datagram.loadbytes() # get a hold of the bytes for the ping so we can modify them and write to a new file.
		counter = counter + 1
		if typeofdatagram == '#MRZ':
			datagram.read()
			# clip the outer beams...
			if clip > 0:
				clipper(datagram, clip)

			#apply the results of the cleaning process...
			setbeamquality(datagram, beamcounter, inlierindex)
			update_progress("Writing cleaned data", counter/recordcount)

			#write out the kmall datagrem with modified beam flags
			barray=bytearray(bbytes)
			for beam in datagram.beams:
				# beam flag offset is 3 bytes into the beam structure so we can now set that flag to whatever we want it to be
				barray [beam.beambyteoffset + 3] = beam.detectionType
				# barray [beam.beambyteoffset + 5] = beam.detectionType
				# barray [beam.beambyteoffset + 7] = beam.detectionType
				# barray [beam.beambyteoffset + 8] = beam.detectionType
				beamcounter += 1
			# now write out the modified byte array
			outfileptr.write(bytes(barray))

		else:
			outfileptr.write(bbytes)

		# if counter == 1000:
		# 	break
		# continue
	return

###############################################################################
def setbeamquality(datagram, beamcounter, inlierindex):
	'''apply the cleaning results to the ping of data'''
	test_set=set(inlierindex)
	for idx, beam in enumerate(datagram.beams):
		if not beamcounter+idx in test_set: 
		# if not beamcounter+idx in inlierindex: 
			beam.detectionType = 2
		else:
			pk=1
	return
###############################################################################
def validateoutliers(inlierraster, outlier_cloud):

	pcd = np.asarray(outlier_cloud.points)
	for row in pcd:
		# py, px = inlierraster.index(row[0], row[1])
		v = inlierraster.sample(row[0], row[1])

	# py, px = inlierraster.index(row[0], row[1])

	return outlier_cloud

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


##################################################################################
def cleanoutlier1(pcd, low, high, target):
	'''clean outliers using binary chop to control how many points we reject'''
	'''use spherical radius to identify outliers and clusters'''
	'''binary chop will aim for target percentage of data deleted rather than a fixed filter level'''
	'''this way the filter adapts to the data quality'''

	currentfilter = (high+low)/2

	#outlier removal by radius
	# http://www.open3d.org/docs/latest/tutorial/geometry/pointcloud_outlier_removal.html?highlight=outlier
	# http://www.open3d.org/docs/latest/tutorial/Advanced/pointcloud_outlier_removal.html
	
	nb_points=3 # the number points inside 
	radius=currentfilter
	#cl: The pointcloud as it was fed in to the model (for some reason, it seems a bit pointless to return this).
	#ind: The index of the points which are NOT outliers
	cl, inlierindex = pcd.remove_radius_outlier(nb_points= nb_points, radius=radius)

	inlier_cloud = pcd.select_by_index(inlierindex, invert=False)
	outlier_cloud = pcd.select_by_index(inlierindex, invert=True)
	# print (inlier_cloud)
	# print (outlier_cloud)
	percentage = (100 * (len(outlier_cloud.points) / len(pcd.points)))
	print ("Percentage rejection %.2f" % (percentage))

	percentage = round(percentage, 1)
	if percentage < target:
		#we have rejected too few, so run again setting the low to the pervious value
		print ("Filter level increasing to reject a few more points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier1(pcd, low, currentfilter, target)
		# percentage = cleanoutlier1(pcd, low, currentfilter, target)
	elif percentage > target:
		#we have rejected too few, so run again setting the low to the pervious value
		print ("Filter level decreasing to reject a few less points...")
		pcd, inlier_cloud, outlier_cloud, inlierindex = cleanoutlier1(pcd, currentfilter, high, target)
		# percentage = cleanoutlier1(pcd, currentfilter, high, target)
	# else:
	return (pcd, inlier_cloud, outlier_cloud, inlierindex)

###############################################################################
def saveastif(outfilename, geo, cloud, resolution=1, fill=False):

	if len(cloud.points)==0:	
		return
		
	NODATA = -999
	pcd = np.asarray(cloud.points)
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
	from rasterio.transform import from_origin
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
	arr = np.full((height+1, width+1), fill_value=NODATA, dtype=float)
	
	from numpy import ma
	arr = ma.masked_values(arr, NODATA)

	for row in pcd:
		px = math.floor((xmax - row[0]) / xres)
		py = math.floor((ymax - row[1]) / yres)
		# py, px = src.index(row[0], row[1])
		arr[py, px] = row[2]
		
	#we might want to fill in the gaps. useful sometimes...
	if fill:
		from rasterio.fill import fillnodata
		arr = fillnodata(arr, mask=None, max_search_distance=xres*2, smoothing_iterations=0)

	src.write(arr, 1)
	src.close()
	# return src

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
		beam.depth = beam.z_reRefPoint_m + datagram.txTransducerDepth_m
		# beam.depth = beam.z_reRefPoint_m - datagram.z_waterLevelReRefPoint_m
		
	npeast = np.fromiter((beam.east for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	npnorth = np.fromiter((beam.north for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	npdepth = np.fromiter((beam.depth for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)
	npq = np.fromiter((beam.rejectionInfo1 for beam in datagram.beams), float, count=len(datagram.beams)) #. Also, adding count=len(stars)

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


###############################################################################
def findFiles2(recursive, filespec, filter):
	'''tool to find files based on user request.  This can be a single file, a folder start point for recursive search or a wild card'''
	matches = []
	if recursive:
		matches = glob(os.path.join(filespec, "**", filter), recursive = True)
	else:
		matches = glob(os.path.join(filespec, filter))
	
	mclean = []
	for m in matches:
		mclean.append(m.replace('\\','/'))
		
	# if len(mclean) == 0:
	# 	print ("Nothing found to convert, quitting")
		# exit()
	return mclean

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
if __name__ == "__main__":
		main()
		# exit()
