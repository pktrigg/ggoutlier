#name:		  	kmallclean
#created:		August 2023
#by:			paul.kennedy@guardiangeomatics.com
#description:   python module to read a Kongsberg KMALL file, create a point cloud, identify outliers, write out a NEW kmall file with flags set

#done
#reading of a kmall file to a point cloud
#pass pcd to open3d
#view pcd file
#find outliers
#save inliers, outliers to a file
#add option to clip on angle

#todo
#rewrite rejected records to a new kmall file
#option to reject n percent of the pcd
#option to use different outlier algorithms
#create tif file from inliers

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
	parser.add_argument('-s', 		action='store', 		default="1",	dest='step', 			help='decimate the data to reduce the output size. [Default: 1]')
	parser.add_argument('-c', 		action='store', 		default="45",	dest='clip', 			help='clip outer beams each side to this max angle. Set to -1 to disable [Default: -1]')
	
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

	print("Loading Point Cloud...")
	pointcloud = kmall.Cpointcloud()

	#create an output file....
	outfilename = fileutils.addFileNameAppendage(filename, "_CLEANED")
	outfileptr = open(outfilename, 'wb')

	r = kmall.kmallreader(filename)

	if args.epsg == '0':
		approxlongitude, approxlatitude = r.getapproximatepositon()
		args.epsg = geodetic.epsgfromlonglat (approxlongitude, approxlatitude)

	#load the python proj projection object library if the user has requested it
	geo = geodetic.geodesy(args.epsg)
	print("EPSGCode for geodetic conversions: %s" % (args.epsg))
	recordcount, starttimestamp, enftimestamp = r.getRecordCount()

	# demonstrate how to load the navigation records into a list.  this is really handy if we want to make a trackplot for coverage
	start_time = time.time() # time the process
	print("Modifying Flags...")
	while r.moreData():
		# read a datagram.  If we support it, return the datagram type and aclass for that datagram
		# The user then needs to call the read() method for the class to undertake a fileread and binary decode.  This keeps the read super quick.
		typeofdatagram, datagram = r.readDatagram()
		bytes = datagram.loadbytes() # get a hold of the bytes for the ping so we can modify them and write to a new file.
		if typeofdatagram == '#MRZ':
			datagram.read()
			# clip the outer beams...
			if clip > 0:
				clipper(datagram, clip)

			x, y, z, q = computebathypointcloud(datagram, geo)
			pointcloud.add(x, y, z, q) # pkpkpk
			update_progress("Extracting Point Cloud", counter/recordcount)
			counter = counter + 1

			#write out the kmall datagrem with modified beam flags
			# for beam in datagram.beams:
			# 	#beam flag offset is 7 bytes into the beam structure so we can now set that flag to whatever we want it to be
			# 	bytes [beam.beambyteoffset +7] = 1
			# 	# now write out the modified byte array
			# 	outfileptr.write(bytes)
		else:
			outfileptr.write(bytes)

		if counter == 1000:
			break
		continue

	print("")
	r.close()

	outfile = os.path.join(os.path.dirname(filename), os.path.basename(filename) + "_R.txt")
	xyz = np.column_stack([pointcloud.xarr,pointcloud.yarr, pointcloud.zarr])
	# xyz[:,2] *= 10.0
	print("Saving point cloud to %s" % (outfile)) 
	# print("Point count to %d" % (len(xyz))) 
	np.savetxt(outfile, (xyz), fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')

	pcd = o3d.geometry.PointCloud()
	pcd.points = o3d.utility.Vector3dVector(xyz)

	obb = pcd.get_oriented_bounding_box()
	obb.color = (0,0,0)

	print("Statistical outlier removal")
	voxel_down_pcd = pcd.voxel_down_sample(voxel_size=0.0002)
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=3.0) # 1.51
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=10, std_ratio=3.0) # 1.89
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=10, std_ratio=1.0) 
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=2.0) # 3.54%
	# cl, ind = voxel_down_pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=1.0) # 9.56

	#outlier removal by radius
	# http://www.open3d.org/docs/latest/tutorial/geometry/pointcloud_outlier_removal.html?highlight=outlier
	nb_points=3
	radius=0.5
	cl, ind = voxel_down_pcd.remove_radius_outlier(nb_points= nb_points, radius=radius)

	inlier_cloud = voxel_down_pcd.select_by_index(ind, invert=False)
	outlier_cloud = voxel_down_pcd.select_by_index(ind, invert=True)
	print (inlier_cloud)
	print (outlier_cloud)
	print ("Percentage rejection %.2f" % (100 * (len(outlier_cloud.points) / len(xyz))))

	outfile = os.path.join(os.path.dirname(filename), os.path.basename(filename) + "_C_Radius" + str(nb_points) + "_" + str(radius) + ".txt")
	# inlier_cloud = voxel_down_pcd.select_by_index(ind, invert=False)
	np.savetxt(outfile, (np.asarray(inlier_cloud.points)), fmt='%.2f, %.3f, %.4f', delimiter=',', newline='\n')
	
	outfilename = os.path.join(outfile + ".tif")
	saveastif(outfilename, geo, inlier_cloud)


	# display_inlier_outlier(voxel_down_pcd, ind)

	# o3d.visualization.draw_geometries([pcd, obb])

	# pc = open3d.io.read_point_cloud(outfile, format='xyz')
	print (pcd)
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

	return


###############################################################################
def saveastif(outfilename, geo, inlier_cloud, resolution=1):

	pcd = np.asarray(inlier_cloud.points)
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
	#gt = (X_topleft, X_resolution, 0, Y_topleft, 0, Y_resolution)
	# transform = Affine.translation(xmin, ymax)
	# transform = rasterio.Affine(xmin, xres, 0, ymax, 0, -yres)
	
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
	) 
	# populate the numpy array with the values....
	arr = np.zeros((height+2, width+2), dtype=float)
	for row in pcd:
		px = math.floor(xmax - row[0])
		py = math.floor(ymax - row[1])
		# py, px = src.index(row[0], row[1])
		arr[py, width - px] = row[2]
	
	src.write(arr, 1)

	src.close()

###############################################################################
def clipper(datagram, clip):
	'''using the datagram, reject if the take off angle is outside the clip limit'''

	for beam in datagram.beams:
		if abs(beam.beamAngleReRx_deg) > clip:
			beam.detectionType = 2 # reject the beam
			# beam.rejectionInfo1 = 

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
