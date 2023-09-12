#name:		  	delivershared.py
#created:		jan 2020
#by:			paul.kennedy@guardiangeomatics.com
#description:   python module of shared methods for various things!
#copyright		Guardian Geomatics Pty Ltd
#				This software is explicitly prohibited by use of any non-guardian employee or subcontractor.
# 
##################
# #DONE
##################
# 09/05/2020 
# add support for multiple cores during tiling process, so the oad of tileing perm mosaic is spread across many cores.
# gdalinfo now uses UUID to permit multicore processing

import sys
import stat
import time
import os
import tempfile
import fnmatch
import math
import re
import logging
import shutil
import csv
import ctypes
from datetime import datetime
from datetime import timedelta
import time
from glob import glob
import shlex
import subprocess
import multiprocessing
import uuid

import lashelper
import gdalhelper

errorlog = []

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'shared'))
import shapefile
import geodetic
import SSDM
import pyproj
import dgnwrite
import fileutils

###############################################################################
# def main():
# 	filename 		= "C:/projects/XTF2LAZTEST/SASI_M228_20191015_035823_A096100VA0090.XTF"
# 	outfilename 	= "C:/projects/XTF2LAZTEST/SASI_M228_20191015_035823_A096100VA0090.laz"
# 	resolution 		= 5
# 	epsg 			='31984'
# 	iclippercentage = 50
# 	oclippercentage = 95
# 	xtf2laz(filename, outfilename, resolution, iclippercentage, oclippercentage, epsg)

# 	injectposition2xtf(filename, positionfilename, outfilename)

###############################################################################
def getcpucount(requestedcpu):
	'''control how many CPU's we use for multi processing'''
	if int(requestedcpu) == 0:
		requestedcpu = multiprocessing.cpu_count()

		stat = MEMORYSTATUSEX()
		ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
		# print("MemoryLoad: %d%%" % (stat.dwMemoryLoad))
		# print("MemoryAvailable: %d%%" % (stat.ullAvailPhys/(1024*1024*1024)))
		availablememoryingigs = stat.ullAvailPhys/(1024*1024*1024)
		# make sure we have enough memory per CPU
		requiredgigspercpu = 4

		maxcpu = max(1, int(availablememoryingigs/ requiredgigspercpu))
		# ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
		# print("MemoryLoad: %d%%" % (stat.dwMemoryLoad))
		# ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
		# print("MemoryLoad: %d%%" % (stat.dwMemoryLoad))
		# ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
		# print("MemoryLoad: %d%%" % (stat.dwMemoryLoad))
		# ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
		# print("MemoryLoad: %d%%" % (stat.dwMemoryLoad))

		if int(requestedcpu) > maxcpu:
			requestedcpu = maxcpu
	return int(requestedcpu)

###############################################################################
def stringexistsinlist(mylist, mystring, mystring2 = ""):
	result = None
	for name in mylist:
		if mystring.strip().lower() in name.strip().lower() and mystring2.strip().lower() in name.strip().lower():
			return name.strip().lower()
	return result
	
###############################################################################
def	makedirs(odir):
	if not os.path.isdir(odir):
		os.makedirs(odir, exist_ok=True)
	odirlog = os.path.join(odir, "log").replace('\\','/')
	if not os.path.isdir(odirlog):
		os.makedirs(odirlog)
	return odirlog

###############################################################################
def	mosaiclaser2(args, rawfiles, tilename, odir, prefix, resamplemode, resolution = 1):
	'''mosaic the LAZ files into a multiple rasters with smaller sub tile sizes'''
	#for each subtile, clip the input files.
	#we then run blast on each individual file to make a raster file (LAZ)
	#we then warp the files together into a single mosaic tif using 
	tilerectangle = loadtilerectangle(tilename)

	ysize = tilerectangle.top - tilerectangle.bottom
	xsize = tilerectangle.right - tilerectangle.left
	subtilecount = max(int(args.stc), 1)
	subtilexsize = xsize / subtilecount
	subtileysize = ysize / subtilecount

	log("mosaic started at %s" % (datetime.now()), False, False)
	#odirlog = makedirs(odir)	
	
	count = 0
	for tileidx in range(subtilecount):
		for tileidy in range(subtilecount):
			log("processing subtile %d of %d " % (count, subtilecount*subtilecount), False)
			p1 = POINT(tilerectangle.left+(tileidx*subtilexsize), tilerectangle.bottom+(tileidy*subtileysize)) 
			p2 = POINT(tilerectangle.left+((tileidx+1)*subtilexsize), tilerectangle.bottom+((tileidy+1)*subtileysize)) 
			rect = RECT(p1,p2)
			stime = datetime.now()

			root = os.path.splitext(tilename)[0]	
			# outfilename = os.path.abspath(os.path.join(odir, prefix + str(len(rawfiles)) + '.tif')).upper()
			outfilename = os.path.abspath(os.path.join(odir, prefix + os.path.basename(root) + "_R" + str(tileidx) + "_C" + str(tileidy) + '.tif')).upper()
			outfilename = outfilename.replace('\\','/')

			if os.path.exists(outfilename):				
				# Handle errors while calling os.remove()
				try:
					os.remove(outfilename)
				except:
					log("Error while deleting file %s, skipping mosaic creation." % (outfilename), True)
					return

			srcfiles = ""
			#odirsubtilemosaic = os.path.join(odir, "subtilemosaic")
			odirclip = os.path.join(odir, "subtile_%d_%d" %(tileidx, tileidy)).replace('\\','/')
			if not os.path.isdir(odirclip):
				os.makedirs(odirclip)

			for num, f in enumerate(rawfiles):
				root = os.path.splitext(f.clippedfilename)[0]	
				testfilename = os.path.join(root + '_G_C.laz').replace('\\','/')
				fileexists = findFiles2(True, odir, os.path.basename(testfilename))
				#only make the clip if the file has not been previoulsy processed
				if len(fileexists) == 0:
					srcfiles = f.clippedfilename + " " + srcfiles
					lashelper.lasclipbb(f.filename, rect, odirclip, resolution, nodata=0, prefix="")
					update_progress("subtiling and gridding file: %d/%d " % (num, len(rawfiles)), num/len(rawfiles))
				else:
					update_progress("skipping existing file: %d/%d " % (num, len(rawfiles)), num/len(rawfiles))

			#get the list of the clipped files for the subtile
			clippedfiles = findFiles2(False, odirclip, "*.laz")		
			if len(clippedfiles) > 0:
				#print(clippedfiles)
				#there are some files, so we can merge, remove duplicates and make a tif...
				filespec = os.path.join(odirclip, "*.laz").replace('\\','/')
				fname = os.path.basename(odirclip) + ".laz"
				mergedfilename = lashelper.lasmerge(filespec, fname, odirclip)
				mergedfilename = lashelper.lasduplicate(mergedfilename)

				#now we have 1 laz file for the tile, we can make the tif and XYZ

				gridfilename = lashelper.lasgrid(mergedfilename, resolution)
				log ("Gridded: %s to %s" % (mergedfilename, gridfilename))
				gdalhelper.pyramid(gridfilename, odir)
				hillshadefilename = lashelper.hillshade(mergedfilename, odirclip, resolution)
				log ("QC Hillshade: %s to %s" % (gridfilename, hillshadefilename))
				gdalhelper.pyramid(hillshadefilename, odir)
				#las2asc(gridfilename)
					
			else:
				log("INFO:subtile folder exists. will use files already in place", False)

			#now we need to merge the candidate files into a single raster.  this can be done with gdalwarp or lasmerge followed by lasduplicate
			#shutil.rmtree(odirclip)

			count = count + 1
		log ("Mosaic Duration: %s" % (datetime.now() - stime))
	return outfilename


###############################################################################
def calculateboundaries(inputfolder, outputfolder, srcfilespec, resolution, fill, epsg):
	'''compute the boundaries area of the input files, and merge into a single output shape file.'''

	boundary = []
	files = findFiles2(False, inputfolder, srcfilespec)
	for filename in files:
		# print ("Coverage: Computing Boundary Polygon for file: %s" %(filename))
		#we need to clip out small degenerate laz files which create invalid boundaries.  make sure they are more than 1*1m bounding box
		rect = lashelper.getlazboundingbox(filename, os.path.dirname(filename))
		if (rect.top-rect.bottom < 1) or (rect.right-rect.left < 1):
			continue
		root = os.path.splitext(filename)[0]
		outfilename = os.path.join(outputfolder, root+"_coverage.shp").replace('\\','/')
		result = lashelper.lasboundary(filename, outfilename, 0, resolution)

#		result = lasboundaries(filename, outputfolder, resolution, fill, epsg, prefix=root+"_coverage")
		polygondetails2 = getpolgonareas(result)

		for idx, p  in enumerate(polygondetails2):
			p.name = os.path.basename(root)
			boundary = boundary + [p]

	###############################################################################
	# WRITE TO SHAPEFILE USING PYSHP
	shapewriter = shapefile.Writer()
	shapewriter.field("ID")
	shapewriter.field("NAME")
	shapewriter.field("AREA", 'N', decimal=3)
	for idx, feature in enumerate(boundary):
		# step1: convert shapely to pyshp using the function above
		converted_shape = shapely_to_pyshp(feature.geometry)
		# step2: tell the writer to add the converted shape
		shapewriter._shapes.append(converted_shape)
		# add a list of attributes to go along with the shape
		shapewriter.record(str(idx), feature.name, feature.geometry.area)

	# save it
	if len(shapewriter.shapes()) > 0:
		coveragefilename = os.path.join(outputfolder, "BATHYMETRY_Boundaries.shp")
		shapewriter.save(coveragefilename)
		prjfilename = coveragefilename.replace('.shp','.prj')
		geodetic.writePRJ(prjfilename, epsg)
		return coveragefilename

###############################################################################
# ###############################################################################
# def calculatecoverage(filesforgriddingfolder, poly, resolution, fill, epsg, outputfolder=""):
# 	'''compute the gaps in the file so we can estimate our coverage rate'''

# 	#we cannot compute gaps in data outside the survey area, so clip each input file to the survey area.  
# 	gaps = []
# 	tiledir = ""

# 	print ("Coverage: Clipping data so we only measure gaps inside the survey area...") 
# 	for filename in filesforgriddingfolder:
# 		clipdir = os.path.join(os.path.dirname(filename),"clip").replace('\\','/')
# 		os.makedirs(clipdir, exist_ok=True)
# 		if len(poly) > 0:
# 			lasclip(filename, poly, clipdir, 0, "", False)
# 		else:	
# 			copyfile(filename, os.path.join(clipdir, os.path.basename(filename)).replace('\\','/'), True)

# 	# as an alternative strategy we tile into smaller chunks, compute the gaps within each tile and then concatenate all the results.  
# 	# this is not perfect but scales.
# 	print ("Coverage: Tile into smaller chunks, compute the gaps within each tile and concatenate all the results. this scales...") 
# 	for filename in filesforgriddingfolder:
# 		tiledir = os.path.join(os.path.dirname(filename),"tile").replace('\\','/')
# 		break

# 	input = os.path.abspath(clipdir).replace('\\','/')+"/*.laz"
# 	tiledir = lastile(input, tiledir, tile_size=5000, prefix="")
# 	tilefiles = findFiles2(False, tiledir, "*.laz")
# 	for filename in tilefiles:
# 		print ("Coverage: Computing gaps for tile: %s" %(filename))
# 		#we need to clip out small degenerate laz files which create invalid boundaries.  make sure they are more than 1*1m bounding box
# 		rect = getlazboundingbox(filename, os.path.dirname(filename))
# 		if (rect.top-rect.bottom < 1) or (rect.right-rect.left < 1):
# 			continue
# 		root = os.path.splitext(filename)[0]
# 		result = lasboundaries(filename, outputfolder, resolution, fill, epsg, prefix=root+"_coverage")
		
# 		polygondetails2 = getpolgonareas(result)

# 		#if there is only 1 polygon it is the entire area, so pop it.
# 		if len(polygondetails2) == 1:
# 			polygondetails2.pop(0)

# 		parentpolygons = findparentpolygon(polygondetails2, resolution)
# 		for idx, p  in enumerate(polygondetails2):
# 			if not idx in parentpolygons:
# 				gaps = gaps + [p]

# 		#we only need to do something if there are polygons...
# 		# if len(polygondetails2) > 0:
# 		# 	#remove the first polygon as it will be the polygon over the entire area.
# 		# 	polygondetails2.sort(key=lambda x: x.area, reverse=True)
# 		# 	polygondetails2.pop(0)
# 		# 	gaps = gaps + polygondetails2

# 		# for outerpoly in polygondetails2:
# 		# 	outergeometry = outerpoly.geometry.simplify(float(resolution)*3).buffer(0)
# 		# 	toberejected = False
# 		# 	if outergeometry.area == 0:
# 		# 		toberejected = True
# 		# 		continue
# 		# 	for innerpoly in polygondetails2:
# 		# 		innergeometry = innerpoly.geometry.simplify(float(resolution)*3).buffer(0)
# 		# 		if outergeometry.area == innergeometry.area:
# 		# 			# print ("skipping self-reference")
# 		# 			continue
# 		# 		if outergeometry.contains(innergeometry):
# 		# 			toberejected = True
# 		# 			continue

# 		# 	if toberejected == False:
# 		# 		gaps = gaps + [outerpoly]
# 		# 	else:
# 		# 		print ("removing poly with area %.3f" % (outergeometry.area))


# 	###############################################################################
# 	# WRITE TO SHAPEFILE USING PYSHP
# 	shapewriter = shapefile.Writer()
# 	shapewriter.field("ID")
# 	shapewriter.field("AREA", 'N', decimal=3)
# 	for idx, feature in enumerate(gaps):
# 		# step1: convert shapely to pyshp using the function above
# 		converted_shape = shapely_to_pyshp(feature.geometry)
# 		# step2: tell the writer to add the converted shape
# 		shapewriter._shapes.append(converted_shape)
# 		# add a list of attributes to go along with the shape
# 		shapewriter.record(str(idx), feature.geometry.area)

# 	# save it
# 	if len(shapewriter.shapes()) > 0:
# 		coveragefilename = os.path.join(outputfolder, "BATHYMETRY_Coverage_Gaps.shp")
# 		shapewriter.save(coveragefilename)
# 		prjfilename = coveragefilename.replace('.shp','.prj')
# 		geodetic.writePRJ(prjfilename, epsg)

# 	# compute the outline around the tif file and all the holes
# 	# input = os.path.abspath(clipdir).replace('\\','/')+"/*.laz"
# 	# result = lasboundaries(input, outputfolder, resolution=1, prefix="coverage")
# 	# polygondetails2 = getpolgonareas(result)
# 	# polygondetails2.sort(key=lambda x: x.area, reverse=True)
# 	# polygondetails2.pop(0)

# 	if len(gaps) == 0:
# 		log ("there are NO GAPS in the survey area so a shape file of gaps cannot be produced")

# 	holestotal = 0
# 	for i, p in enumerate(gaps):
# 		log ("Gap %d Area: %.5fm²" %(i, p.area))
# 		holestotal = holestotal + p.area

# 	if len(poly) > 0:
# 		surveyareapolygondetails = getpolgonareas(poly)
# 		#print the results...
# 		for i, p in enumerate(surveyareapolygondetails):
# 			print ("SurveyPolygon: Area: %.5fm²" %(p.area))
# 		#compute the 
# 		print ("Sum of all holes: %.5fm²" %(holestotal))
# 		for i, p in enumerate(surveyareapolygondetails):
# 			if p.area > 0:
# 				print ("Gaps as a percentage of entire Survey Area: %.4f " %(holestotal/p.area*100))

###############################################################################
def cleanshapefile(filename, epsg, resolution=0):
	'''clean a shapefile to remove duplicate points from gemetry'''
	#https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python

	if not os.path.exists(filename):
		print ("oops, shape file does not exist %s" % (filename))
		return	
	root, ext = os.path.splitext(filename)
	outputfolder = os.path.dirname(filename)
	boundary = []

	polygondetails2 = getpolgonareas2(filename, resolution)
	polygondetails3 = removeinternalpolygons(polygondetails2)
	for idx, p  in enumerate(polygondetails3):
		p.name = os.path.basename(root)
		boundary = boundary + [p]


	###############################################################################
	# WRITE TO SHAPEFILE USING PYSHP
	shapewriter = shapefile.Writer()
	shapewriter.field("ID")
	shapewriter.field("NAME")
	shapewriter.field("AREA", 'N', decimal=3)
	for idx, feature in enumerate(boundary):
		# step1: convert shapely to pyshp using the function above
		converted_shape = shapely_to_pyshp(feature.geometry)
		# step2: tell the writer to add the converted shape
		shapewriter._shapes.append(converted_shape)
		# add a list of attributes to go along with the shape
		shapewriter.record(str(idx), feature.name, feature.geometry.area)

	# save it
	if len(shapewriter.shapes()) > 0:
		coveragefilename = os.path.join(outputfolder, "MYcleaned_Boundaries.shp")
		shapewriter.save(coveragefilename)
		prjfilename = coveragefilename.replace('.shp','.prj')
		geodetic.writePRJ(prjfilename, epsg)
		return coveragefilename


	# #remember the path.  shapely changes it!
	# mypath = os.environ['PATH']

	# import shapely
	# shapelypath = os.path.dirname(os.path.abspath(shapely.__file__))
	# os.environ["PATH"] += os.pathsep + shapelypath+"\\DLLs"
	# # from shapely.geometry.polygon import Polygon
	# from shapely.geometry import shape

	# details=[]
	# if not os.path.exists(filename):
	# 	return details

	# e = shapefile.Editor(filename)

	# # e.shapes()[0].__geo_interface__
	# # for s in e.shapes(): #Iterate through shapes in shapefile
	# # 	f = s.__geo_interface__
	# # 	#get hold of the area.  this is important
	# # 	shp_geom = shape(f) 
	# 	# shp_geom = shp_geom.buffer(0)
	# 	# if shp_geom.area > 0:
	# 		# converted_shape = shapely_to_pyshp(shp_geom)
	# 		# s = converted_shape
	# # for s in enumerate(e._shapes):
	# # 		geom = e._shapes[i]
	# # 		geom = geom.geometry.buffer(0)
	# # 		# Change the y value to 5
	# # 		geom.points[0][1] = 5
	# # 		# Change the z value to 9
	# # 		if hasattr(geom, "z"):
	# # 			geom.z = (9,)
	# # 		else:
	# # 			geom.points[0][2] = 9
	
	# root, ext = os.path.splitext(filename + "pk")
	# e.save(root)


	# sf = shapefile.Reader(filename)
	# #catch empty shape files
	# if sf.numRecords == 0:
	# 	return details

	# # shapes = sf.shapes()
	# # for shape in shapes:
	# sf.shapes()[0].__geo_interface__
	# for s in sf.shapes(): #Iterate through shapes in shapefile
	# 	f = s.__geo_interface__
	# 	#get hold of the area.  this is important
	# 	shp_geom = shape(f) 
	# 	p = POLYGON("", shp_geom.area, shp_geom)
	# 	if shp_geom.area > 0:
	# 		details.append(p)

	# #reset the path to how it used to be
	# os.environ['PATH'] = mypath
	# return details


	# parentidx = []
	# for polya in polygondetails2:
	# 	geometrya = polya.geometry.buffer(0)
	# 	if geometrya.area == 0:
	# 		continue
	# 	for idx, polyb in enumerate(polygondetails2):
	# 		geometryb = polyb.geometry.buffer(0)
	# 		if geometrya.area == geometryb.area:
	# 			continue
	# 		if geometrya.within(geometryb):
	# 			if not idx in parentidx:
	# 				#if it is inside a polygon, its a hole to mark the parent for deletion
	# 				parentidx.append(idx)
	# return parentidx

###############################################################################
def removeinternalpolygons(polygondetails2):
	'''test each polygon and if there is a parent, set its status so it can be killed off.  these are the holes in the data'''
	externalpolygons = []
	for polya in polygondetails2:
		geometrya = polya.geometry.buffer(0)
		if geometrya.area == 0:
			continue
		for idx, polyb in enumerate(polygondetails2):
			geometryb = polyb.geometry.buffer(0)
			if geometrya.area == geometryb.area:
				continue
			if geometrya.within(geometryb):
				polya.isinsideanotherpolygon=True
	
	for poly in polygondetails2:
		#now only export out the external polygons
		if poly.isinsideanotherpolygon == False:
			p = POLYGON("", poly.geometry.area, poly.geometry)
			externalpolygons.append(p)
	return externalpolygons

###############################################################################
def findparentpolygon(polygondetails2, resolution=1):
	'''test each polygon and if there is a parent, return the ID so it can be killed off'''
	parentidx = []
	for polya in polygondetails2:
		geometrya = polya.geometry.buffer(0)
		if geometrya.area == 0:
			continue
		for idx, polyb in enumerate(polygondetails2):
			geometryb = polyb.geometry.buffer(0)
			if geometrya.area == geometryb.area:
				continue
			if geometrya.within(geometryb):
				if not idx in parentidx:
					#if it is inside a polygon, its a hole to mark the parent for deletion
					parentidx.append(idx)
	return parentidx

###############################################################################
def getgaps(filename):
	'''compute the areas of a polygon from a shape file which are gaps.  these are dfined as interior rings and are clockwise in nature'''
	#https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python

	#remember the path.  shapely changes it!
	mypath = os.environ['PATH']

	import shapely
	shapelypath = os.path.dirname(os.path.abspath(shapely.__file__))
	os.environ["PATH"] += os.pathsep + shapelypath+"\\DLLs"
	# from shapely.geometry.polygon import Polygon
	from shapely.geometry import shape
	from shapely import geometry

	details=[]
	if not os.path.exists(filename):
		return details

	sf = shapefile.Reader(filename)
	#catch empty shape files
	if sf.numRecords == 0:
		return details

	# sf.shapes()[0].__geo_interface__
	for s in sf.shapes(): #Iterate through shapes in shapefile
		f = s.__geo_interface__
		#a polygon object can contain many parts, typically the inner rings.  split up the polygon into the individual polygons so ew can delete the largest one.
		for thispart in f["coordinates"]:
				pl = shapely.geometry.Polygon(thispart)
				p = POLYGON("", pl.area, pl)
				if p.area > 0:
					details.append(p)

	#reset the path to how it used to be
	os.environ['PATH'] = mypath

	#sort the polygons into ascending orderby size
	details.sort(key=lambda x: x.area, reverse=False)

	#delete the largest one, is it is the outline around entire survey
	print("extrior polygon represents survey area.  Size is: %.1f" % (details[-1].area))
	details.pop()
	return details

###############################################################################
# def getpolgonareas2(filename, resolution=0):
# 	'''compute the areas of a polygon from a shape file'''
# 	#https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python

# 	#remember the path.  shapely changes it!
# 	mypath = os.environ['PATH']

# 	import shapely
# 	from shapely.validation import explain_validity

# 	shapelypath = os.path.dirname(os.path.abspath(shapely.__file__))
# 	os.environ["PATH"] += os.pathsep + shapelypath+"\\DLLs"
# 	# from shapely.geometry.polygon import Polygon
# 	from shapely.geometry import shape

# 	details=[]
# 	if not os.path.exists(filename):
# 		return details

# 	sf = shapefile.Reader(filename)
# 	#catch empty shape files
# 	if sf.numRecords == 0:
# 		return details

# 	# shapes = sf.shapes()
# 	# for shape in shapes:
# 	sf.shapes()[0].__geo_interface__
# 	for s in sf.shapes(): #Iterate through shapes in shapefile
# 		f = s.__geo_interface__
# 		#get hold of the area.  this is important
# 		shp_geom = shape(f) 
# 		# shp_geomclean = shp_geom.buffer(0)	
# 		validity = explain_validity(shp_geom) # 'Self-intersection[234107.500043162 581106.999994076]'
# 		if validity.lower() != "valid geometry":
# 			points=shapely.geometry.MultiPoint(shp_geom.exterior.coords[1:])
# 			shp_geom=points.convex_hull
			
# 		shp_geomclean = shp_geom.simplify(resolution).buffer(0)

# 		print(explain_validity(shp_geom))

# 		p = POLYGON("", shp_geomclean.area, shp_geomclean)
# 		if shp_geomclean.area > 0:
# 			details.append(p)

# 	#reset the path to how it used to be
# 	os.environ['PATH'] = mypath
# 	return details

###############################################################################
def getpolgonareas(filename):
	'''compute the areas of a polygon from a shape file'''
	#https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python

	#remember the path.  shapely changes it!
	mypath = os.environ['PATH']

	import shapely
	shapelypath = os.path.dirname(os.path.abspath(shapely.__file__))
	os.environ["PATH"] += os.pathsep + shapelypath+"\\DLLs"
	# from shapely.geometry.polygon import Polygon
	from shapely.geometry import shape

	details=[]
	if not os.path.exists(filename):
		return details

	sf = shapefile.Reader(filename)
	#catch empty shape files
	if sf.numRecords == 0:
		return details

	# shapes = sf.shapes()
	# for shape in shapes:
	sf.shapes()[0].__geo_interface__
	for s in sf.shapes(): #Iterate through shapes in shapefile
		f = s.__geo_interface__
		#get hold of the area.  this is important
		shp_geom = shape(f) 
		p = POLYGON("", shp_geom.area, shp_geom)
		if shp_geom.area > 0:
			details.append(p)

	#reset the path to how it used to be
	os.environ['PATH'] = mypath
	return details

###############################################################################
def shapely2shpfile(POLYGON, filename, epsg="4326"):
	# WRITE TO SHAPEFILE USING PYSHP
	shapewriter = shapefile.Writer()
	# shapewriter.field("ID","C", 255)
	shapewriter.field("dataimg",	 	"C", size=255)					#SSDMSurvey_TrackLines class
	shapewriter.field("fileimg", 		"C", size=255)					#SSDMSurvey_TrackLines class
	for feature in POLYGON:
		# step1: convert shapely to pyshp using the function above
		converted_shape = shapely_to_pyshp(feature.geometry)
		# step2: tell the writer to add the converted shape
		shapewriter._shapes.append(converted_shape)
		# add a list of attributes to go along with the shape
		# root, ext = os.path.splitext(feature.name)
		
		#replace // and the /./ by normalising
		fname = os.path.normpath(feature.name.upper())
		#get the file name and extension from the base.
		root, ext = os.path.splitext(os.path.basename(fname))
		#the client wants the clipped files rathe than the original filed (this is a crime), so provide this...
		root = root + "_clipped"
		#get the foldername
		dirname = os.path.dirname(fname)
		strdataimg = os.path.basename(root).upper() + ext
		strfileimg = os.path.normpath(os.path.join(dirname, "CLIPPED", root.upper()+ext.upper()))
		shapewriter.record(strdataimg, strfileimg)
		# shapewriter.record(os.path.basename(feature.name).upper(), feature.name.upper())

	# save it
	if len(shapewriter.shapes()) > 0:
		shapewriter.save(filename)
		prjfilename = filename.replace('.SHP','.PRJ')
		geodetic.writePRJ(prjfilename, epsg)


###############################################################################
def string_escape(s, encoding='utf-8'):
	return (s.encode('latin1')         # To bytes, required by 'unicode-escape'
		.decode('unicode-escape') # Perform the actual octal-escaping decode
		.encode('latin1')         # 1:1 mapping back to bytes
		.decode(encoding))  

# THIS FUNCTION CONVERTS A GEOJSON GEOMETRY DICTIONARY TO A PYSHP SHAPE OBJECT
###############################################################################
def shapely_to_pyshp(shapelygeom):
	# first convert shapely to geojson
	try:
		shapelytogeojson = shapely.geometry.mapping
	except:
		import shapely.geometry
		shapelytogeojson = shapely.geometry.mapping
	geoj = shapelytogeojson(shapelygeom)
	# create empty pyshp shape
	record = shapefile._Shape()
	# set shapetype
	if geoj["type"] == "Null":
		pyshptype = 0
	elif geoj["type"] == "Point":
		pyshptype = 1
	elif geoj["type"] == "LineString":
		pyshptype = 3
	elif geoj["type"] == "Polygon":
		pyshptype = 5
	elif geoj["type"] == "MultiPoint":
		pyshptype = 8
	elif geoj["type"] == "MultiLineString":
		pyshptype = 3
	elif geoj["type"] == "MultiPolygon":
		pyshptype = 5
	record.shapeType = pyshptype
	# set points and parts
	if geoj["type"] == "Point":
		record.points = geoj["coordinates"]
		record.parts = [0]
	elif geoj["type"] in ("MultiPoint","Linestring"):
		record.points = geoj["coordinates"]
		record.parts = [0]
	elif geoj["type"] in ("Polygon"):
		record.points = geoj["coordinates"][0]
		record.parts = [0]
	elif geoj["type"] in ("MultiPolygon","MultiLineString"):
		index = 0
		points = []
		parts = []
		for eachmulti in geoj["coordinates"]:
			points.extend(eachmulti[0])
			parts.append(index)
			index += len(eachmulti[0])
		record.points = points
		record.parts = parts
	return record

# ###############################################################################
# def appendtoshapefile(filename, shape, record=""):
# 	# We must include a file extension in
# 	# this case because the file name
# 	# has multiple dots and pyshp would get
# 	# confused otherwise.
# 	# Create a shapefile reader
# 	r = shapefile.Reader(filename)
# 	# Create a shapefile writer
# 	# using the same shape type
# 	# as our reader
# 	w = shapefile.Writer(r.shapeType)
# 	# Copy over the existing dbf fields
# 	w.fields = list(r.fields)
# 	# Copy over the existing dbf records
# 	w.records.extend(r.records())
# 	# Copy over the existing polygons
# 	w._shapes.extend(r.shapes())

# 	shapewriter._shapes.append(converted_shape)

# 	# Add a new polygon
# 	##w.poly(parts=[[[-104,24],[-104,25],[-103,25],[-103,24],[-104,24]]])
# 	# Add a new dbf record for our polygon making sure we include
# 	# all of the fields in the original file (r.fields)
# 	##w.record("STANLEY","TD","091022/1500","27","21","48","ep")
# 	# Overwrite the old shapefile or change the name and make a copy 
# 	w.save(file_name)


###############################################################################
def tif2outline(filename, odir, resolution, nodata, warp, prefix, replace=False):
		'''read a tif file and create a text file polygon representation around it using lastools'''
		root = os.path.splitext(filename)[0]
		if os.path.exists(root+".txt"):
			return root+".txt"
	
		#convert the tif files to single band files with a lesser resolution so they are quick to make...
		result1 = gdalhelper.createsingleband(filename, odir, resolution, nodata, warp, prefix, replace)

		#convert the tif files to laz so we can make boundary polygons efficiently...
		result2 = lashelper.demzip(result1, odir, nodata, prefix, replace)
		
		#create the polygon txt files...
		outfilename = os.path.join(odir, prefix + root + '_boundary.txt').replace('\\','/')
		result3 = lashelper.lasboundary(result2, outfilename, nodata, resolution, replace)

		#give the output a sensible name
		root = os.path.splitext(os.path.basename(filename))[0]
		outlinefilename = os.path.join(os.path.dirname(result3), root+".txt").replace('\\','/')
		
		try:
			deletefile(outlinefilename)
			os.rename(result3, outlinefilename)
		except:			
			log("Error while renaming file %s" % (result3), True)

		#clean up
		try:
			os.remove(result1)
			os.remove(result2)
		except:			
			log("Error while deleting file %s" % (result1), True)

		return outlinefilename

###############################################################################
def tif2shp(filename, odir, resolution, nodata, warp, prefix, replace=False):
		'''read a tif file and create a text file polygon representation around it using lastools'''
		# root, ext = os.path.splitext(filename)
		# if os.path.exists(root+".shp"):
		# 	return root+".shp"
	
		#convert the tif files to single band files with a lesser resolution so they are quick to make...
		result1 = gdalhelper.createsingleband(filename, odir, resolution, nodata, warp, prefix, replace)

		#convert the tif files to laz so we can make boundary polygons efficiently...
		result2 = lashelper.demzip(result1, odir, nodata, prefix, replace)
		
		#create the polygon txt files...
		outfilename = os.path.join(odir, prefix + '.shp').replace('\\','/')
		result3 = lashelper.lasboundary(result2, outfilename, nodata, resolution, replace)

		return result3


###############################################################################
def txt2shp(args, outfilename, sourcefiles):

	if args.rebuild:
		deleteexistingshapefile(outfilename)
	else:	
		if os.path.exists(outfilename):
			log ("shapefile exists, skipping: %s" % (outfilename))
			return

	shp = SSDM.createCoveragePetrobrasSHP(outfilename)

	polygonspershape = 10000
	for idx, f in enumerate(sourcefiles):
		boundary2shp(f, shp, outfilename, "")
		update_progress("exporting to shape files: %d/%d " % (idx, len(sourcefiles)), idx/len(sourcefiles))

		if idx % polygonspershape == 0:
			if len(shp.records) > 0:
				saveshp(args, shp, outfilename)
				outfilename = fileutils.createOutputFileName(outfilename, ext="")
				shp = SSDM.createCoveragePetrobrasSHP(outfilename)

	print ("Shapefile records created: %d" % (len(shp.records)))
	if len(shp.records) > 0:
		saveshp(args, shp, outfilename)

###############################################################################
def loadtilerectangle(tilename):
	tilefile = TILEFILE("tiles.csv")
	for tile in tilefile.tiles:
		if tile.name.lower() in tilename.lower() :
			return tile.rectangle
	return None

###############################################################################
def clip(sourcefiles, tilename, odir, prefix, resamplemode, args):
	'''clip the files to the tile rectangle '''
	clippedfiles = []
	log("clip started for tilename %s" % (tilename), False, True)
	tilerectangle = loadtilerectangle(tilename)

	# would be good if we can multi process this!!!!  just need to figure out how to remember the output filename
	# if int(args.cpu) == 0:
	# 	print("Number of cpu process : ", multiprocessing.cpu_count())
	# 	args.cpu = multiprocessing.cpu_count()
	# # cpu = multiprocessing.cpu_count()
	# pool = multiprocessing.Pool(max(int(args.cpu),1))
	# log ("def clip MP processing with %d cpu's" % (max(int(args.cpu),1)), False, False)

	# for raster in sourcefiles:
	# 	if '_clipped' in raster.filename:
	# 		raster.clippedfilename = raster.filename
	# 		clippedfiles.append(raster.filename)
	# 		#continue
	# 	else:
	# 		root = os.path.splitext(raster.filename)[0]
	# 		clipfilename = os.path.abspath(os.path.join(odir, os.path.basename(root) + prefix + '.tif')).replace('\\','/')
	# 		raster.clippedfilename = clipfilename
	# 		pool.apply_async(cliptotile, args=(raster.filename, raster.clippedfilename, tilerectangle, odir, root, resamplemode, args, False) )
	# pool.close()
	# pool.join()

	for raster in sourcefiles:		
		#source has a reference to a clipped file, so we do not need to clip it again...
		if '_clipped' in raster.filename:
			raster.clippedfilename = raster.filename
			clippedfiles.append(raster.filename)
			continue
		else:
			#clip the file...
			root = os.path.splitext(raster.filename)[0]
			raster.clippedfilename = os.path.abspath(os.path.join(odir, os.path.basename(root) + prefix + '.tif')).replace('\\','/')
			clipfilename = gdalhelper.cliptotile(raster.filename, raster.clippedfilename, tilerectangle, odir, prefix, resamplemode, args, False)
			raster.clippedfilename = clipfilename

	# log("clipping finished at %s" % (datetime.now()), False)
	# clippedfiles = []
	# for f in sourcefiles:
	# 	clippedfiles.append(f.clippedfilename)
	return sourcefiles

###############################################################################
def saveshp(args, shp, outfilename):
	#log ("Saving coverage: %s" % outfilename)
	shp.save(outfilename)
	# now write out a prj file so the data has a spatial Reference
	filename = outfilename.lower().replace('.shp','.prj').upper()
	geodetic.writePRJ(filename, args.epsg)
	# if args.dgn:
	# 	dgnwrite.convert2DGN(trackPointFileName)
	# subprocess.run(['open', outfilename], check=True)
	# if args.openresult:
	# 	os.startfile(outfilename, 'open')
###############################################################################
def getsegmentfromsgyfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	words = root.split("-")
	segment = "SEGM"
	if len(words) > 1:
		segment = words[1]
	return segment

###############################################################################
def getlinenamefromsgyfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	words = root.split("-")
	result = "NAME"
	if len(words) > 3:
		result = words[2].upper()
		linename, linenumber = linenamesplit(words[2].upper())
		if len(linenumber) > 0:
			result = "%s%03d" % (linename, int(linenumber))
		else:
			result = linename
	return result

###############################################################################
def deleteunrequiredfiles(rootfolder, filter, recursive="True"):
	'''recursively delete files which are not required'''
	matches = findFiles2(recursive, rootfolder, filter)
	for filename in matches:
		if os.path.isfile(filename):
			deletefile(filename)
		# try:
		# 	os.remove(filename)
		# except OSError:
		# 	print("Error while deleting file")

###############################################################################
def deletefile(filename):
	if os.path.exists(filename):
		try:			
			os.remove(filename)
		except:	
			return
			#log("file is locked, cannot delete: %s " % (filename))

###############################################################################
def deleteexistingshapefile(outfilename):
	
	if os.path.exists(outfilename):
		try:			
			os.remove(outfilename)
		except:			
			log("file is locked, cannot delete: %s " % (outfilename))

	filename = outfilename.lower().replace(".shp", ".dbf")
	if os.path.exists(filename):
		try:			
			os.remove(filename)
		except:			
			log("file is locked, cannot delete: %s " % (filename))

	filename = outfilename.lower().replace(".shp", ".prj")
	if os.path.exists(filename):
		try:			
			os.remove(filename)
		except:			
			log("file is locked, cannot delete: %s " % (filename))

	filename = outfilename.lower().replace(".shp", ".shx")
	if os.path.exists(filename):
		try:			
			os.remove(filename)
		except:			
			log("file is locked, cannot delete: %s " % (filename))

###############################################################################
def boundary2shp(raster, shp, shpfilename, tiffile):
	'''create a shape file witht the correct attribute link in the shapefile'''
	
	#multiple polygons are separaeted like this...
	# 388425.5,7514390.5,130
	# 388444.5,7514429.5,107
	# 388448.5,7514437.5,130
	# #
	# 388349.5,7514234.5,110
	# 388350.5,7514234.5,97
	# 388351.5,7514234.5,101

	if not os.path.exists(raster.outlinefilename):
		return

	poly = []

	f = open(raster.outlinefilename, 'r')
	for line in f:
		if line.startswith('#'):
			if len(poly) > 0:
				#write out the polygon
				if len(poly) > 3:
					writepolygon(shp, poly, os.path.abspath(raster.filename))
					#writepolygon(shp, poly, os.path.abspath(raster.outfilename), raster.sortorder)
				poly = []
		else:
			words = line.split(',')
			x = float(words[0])
			y = float(words[1])
			poly.append([x,y])

	f.close()

	return
###############################################################################
def boundary2shpv2(filename, shp, shpfilename, tiffile):
	'''create a shape file with the correct attribute link in the shapefile'''
	
	#multiple polygons are separaeted like this...
	# 388425.5,7514390.5,130
	# 388444.5,7514429.5,107
	# 388448.5,7514437.5,130
	# #
	# 388349.5,7514234.5,110
	# 388350.5,7514234.5,97
	# 388351.5,7514234.5,101

	if not os.path.exists(filename):
		return

	poly = []

	f = open(filename, 'r')
	for line in f:
		if line.startswith('#'):
			if len(poly) > 0:
				#write out the polygon
				if len(poly) > 3:
					writepolygon(shp, poly, os.path.abspath(tiffile))
					#writepolygon(shp, poly, os.path.abspath(raster.outfilename), raster.sortorder)
				poly = []
		else:
			words = line.split(',')
			x = float(words[0])
			y = float(words[1])
			poly.append([x,y])

	f.close()

	return

###############################################################################
def area(coords):
	t=0
	for count in range(len(coords)-1):
		y = coords[count+1][1] + coords[count][1]
		x = coords[count+1][0] - coords[count][0]
		z = y * x
		t += z
	
	#print(abs(t/2.0))
	return abs(t/2.0)
#a=[(5.09,5.8), (1.68,4.9), (1.48,1.38), (4.76,0.1), (7.0,2.83), (5.09,5.8)]
#print _area_(a)

###############################################################################
def writepolygon(coveragePoly, poly, fileName, idx=1):
	
	if len(poly) == 0:
		return

	if area(poly) < 2:
		return

	parts = []

	parts.append(poly)
	coveragePoly.poly(parts=parts)

	# write out the shape file FIELDS data
	# coveragePoly.record(os.path.basename(fileName), recDate)
	#as per client requirement point to the files only from the 01_levant part.
	#coveragePoly.record(fileName, fileName, fileimg)
	
	# recDate = datetime.now().strftime("%Y%m%d")
	dataimg = "file://" + fileName
	dataimg = os.path.basename(fileName)
	if '01_levant' in fileName.lower():
		fileimg = fileName[fileName.lower().find('01_levant'):]
	else:
		fileimg = fileName.lower()
	fileimg = os.path.normpath(fileimg)
	# coveragePoly.record(dataimg, fileimg, idx)
	coveragePoly.record(dataimg.upper(), fileimg.upper())
	
	return coveragePoly

###############################################################################
def injectposition2xtf(filename, positionfilename, outfilename):
	'''inject a position.txt file into an xtf file using infinitytool'''

	odirlog = makedirs(os.path.dirname(outfilename))

	cmd = "c:/infinitytool/tool/toolinjectpositonsintoxtf.exe"
	cmd = cmd + " -i %s" % (filename) 					#input filename
	cmd = cmd + " -p %s" % (positionfilename)			#position file name from where we ge tth erevised positons.
	cmd = cmd + " -o %s" % (outfilename)				#output name
	args = shlex.split(cmd)			
		
	fileout = os.path.join(odirlog, "xtf2laz.stdout.txt")
	fout = open(fileout, 'w')
	fileerr = os.path.join(odirlog, "xtf2laz.stderr.txt")
	ferr = open(fileerr, 'w')

	proc = subprocess.Popen(args, stdout=fout, stderr=ferr)
	proc.wait()

	return outfilename
# ###############################################################################
# def xtf2laz(filename, outfilename, resolution, iclippercentage=50, oclippercentage=95, epsg='31984'):
# 	'''convert an xtf file to a laz file using infinitytool'''

# 	odirlog = makedirs(os.path.dirname(outfilename))

# 	cmd = "c:/infinitytool/tool/toolxtf2laz.exe"
# 	cmd = cmd + " -i %s" % (filename) 					#input filename
# 	cmd = cmd + " -o %s" % (outfilename)				#output name
# 	cmd = cmd + " -x"									# apply X lever arm
# 	cmd = cmd + " -y"									# apply Y lever arm
# 	cmd = cmd + " -p"									# pitch correction
# 	cmd = cmd + " -iclip %s " % (str(iclippercentage))	# inner clip percentage
# 	cmd = cmd + " -oclip %s " % (str(oclippercentage))	# inner clip percentage
# 	cmd = cmd + " -r %s" % (str(resolution)) 
# 	args = shlex.split(cmd)			
		
# 	fileout = os.path.join(odirlog, "xtf2laz.stdout.txt")
# 	fout = open(fileout, 'w')
# 	fileerr = os.path.join(odirlog, "xtf2laz.stderr.txt")
# 	ferr = open(fileerr, 'w')

# 	proc = subprocess.Popen(args, stdout=fout, stderr=ferr)
# 	proc.wait()

# 	return outfilename

###############################################################################
def xtf2laz(filename, outfilename, resolution, iclippercentage=50, oclippercentage=95, epsg='31984'):
	'''convert an xt file to a laz file using infinitytool'''

	odirlog = makedirs(os.path.dirname(outfilename))

	cmd = "c:/infinitytool/tool/toolxtf2laz.exe"
	cmd = cmd + " -i %s" % (filename) 					#input filename
	cmd = cmd + " -o %s" % (outfilename)				#output name
	cmd = cmd + " -x"									# apply X lever arm
	cmd = cmd + " -y"									# apply Y lever arm
	cmd = cmd + " -p"									# pitch correction
	cmd = cmd + " -iclip %s " % (str(iclippercentage))	# inner clip percentage
	cmd = cmd + " -oclip %s " % (str(oclippercentage))	# inner clip percentage
	cmd = cmd + " -r %s" % (str(resolution)) 
	args = shlex.split(cmd)			
		
	fileout = os.path.join(odirlog, "xtf2laz.stdout.txt")
	fout = open(fileout, 'w')
	fileerr = os.path.join(odirlog, "xtf2laz.stderr.txt")
	ferr = open(fileerr, 'w')

	proc = subprocess.Popen(args, stdout=fout, stderr=ferr)
	proc.wait()

	return outfilename

##############################################################################
##############################################################################
##############################################################################
def ficha(config, filename):
	'''fix up the CATHX laser raw data files to PB specification for 6-1-3-7_raw_laser_profiler'''

	'''if the words of the file are in the format code_auvpassmission_date rearrange into code_date_auvpassmission '''
	#  "D:\FortnightDeliverables_20190725_20190807\WideArea\20190730_6100_015_P11_S01_AML_014\2_post_dive_48h\6-1-2-13_coverage_all_sensors_DGN\MBES\CMBES_9VA00014_20190730.dgn"
	# IN  FICHA_SON_9VA00014_20190730.pdf
	# OUT FICHA_SON_20190730_9VA.pdf

	root, ext = os.path.splitext(filename)
	words = root.upper().split("_")
	result = filename
	if len(words) == 4:
		# only mangle files if they have the fly height
		if '9VA' in words[2]:
			# words[2] = "6100VA%s" % (words[2][4:])
			words[2] = "6100VA"
		elif '4VA' in words[2]:
			# words[2] = "4900VA%s" % (words[2][4:])
			words[2] = "4900VA"
		elif '6100VA' in words[2]:
			words[2] = words[2]
		elif '6100VA' in words[2]:
			words[2] = words[2]
		
		
		if 'VB' in words[2] or 'VA' in words[2]:
			result = "%s_%s_%s_%s%s" %(words[0], words[1], words[3], words[2], ext)
	return result

##############################################################################
def foto(config, filename):
	'''fix up the CATHX laser raw data files to PB specification for 6-1-3-7_raw_laser_profiler'''

	result 			= filename
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	# surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	area 			= decode(config, "{area}")
	imagedate 		= getdatefromcathxphotofilename(filename)

	result = "FT_%s_%s%s%s%s%s" % (imagedate, area, auv, passtypeshort, divenumber, ext)
#	result = "FT_%s_%s_%s%s%s%s" % (area, imagedate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def laser(config, filename):
	'''fix up the CATHX laser raw data files to PB specification for 6-1-3-7_raw_laser_profiler'''
	'''for zip files...'''
	'''from this...cathx-0006-ABDC054-0001-00.cathxlaser.zip'''
	'''to this...BLAS_0006_ABDC054_20190828_ABD4900VB0038.zip'''
	'''for laz files...'''

	result 			= filename
	ext 		= os.path.splitext(os.path.expanduser(filename.lower()))[1]

	if 'zip' in ext:
		auv 			= decode(config, "{auv}")
		surveydate 		= decode(config, "{surveydate}")
		passtypeshort 	= decode(config, "{passtypeshort}")
		divenumber 		= getdivenumber(config)
		area 			= decode(config, "{area}")
		sequence 		= getsequencefromlaserfilename(filename)
		subsequence 	= getsubsequencefromlaserfilename(filename)
		linename 		= getlinenamefromlaserfilename(filename)

		divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

		# root = root.replace('.cathxlaser', '')
		# root = root.replace('cathx', '')
		# root = root.replace('-0001-00', '')
		# root = root.replace('-', '_')
		result = "BLAS_%s_%s_%s_%s_%s" % (sequence, subsequence, linename, surveydate, divecode)

		#original
		# result = "BLAS_%s_%s_%s_%s" % (sequence, linename, surveydate, divecode)

		# result = "BLAS_%s_%s_%s_%s_000000_%s%s%s%s" % (sequence, area, linename, surveydate, auv, passtypeshort, divenumber, ext)

		# root = root.replace('__', '_')
	
	if 'laz' in ext:
		'''from this...CATHX_D2019-08-29T00-36-52-235Z_None_X0.laz'''
		'''to this... '''
		auv 			= decode(config, "{auv}")
		# surveydate 		= decode(config, "{surveydate}")
		#revised survey date source
		surveydate 		= getdatefromlazfilename(filename)
		passtypeshort 	= decode(config, "{passtypeshort}")
		divenumber 		= getdivenumber(config)
		area 			= decode(config, "{area}")
		sequence 		= getsequencefromlazfilename(filename)
		linename 		= "noline" #getlinenamefromlaserfilename(filename)
		divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

		# root = root.replace('.cathxlaser', '')
		# root = root.replace('cathx', '')
		# root = root.replace('-0001-00', '')
		# root = root.replace('-', '_')
		result = "BLAS_%s_%s_%s_%s" % (sequence, linename, surveydate, divecode)
	

	return result

##############################################################################
def cmbesshp(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CMBES_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def clasdgn(config, filename):
	'''fix up the LIDAR coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CLAS_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result
##############################################################################
def cftdgn(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CFT_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def csonbdgn(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CSONB_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result
##############################################################################
def cmbesdgn(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CMBES_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result
##############################################################################
def csondgn(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CSON_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result
##############################################################################
def posplotshp(config, filename):
	'''fix up the post processed filename ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	#we have a line shp file so deal with it...
	if 'line' in filename.lower():
		result = "NAVLN_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	#we have a line shp file so deal with it...
	if 'line' in filename.lower():
		result = "NAVLN_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)


	return result
##############################################################################
def posplot(config, filename):
	'''fix up the post processed filename ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "POSP_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result
##############################################################################
def preplot(config, filename):
	'''fix up the proposed survey shape file'''

	result 			= filename
	if 'proposed_survey' in filename.lower():
		ext 		= os.path.splitext(os.path.expanduser(filename))[1]
		auv 			= decode(config, "{auv}")
		surveydate 		= decode(config, "{surveydate}")
		passtypeshort 	= decode(config, "{passtypeshort}")
		divenumber 		= getdivenumber(config)

		result = "PREP_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def plotdgn(config, filename):
	'''fix up the proposed survey DGN file'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	surveydate 		= decode(config, "{surveydate}")
	vessel  		= decode(config, "{vessel}")
	result = "PLOT_%s%s%s" % (vessel, surveydate, ext)

	return result

##############################################################################
def offset(config, filename):
	'''add the suffix to the processed offset files'''

	result 			= filename
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]

	# deal with the arms.ini file...
	if 'arms.ini' in filename.lower():
		result = "OFFSET_INI_%s%s%s%s" % (auv, passtypeshort, divenumber, ext)

	# this is an AUV so handle it...
	if 'auv' in filename.lower():
		result = "OFFSET_AUV_%s%s%s%s" % (auv, passtypeshort, divenumber, ext)

	# for Island Pride
	if 'ip' in filename.lower():
		result = "OFFSET_IP_%s%s%s%s" % (auv, passtypeshort, divenumber, ext)


	return result

##############################################################################
def rawnavigation(config, filename):
	'''add the suffix to the raw navigation files'''

	result 			= filename
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	if 'topsideauvpos.txt' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "%s_%s%s%s%s" % (root, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def processednavigation(config, filename):
	'''add the suffix to the processed navigation files'''

	result 			= filename
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	if 'attitude_smooth.txt' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "%s_%s%s%s%s" % (root, auv, passtypeshort, divenumber, ext)

	if 'navlab_smooth.bin' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "%s_%s%s%s%s" % (root, auv, passtypeshort, divenumber, ext)

	if 'smoothedpressure.txt' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "%s_%s%s%s%s" % (root, auv, passtypeshort, divenumber, ext)

	if 'position_smooth.txt' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "%s_%s%s%s%s" % (root, auv, passtypeshort, divenumber, ext)

	if '.pdf' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "NAVLAB_RESULT_%s%s%s%s" % (auv, passtypeshort, divenumber, ext)

	if '.ppt' in filename.lower():
		root, ext = os.path.splitext(os.path.expanduser(filename))
		result = "NAVLAB_RESULT_%s%s%s%s" % (auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def processedtrack(config, filename):
	'''rename the processed trackplot shape file'''
	# IN:  'POSP_TrackLine_6100VA_20180708.dgn'
	# OUT: 'NAV_LN_20190708_6100VA0004.dgn'

	result = filename
	# area 			= decode(config, "{area}")
	surveydate 		= decode(config, "{surveydate}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]

	if 'trackline' in filename.lower():
		ext = os.path.splitext(os.path.expanduser(filename))[1]
		result = "NAV_LN_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	if 'trackpoint' in filename.lower():
		ext = os.path.splitext(os.path.expanduser(filename))[1]
		result = "NAV_PT_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def processedembackscatter(config, filename):
	'''rename the EM2040 Backscatter .TIF files as per PB specifications '''
	# 7.4. Backscatter processed data (preliminarily)
	# BACK - Area Acronym _ Type of line_ Line number _ Date _ Time _ AUV Number _ FlightHeight VA / VB or HF / LH _ Diving number.tif
	# OUT: BACK_B004_20190730_204850_AML4VB00015.tif
	# IN:  EM2040-0032-p6-20190709-215607_P.db_DepthTiff.tif

	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	line = getlinenamefromallfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromallfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "BACK_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result


	# result 			= filename
	# area 			= decode(config, "{area}")
	# type 			= gettypefromallfilename(filename)
	# #filedate 		= getdatefromallfilename(filename)
	# auv 			= decode(config, "{auv}")
	# passtypeshort 	= decode(config, "{passtypeshort}")
	# divenumber 		= getdivenumber(config)
	# sequence 		= getsequencefrombackscattertiffilename(filename)
	# line			= getlinenamefrombackscattertiffilename(filename)
	# root, ext 		= os.path.splitext(os.path.expanduser(filename))

	# try:
	# 	filedate 		= datetime.strftime(getdatefromallfilename(filename), "%Y%m%d_%H%M%S")
	# 	result 			= "BACK_%s_%s_%s%s%s%s%s" % (line, filedate, area, auv, passtypeshort, divenumber, ext)
	# except:
	# 	log("file cannot be renamed, passing through directly %s" %(filename))
	# return result

##############################################################################
def back(config, filename):
	'''rename the MBES backscatter files as per PB specifications 7.4 dados processados backscatter preliminary'''
	# IN: 'BACK_9VA_20190730_amlb001.tiff'
	# OUT BACK_AMLB004_20190730_204850_4VB00015.tif

	result 			= filename
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	line 			= getlinenamefrombackfilename(filename)
	ext = os.path.splitext(os.path.expanduser(filename))[1]

	try:
		filedate 		= datetime.strftime(getdatefrombackfilename(filename), "%Y%m%d_%H%M%S")
		result	= "BACK_%s_%s_%s%s%s%s" % (line, filedate, auv, passtypeshort, divenumber, ext)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		result	= "BACK_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def sonb(config, filename):
	'''rename the.TIF files as per PB specifications 8.4 dados processados de batimetria interferometrica preliminary'''
	# IN: 'sssb-20190801-024949-amlm025-3-S-SC-000.tif'
	# OUT: SONB_AMLB004_20190730_204850_P_4VA00016.tif


	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	channel 		= getchannelfromsasfilename(filename)
	line = getlinenamefromsasfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SONB_%s_%s_%s_%s" % (line, filedate, channel, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result


	# result 			= filename
	# auv 			= decode(config, "{auv}")
	# passtypeshort 	= decode(config, "{passtypeshort}")
	# divenumber 		= getdivenumber(config)
	# line 			= getlinenamefromsasfilename(filename)
	# channel 		= getchannelfromsasfilename(filename)
	# root, ext 		= os.path.splitext(os.path.expanduser(filename))

	# try:
	# 	filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
	# 	result	= "SONB_%s_%s_%s_%s%s%s%s" % (line, filedate, channel, auv, passtypeshort, divenumber, ext)
	# except:
	# 	log("file cannot be renamed, passing through directly %s" %(filename))
	# 	result	= "SONB_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	# return result

##############################################################################
def processedmbes(config, filename):
	'''rename the EM2040 .TIF files as per PB specifications 7.3 Dados processados'''
	# MBES - Area Acronym _ Type of line_ Line number _ Date _ Time _ AUV Number _ FlightHeight_ Diving number.tif
	# Exemplo: MBES_B004_20190730_204850_AML4VB00015.tif

	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	line = getlinenamefromallfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromallfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "MBES_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result

	# auv 			= decode(config, "{auv}")
	# passtypeshort 	= decode(config, "{passtypeshort}")
	# divenumber 		= getdivenumber(config)
	# root, ext 		= os.path.splitext(os.path.expanduser(filename))

	# words = root.split("-")
	# if len(words) > 4:
	# 	root = "%s_%s_%s_%s_%s" % (words[0], words[1], words[2], words[3], words[4][:6])


	# # if root[-2:] == '_P':
	# # 	root = root.replace('_P', '')
	
	# root 			= root.replace('EM2040', 'MBES')
	# result			= "%s_%s%s%s%s" % (root, auv, passtypeshort, divenumber, ext)

	# return result



##############################################################################
def surveyreport(config, filename):
	'''rename the survey reports to the PB standard'''
	# IN:  'SBP-0008-amlb004-20190730-204850_P.E.sgy'
	# OUT: 'AV-CC20880115.pdf'

	result 			= filename
	fortnight 		= decode(config, "{fortnight}")
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]

	result 			= "AV-%s%s" % (fortnight, ext)

	return result

##############################################################################
def controle(config, filename):
	'''rename the survey reports to the PB standard'''
	# IN:  anything
	# OUT: 'CONTROLE-A09-CC20880115.pdf'

	result 			= filename
	area 			= decode(config, "{area}")
	fortnight 		= decode(config, "{fortnight}")

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]

	result 			= "CONTROLE_%s_%s%s" % (area, fortnight, ext)

	return result

# ##############################################################################
# def sbpi(config, filename):
# 	'''rename the .sgy files as per PB specifications 12.1 Atributo analítico imaginário/ Imaginary analytical attribute'''
# 	# IN:  SBP-0008-amlb004-20190730-204850_P.I.sgy
# 	# OUT: SBPI_0008_AMLB004_20190730_204850_9VA00014.sgy

# 	segment 		= getsegmentfromsgyfilename(filename)
# 	auv 			= decode(config, "{auv}")
# 	passtypeshort 	= decode(config, "{passtypeshort}")
# 	divenumber 		= getdivenumber(config)
# 	line 			= getlinenamefromsgyfilename(filename)
# 	root, ext 		= os.path.splitext(os.path.expanduser(filename))

# 	try:
# 		filedate 		= datetime.strftime(getdatefromsgyfilename(filename), "%Y%m%d_%H%M%S")
# 		result 			= "SBPI_%s_%s_%s_%s%s%s%s" % (segment, line, filedate, auv, passtypeshort, divenumber, ext)
# 	except:
# 		log("file cannot be renamed, passing through directly %s" %(filename))
# 		result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

# 	return result


##############################################################################
def rawmbes(config, filename):
	'''rename the EM2040 .all files as per PB specifications 7.1 dados brutos'''
	# EM2040  _  Segment of the line  _  Area Acronym  _  Type of line  _  Line number  _  Date  _  Time  _  AUV Number  _  FlightHeight  _  Diving number.all
	# IN: EM2040_0043_AMLB004_20190730_204850.all
	# OUT: 'EM2040_b004_20190730_204850_AML6100VA00014.all'

	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	line = getlinenamefromallfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromallfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "EM2040_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	# if root[-2:] == '_P':
	# 	root = root.replace('_P', '')
	# root 			= root.replace("-", "_")
	# result 			= "%s_%s" % (root, divecode)
	
	return result



##############################################################################
def sbpi(config, filename):
	'''rename the .sgy files as per PB specifications 12.1 Atributo analítico imaginário/ Imaginary analytical attribute'''
	# IN:  SBP-0008-amlb004-20190730-204850_P.I.sgy
	# OUT: SBPI_0008_AMLB004_20190730_204850_9VA00014.sgy

	result 			= filename
	area 			= decode(config, "{area}")
	# segment 		= getsegmentfromsgyfilename(filename)
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	line 			= getlinenamefromsgyfilename(filename)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	try:
		filedate 		= datetime.strftime(getdatefromsgyfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SBPI_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def sbpr(config, filename):
	'''rename the .sgy files as per PB specifications 12.1 Atributo analítico imaginário/ Imaginary analytical attribute'''
	# IN:  SBP-0008-amlb004-20190730-204850_P.I.sgy
	# OUT: 'SBPR_0008_AMLB004_20190730_204850_9VA00014.sgy'

	result 			= filename
	area 			= decode(config, "{area}")
	# segment 		= getsegmentfromsgyfilename(filename)
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	line 			= getlinenamefromsgyfilename(filename)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	try:
		filedate 		= datetime.strftime(getdatefromsgyfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SBPR_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPR_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def sbpe(config, filename):
	'''rename the .sgy files as per PB specifications 12.1 Atributo analítico imaginário/ Imaginary analytical attribute'''
	# IN:  'SBP-0008-amlb004-20190730-204850_P.E.sgy'
	# OUT: 'SBPE_0008_AMLB004_20190730_204850_9VA00014.sgy'

	result 			= filename
	area 			= decode(config, "{area}")
	# segment 		= getsegmentfromsgyfilename(filename)
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	line 			= getlinenamefromsgyfilename(filename)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	try:
		filedate 		= datetime.strftime(getdatefromsgyfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SBPE_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPE_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def son(config, filename):
	'''rename the SAS sonar files as per PB specifications '''
	# FROM 'sasi-20190731-000037-amlm023-1-L-WN-000_xtf-CH12.tiff'
	# TO SON_AMLB004_20190730_294859_4VA00016tif

	# FROM: sasi-20190730-205933-amlb005-1-L-WN-000.xtf
	# TO: SON_AMLB005_20190730_205933_9VA00014.xtf


	result 			= filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	line = getlinenamefromsasfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SON_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result


	# # area 			= decode(config, "{area}")
	# # type 			= gettypefromallfilename(filename)
	# auv 			= decode(config, "{auv}")
	# passtypeshort 	= decode(config, "{passtypeshort}")
	# divenumber 		= getdivenumber(config)
	# line 			= getlinenamefromsasfilename(filename)
	# root, ext 		= os.path.splitext(os.path.expanduser(filename))

	# if root[-2:] == '_P':
	# 	root = root.replace('_P', '')

	# try:
	# 	filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
	# 	result	= "SON_%s_%s_%s%s%s%s" % (line, filedate, auv, passtypeshort, divenumber, ext)
	# except:
	# 	log("file cannot be renamed, passing through directly %s" %(filename))
	# 	result	= "SON_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	# return result

##############################################################################
def sssb(config, filename):
	'''rename the SAS sonar files as per PB specifications 8.2 dados brutos de batimetria interferometrica'''
	# FROM 'sssb-20190730-204849-amlb004-1-P-SC-000.all'
	# TO 'sssb_amlb004_20190730-204849_P_4V00016.all'


	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	channel 		= getchannelfromsasfilename(filename)
	line = getlinenamefromsasfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SSSB_%s_%s_%s_%s" % (line, filedate, channel, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result


	# # area 			= decode(config, "{area}")
	# result 			= filename
	# auv 			= decode(config, "{auv}")
	# passtypeshort 	= decode(config, "{passtypeshort}")
	# divenumber 		= getdivenumber(config)
	# line 			= getlinenamefromsasfilename(filename)
	# channel 		= getchannelfromsasfilename(filename)
	# root, ext 		= os.path.splitext(os.path.expanduser(filename))

	# try:
	# 	filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
	# 	result	= "SSSB_%s_%s_%s_%s%s%s%s" % (line, filedate, channel, auv, passtypeshort, divenumber, ext)
	# except:
	# 	log("file cannot be renamed, passing through directly %s" %(filename))
	# 	result	= "SSSB_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	# return result

##############################################################################
def bson(config, filename):
	'''rename the SAS sonar files as per PB specifications 8.1 dados brutos de sonar'''
	# FROM 'sasi-20190730-204849-amlb004-1-L-WN-000.xtf'
	# TO SASI_AMLB004_20190730_204850_4VA00016.xtf

	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)
	# channel 		= getchannelfromsasfilename(filename)
	line 			= getlinenamefromsasfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "SASI_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	return result

	# auv 			= decode(config, "{auv}")
	# passtypeshort 	= decode(config, "{passtypeshort}")
	# divenumber 		= getdivenumber(config)
	# line 			= getlinenamefromsasfilename(filename)
	# root, ext 		= os.path.splitext(os.path.expanduser(filename))

	# if root[-2:] == '_P':
	# 	root = root.replace('_P', '')

	# try:
	# 	filedate 		= datetime.strftime(getdatefromsasfilename(filename), "%Y%m%d_%H%M%S")
	# 	result	= "SASI_%s_%s_%s%s%s%s" % (line, filedate, auv, passtypeshort, divenumber, ext)
	# except:
	# 	log("file cannot be renamed, passing through directly %s" %(filename))
	# 	result	= "%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)

	# return result

##############################################################################
def rawsas(config, filename):
	'''rename the SASI .all files as per PB specifications for dados brutos'''
	# HISAS  _  Segment of the line  _  Type of line  _  Line number  _  Date  _  Time  _  Area Acronym _ AUV Number  _  FlightHeight  _  Diving number.all
	# IN:  hisas-0014-amlb004-20190730-204850.all
	# OUT: 'hisas_b004_20190730_204850_AML6100VA00014.all'

	result = filename
	area 			= decode(config, "{area}")
	auv 			= decode(config, "{auv}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)
	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	divecode 		= "%s%s%s%s%s" % (area, auv, passtypeshort, divenumber, ext)

	line = getlinenamefromallfilename(filename)

	try:
		filedate 		= datetime.strftime(getdatefromallfilename(filename), "%Y%m%d_%H%M%S")
		result 			= "HISAS_%s_%s_%s" % (line, filedate, divecode)
	except:
		log("file cannot be renamed, passing through directly %s" %(filename))
		# result	= "SBPI_%s_%s%s%s%s" % (filename, auv, passtypeshort, divenumber, ext)


	# if root[-2:] == '_P':
	# 	root = root.replace('_P', '')
	# root 			= root.replace("-", "_")
	# result 			= "%s_%s" % (root, divecode)
	
	return result

###############################################################################
def getsegmentfromallfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	segment = "SEGM"
	if len(words) > 1:
		segment = words[1]
	return segment

###############################################################################
def linenamesplit(s):
	head = s.rstrip('0123456789')
	tail = s[len(head):]
	return head, tail

###############################################################################
def getlinenamefromsasfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "NAME"
	if len(words) > 4:
		result = words[3].upper()

		linename, linenumber = linenamesplit(words[3].upper())
		if len(linenumber) > 0:
			result = "%s%03d" % (linename, int(linenumber))
		else:
			result = linename

	return result

###############################################################################
def getlinenamefromclientfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "NAME"
	if len(words) > 4:
		result = words[1].upper()
	return result

###############################################################################
def getlinenamefrombackfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("_")
	result = "NAME"
	if len(words) > 2:
		result = words[3].upper()
	return result

###############################################################################
def getchannelfromsasfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "CHANNEL"
	if len(words) > 6:
		result = words[5].upper()
	return result

###############################################################################
def gettypefromallfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	words = root.split("-")
	result = "TYPE"
	if len(words) > 2:
		result = words[2].upper()

		# if the word is not any of these, it is a mainline...
		# "M" mainline for main lines;
		# "B" Boxing for orthogonal lines to the equipment;
		# "X" Cross for the transverse lines;
		# "R" re-run for refazing lines;
		# "I" infill for gap coverage;
		# "RM" Main line re-run for the main refazing lines;
		# "RX" Cross line re-run for the transverse lines of refazing;

		# "START" to start the dive;
		# "DIVE" descent into the dive;
		# "RI" run in for entering the line;
		# "RO" run out for line output;
		# "LOITER" for mission termination or when the AUV is in standby mode.
		# Note: The data with the line types "START", "DIVE", "RI", "RO" and "LOITER" should be allocated in a folder with the suffix "_FORA_ESCOPO". (out scope)
		types = ["start", "dive", "ri", "ro", "loiter"]
		for type in types:
			if result.lower() == type.lower():
				result = type
				return result
		# we did not find a match, so it is a main line
		result = "M"

	return result

# ###############################################################################
def getdatefromcathxphotofilename(filename):
	'''extract the DATE from the cathx image file name'''
	# IN: image_D2019-07-09T09-35-33-317515Z_1.tiff
	#  OUT: 20190709-093533
	result = "DATE"
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("image_", "")
	root = root.replace("_", "-")
	words = root.split("-")
	if len(words) > 5:
		result = "%s%s%s_%s%s%s%s_%s" % (words[0][-4:], words[1], words[2][:2], words[2][-2:], words[3], words[4], words[5][:3], words[6])
		#result = "%s%s%s_%s%s%s_%s" % (words[0][-4:], words[1], words[2][:2], words[2][-2:], words[3], words[4], words[5][:1])
	return result

# ###############################################################################
# def getdatefromcathxphotofilename(filename):
# 	'''extract the DATE from the cathx image file name'''
# 	# IN: image_D2019-07-09T09-35-33-317515Z_1.tiff
# 	#  OUT: 20190709-093533
# 	result = "DATE"
# 	root = os.path.splitext(os.path.expanduser(filename))[0]
# 	root = root.replace("_", "-")
# 	words = root.split("-")
# 	if len(words) > 5:
# 		result = "%s%s%s_%s%s%s_%s" % (words[0][-4:], words[1], words[2][:2], words[2][-2:], words[3], words[4], words[5][:1])
# 	return result

###############################################################################
def getlinenamefromlaserfilename(filename):
	'''extract the LINENAME from the cathx laser zip file '''
	# IN: cathx-0005-RO-0001-00.cathxlaser.zip
	#  OUT: RO
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "SEQUENCE"
	if len(words) > 3:
		result = words[2]

		linename, linenumber = linenamesplit(words[2].upper())
		if len(linenumber) > 0:
			result = "%s%03d" % (linename, int(linenumber))
		else:
			result = linename
	return result

###############################################################################
def getsequencefrombackscattertiffilename(filename):
	'''extract the LINENAME from the file '''
	# IN: 'EM2040-0010-dv0-20190709-195049_P_BS.tiff'
	#  OUT: 0010
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "SEQ"
	if len(words) > 3:
		result = words[1]
	return result

###############################################################################
def getlinenamefrombackscattertiffilename(filename):
	'''extract the LINENAME from the file '''
	# IN: 'EM2040-0010-dv0-20190709-195049_P_BS.tiff'
	#  OUT: dv0
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "LINE"
	if len(words) > 3:
		result = words[2]
	return result

###############################################################################
def getsequencefromlaserfilename(filename):
	'''extract the sequence from the cathx laser zip file '''
	# IN: cathx-0005-RO-0001-00.cathxlaser.zip
	#  OUT: 0005
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace(".", "-")
	root = root.replace("_", "-")
	words = root.split("-")
	result = "SEQUENCE"
	if len(words) > 2:
		result = words[1]
	return result

###############################################################################
def getsubsequencefromlaserfilename(filename):
	'''extract the sub sequence (the second last word from the cathx laser zip file '''
	# IN: cathx-0005-RO-0001-00.cathxlaser.zip
	#  OUT: 0005
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace(".", "-")
	root = root.replace("_", "-")
	words = root.split("-")
	result = "SEQUENCE"
	if len(words) > 2:
		result = words[4]
	return result

# ###############################################################################
def getdatefromlazfilename(filename):
	'''extract the DATE from the cathx laser file name'''
	# IN: CATHX_D2019-08-29T00-36-52-235Z_None_X0.laz'''
	# OUT: 20190829
	result = "DATE"
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.lower().replace("cathx_", "")
	root = root.replace("_", "-")
	root = root.lower().replace("d", "")
	root = root.lower().replace("t", "-")
	words = root.split("-")
	if len(words) > 5:
		result = "%s%s%s" % (words[0], words[1], words[2])
		#result = "%s%s%s_%s%s%s_%s" % (words[0][-4:], words[1], words[2][:2], words[2][-2:], words[3], words[4], words[5][:1])
	return result

###############################################################################
def getsequencefromlazfilename(filename):
	'''extract the sequence from the cathx laser zip file '''
	# IN: CATHX_D2019-08-29T00-36-52-235Z_None_X0.laz'''
	# OUT: 00365223
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.lower().replace("cathx", "")
	root = root.lower().replace("t", "-")
	words = root.split("-")
	result = "SEQUENCE"
	if len(words) > 2:
		result = "%s%s%s%s" % (words[3], words[4], words[5], words[6][:3])
	return result

###############################################################################
def getdatefromsgyfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "DATE_TIME"
	if len(words) > 4:
		result = "%s_%s" % (words[3], words[4])
		result = datetime.strptime(words[3] + " " + words[4][:6], "%Y%m%d %H%M%S")
	return result

###############################################################################
def getdatefromallfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "DATE_TIME"
	if len(words) > 4:
		result = "%s_%s" % (words[3], words[4])
		result = datetime.strptime(words[3] + " " + words[4][:6], "%Y%m%d %H%M%S")
	return result

###############################################################################
def getdatefromsasfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	root = root.replace("_", "-")
	words = root.split("-")
	result = "DATE_TIME"
	if len(words) > 4:
		try:
			result = datetime.strptime(words[1] + " " + words[2], "%Y%m%d %H%M%S")
		except:
			result = root
	return result

###############################################################################
def getdatefrombackfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	words = root.split("_")
	result = "DATE_TIME"
	if len(words) > 2:
		result = datetime.strptime(words[2], "%Y%m%d")
	return result

###############################################################################
def splitall(path):
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts

###############################################################################
def getlinenamefromallfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	words = root.split("-")
	result = "NAME"
	if len(words) > 3:
		result = words[2].upper()

		linename, linenumber = linenamesplit(words[2].upper())
		if len(linenumber) > 0:
			result = "%s%03d" % (linename, int(linenumber))
		else:
			result = linename
	return result

##############################################################################
def replacecodesbyfunction(config, string, renamefunction):
	'''replace the codes in the file name with the values in the config files using a function to control each type.  this way we can undertake very complex logic / rules for each scenario'''
	
	result = string

	if renamefunction.lower() == 'ficha':
		result = ficha(config, string)
	elif renamefunction.lower() == 'classhp':
		result = classhp(config, string)
	elif renamefunction.lower() == 'cmbesshp':
		result = cmbesshp(config, string)
	elif renamefunction.lower() == 'csonbshp':
		result = csonbshp(config, string)
	elif renamefunction.lower() == 'csonshp':
		result = csonshp(config, string)
	elif renamefunction.lower() == 'clasdgn':
		result = clasdgn(config, string)
	elif renamefunction.lower() == 'cftshp':
		result = cftshp(config, string)
	elif renamefunction.lower() == 'cftdgn':
		result = cftdgn(config, string)
	elif renamefunction.lower() == 'csonbdgn':
		result = csonbdgn(config, string)
	elif renamefunction.lower() == 'cmbesdgn':
		result = cmbesdgn(config, string)
	elif renamefunction.lower() == 'csondgn':
		result = csondgn(config, string)
	elif renamefunction.lower() == 'posplotshp':
		result = posplotshp(config, string)
	elif renamefunction.lower() == 'posplot':
		result = posplot(config, string)
	elif renamefunction.lower() == 'foto':
		result = foto(config, string)
	elif renamefunction.lower() == 'laser':
		result = laser(config, string)
	elif renamefunction.lower() == 'preplot':
		result = preplot(config, string)
	elif renamefunction.lower() == 'plotdgn':
		result = plotdgn(config, string)
	elif renamefunction.lower() == 'offset':
		result = offset(config, string)
	elif renamefunction.lower() == 'rawnavigation':
		result = rawnavigation(config, string)
	elif renamefunction.lower() == 'processednavigation':
		result = processednavigation(config, string)
	elif renamefunction.lower() == 'processedtrack':
		result = processedtrack(config, string)
	elif renamefunction.lower() == 'processedbackscatter':
		result = processedembackscatter(config, string)
	elif renamefunction.lower() == 'processedmbes':
		result = processedmbes(config, string)
	elif renamefunction.lower() == 'rawsas':
		result = rawsas(config, string)
	elif renamefunction.lower() == 'rawmbes':
		result = rawmbes(config, string)
	elif renamefunction.lower() == 'sonb':
		result = sonb(config, string)
	elif renamefunction.lower() == 'bson':
		result = bson(config, string)
	elif renamefunction.lower() == 'sssb':
		result = sssb(config, string)
	elif renamefunction.lower() == 'son':
		result = son(config, string)
	elif renamefunction.lower() == 'back':
		result = back(config, string)
	elif renamefunction.lower() == 'cmbesshp':
		result = cmbesshp(config, string)
	elif renamefunction.lower() == 'sbpi':
		result = sbpi(config, string)
	elif renamefunction.lower() == 'sbpr':
		result = sbpr(config, string)
	elif renamefunction.lower() == 'sbpe':
		result = sbpe(config, string)
	elif renamefunction.lower() == 'surveyreport':
		result = surveyreport(config, string)
	elif renamefunction.lower() == 'controle':
		result = controle(config, string)
	elif renamefunction.lower() == 'passthru':
		result = string
	else:
		log("Oops, no transform function in place.  please check the codes in the mapping table.  we need to have a function with exactly the same name as the code in the mapping.csv")
	# if renamefunction == 'rawhisas':
	# 	result = rawem2040(config, string)
	return result

###############################################################################
def cftshp(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CFT_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def getdivenumber(config):
	dive = "DIVE"
	try:
		dive = "%04d" % (int(decode(config, "{divenumber}")))
		return dive
	except:
		return dive

##############################################################################
def csonbshp(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CSONB_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def csonshp(config, filename):
	'''fix up the sonar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CSON_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

##############################################################################
def classhp(config, filename):
	'''fix up the lidar coverage ignoring EVERYTHING the user specified'''

	ext 		= os.path.splitext(os.path.expanduser(filename))[1]
	auv 			= decode(config, "{auv}")
	surveydate 		= decode(config, "{surveydate}")
	passtypeshort 	= decode(config, "{passtypeshort}")
	divenumber 		= getdivenumber(config)

	result = "CLAS_%s_%s%s%s%s" % (surveydate, auv, passtypeshort, divenumber, ext)

	return result

###############################################################################
def getdatefromprocessedfilename(filename):
	root = os.path.splitext(os.path.expanduser(filename))[0]
	words = root.split("_")
	result = "DATE_TIME"
	if len(words) > 4:
		result = "%s_%s" % (words[2], words[3])
		# result = datetime.strptime(words[2] + " " + words[3], "%Y%m%d %H%M%S")
	return result

###############################################################################
def findFiles2(recursive, filespec, filter):
	'''tool to find files based on user request.  This can be a single file, a folder start point for recursive search or a wild card'''
	matches = []
	
	if not os.path.exists(filespec):
		return matches

	if recursive:
		matches = glob(os.path.join(filespec, "**", filter), recursive = True)
	else:
		matches = glob(os.path.join(filespec, filter))
	# if len(matches) == 0:
		# print ("Nothing found to copy to %s")
		# logging.info ("Nothing found to copy")
		# exit()
	# print (matches)
	#now fix up the bad crazy slashes so matching text works....
	cleanmatches = []
	for m in matches:
		cleanmatches.append(m.replace('\\','/'))
	return cleanmatches

###############################################################################
def copyfile(srcfile, dstfile, replace):
	'''Copy a file safely'''	
	
	# log ("Copying %s to %s" %(srcfile, dstfile))

	if not os.path.exists(srcfile):
		log ("source file does not exist, skipping : %s" % (srcfile))
		return 0

	if os.path.isfile(dstfile) and replace:
		# Handle errors while calling os.remove()
		try:
			os.remove(dstfile)
		except:			
			log("Error while deleting file %s. Maybe its in use?" % (dstfile), True)

		# Handle errors while calling os.ulink()
		# try:
		# 	os.unlink(dstfile)
		# except:
		# 	log("Error while deleting file %s. Maybe its in use?" % (dstfile), True)

	if os.path.exists(dstfile):
		log ("destination file exists, skipping : %s" % (dstfile), False, False)
		return 0

	# the file does not exist so copy it.
	try:
		shutil.copy(srcfile, dstfile)
		return 1
	except:
		log("Error while copying file 1st attempt%s" % (dstfile), True)
		#try 2 
		try:
			shutil.copy(srcfile, dstfile)
			return 1
		except:
			log("Error while copying file 2nd attempt%s" % (dstfile), True)
			#try 3
			try:
				shutil.copy(srcfile, dstfile)
				return 1
			except:
				log("***Error while copying 3rd attempt, skipping file %s" % (dstfile), True)
		return 0

###############################################################################
def	contractcode2sourcefolder(args, config, contractcode):
	'''find the 48 hour folder path name based on the contract code '''

	#load the mapping table.
	mappingtable 	= loadmappingtable(args)
	#now for each item in the contract, count the number of files...
	for mapper in mappingtable:
		folder 	= mapper[1].strip() #get the 48 FOUR FOLDER
		sheetcontractcode	= str(mapper[4]).strip()
		enabled				= str(mapper[6]).strip()
		renamefunction		= str(mapper[5]).strip()

		# the contract item in the mapping table is the same as the checklist, so count the files.
		# there could be more than 1 item in the mapping table
		if contractcode == sheetcontractcode:
			folder = replacecodes(config, folder)
			return (folder, enabled, renamefunction)
	
	log("ERROR: contract code %s is not found.  quitting.  Please ensure the mapping table has an entry for this contract code" % (contractcode))
	exit(0)

###############################################################################
def	contractcode2folder(args, config, contractcode):
	'''find the LEVANTAMENTO folder path based on the contract code '''

	# contractcode = ""
	#load the mapping table.
	mappingtable 	= loadmappingtable(args)
	#now for each item in the contract, count the number of files...
	for mapper in mappingtable:
		destinationfolder 	= mapper[2].strip() #get the LEVANTAMENTO FOLDER
		sheetcontractcode	= str(mapper[4]).strip()
		enabled				= str(mapper[6]).strip()
		# the contract item in the mapping table is the same as the checklist, so count the files.
		# there could be more than 1 item in the mapping table
		if contractcode == sheetcontractcode:
			destinationfolder = replacecodes(config, destinationfolder)
			return (destinationfolder, enabled)

	log("ERROR: contract code %s is not found.  quitting.  Please ensure the mapping table has an entry for this contract code" % (contractcode))
	exit(0)

###############################################################################
def	contractcode2levantamentofolder(args, config, contractcode):
	'''find the LEVANTAMENTO folder path based on the contract code '''

	#load the mapping table.
	mappingtable 	= loadmappingtable(args)
	#now for each item in the contract, count the number of files...
	for mapper in mappingtable:
		destinationfolder 	= mapper[2].strip() #get the LEVANATMENTO FOLDER
		sheetcontractcode	= str(mapper[4]).strip()
		renamefunction		= str(mapper[5]).strip()
		enabled				= str(mapper[6]).strip()
		# the contract item in the mapping table is the same as the checklist, so count the files.
		# there could be more than 1 item in the mapping table
		if contractcode == sheetcontractcode:
			destinationfolder = replacecodes(config, destinationfolder)
			destinationfolder = os.path.join(decode(config, "{outputfolder}"), destinationfolder)
			destinationfolder = os.path.normpath(destinationfolder)
			return (destinationfolder, enabled, renamefunction)

	log("ERROR: contract code %s is not found.  quitting.  Please ensure the mapping table has an entry for this contract code" % (contractcode))
	exit(0)

###############################################################################
def countfiles(srcfolder):
	if os.path.exists(srcfolder):
		#result = len([name for name in os.listdir(srcfolder) if os.path.isfile(name)])
		#onlyfiles = next(os.walk(srcfolder))[2]
		#result = len(onlyfiles)
		result = sum(len(files) for _, _, files in os.walk(srcfolder))
	else:
		result = "No"
	return result

###############################################################################
def	log(msg, error = False, printmsg=True):
		if printmsg:
			print (msg)
		if error == False:
			logging.info(msg)
		else:
			logging.error(msg)

###############################################################################
def	loadmappingtable(args):
	''' '''
	#look for mapping table in the mission folder first.  if not go to the python folder.
	filename = os.path.join(os.path.dirname(args.config), args.mappingtable)
	if not os.path.isfile(filename):
		print("!!oops, File: %s not found, trying delivery folder..." % (filename))
		filename = os.path.join(os.path.dirname(__file__), args.mappingtable)
	#if it is still not there, just accept the filename user proveded as the -mappingtable argument
	if not os.path.isfile(filename):
		print("!!oops, File: %s not found in delivery folder, trying script argument..." % (filename))
		filename = args.mappingtable

	if os.path.isfile(filename):
		with open(filename, 'r') as f:
			reader = csv.reader(f)
			mappingtable = list(reader)
			# remove the header
			mappingtable.pop(0)
			print("File: %s loaded, ok" % (filename))
			return mappingtable
	else:
		print("!!oops, File: %s not found, quitting!!" % (filename))
		exit()
	
###############################################################################
def loadconfig(args):
	'''read configuration variables from config file'''
	# filename = os.path.join(os.path.dirname(__file__), args.config)
	filename = args.config

	#if the user has specified just the mission, add the current working folder to the path
	if not os.path.isfile(filename):
		filename = os.path.join(os.path.dirname(__file__), filename)

	config = dict()

	if os.path.isfile(filename):
		config["root"] = os.path.dirname(filename)
		with open(filename, 'r') as f:
			for line in f:
				if len(line) < 5:
					continue
				if line.startswith("#"):
					continue
				words = line.split(",")
				if len(words) > 1:
					#add the root path to the various folders we need to work on
					if "{project}" in words[0]:
						config[words[0]] = os.path.abspath(os.path.join(config["root"], words[1].strip()))
					elif "{inputfolder}" in words[0]:
						config[words[0]] = os.path.abspath(os.path.join(config["root"], words[1].strip()))
					elif "{outputfolder}" in words[0]:
						config[words[0]] = os.path.abspath(os.path.join(config["root"], words[1].strip()))
					else:
						config[words[0]] = words[1].strip()

	else:
		log("!!oops, config table not found, quitting!!")
		exit()
	return config
	
##############################################################################
def replacecodes(config, string):

	#result = string
	'''replace the codes in the file name with the values in the config files '''
	for name, value in config.items():
		name = name.lower()
		string = string.lower()
		# if "sbp" in string:
		# 	print ("debug")
		if name in string:
		# 	string = string.replace(name, value)
		# if "(2)" in string:
		# 	print(string)
			string = string.replace(name, value)
			#string = re.sub(re.escape(name), value, string, flags=re.IGNORECASE)
	return string

##############################################################################
def creation_date(path_to_file):
	"""
	Try to get the date that a file was created, falling back to when it was
	last modified if that isn't possible.
	See http://stackoverflow.com/a/39501288/1709587 for explanation.
	"""
	import platform
	if platform.system() == 'Windows':
		return os.path.getctime(path_to_file)
	else:
		stat = os.stat(path_to_file)
		try:
			return stat.st_birthtime
		except AttributeError:
			# We're probably on Linux. No easy way to get creation dates here,
			# so we'll settle for when its content was last modified.
			return stat.st_mtime

##############################################################################
def decode(config, code):
	try:
		return config[code]
	except :
		log("oops code coude not be found, using code as value instead (this should not happen, please check the config file for this code:%s" %(code), True)
		return code

###############################################################################
def from_timestamp(unixtime):
	return datetime(1970, 1 ,1) + timedelta(seconds=unixtime)

###############################################################################
def to_timestamp(recordDate):
	return (recordDate - datetime(1970, 1, 1)).total_seconds()

###############################################################################
def update_progress(job_title, progress):
	length = 20 # modify this to change the length
	block = int(round(length*progress))
	msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
	if progress >= 1: msg += " DONE\r\n"
	sys.stdout.write(msg)
	sys.stdout.flush()

###############################################################################
def replacesuffix(filename, newsuffix):
	newfilename = os.path.splitext(filename)[0] + newsuffix
	return newfilename

###############################################################################
def deletefolder(foldername):
	'''Delete an empty directory'''
	#  Delete an empty directory using os.rmdir() and handle exceptions
	# try:
	# 	os.rmdir(foldername)
	# except:
	# 	print('Error while deleting directory')

	# Delete all contents of a directory using shutil.rmtree() and  handle exceptions
	# try:
	# 	shutil.rmtree(foldername)
	# except:
	# 	print('Error while deleting directory')

	# Delete all contents of a directory and ignore errors
	shutil.rmtree(foldername, ignore_errors=True)

	# Delete all contents of a directory and handle errors
	# shutil.rmtree(foldername, onerror=handleError )
 
###############################################################################
def handleError(func, path, exc_info):
	print('Handling Error for file ' , path)
	print(exc_info)
	# Check if file access issue
	if not os.access(path, os.W_OK):
		print('Hello')
		# Try to change the permision of file
		os.chmod(path, stat.S_IWUSR)
		# call the calling function again
		func(path)
###############################################################################
class RASTERSOURCE:
	'''load Petrobras tile conventions '''
	def __init__(self, filename):
		self.filename 			= filename
		self.clippedfilename 	= ""
		self.outlinefilename 	= ""
		self.outfilename 	= ""
		self.tilename 			= ""
		self.sortorder 			= ""
		self.rect 				= ""

	def decode(self, jdict):
		self.filename 			= jdict['filename']
		self.clippedfilename 	= jdict['clippedfilename']
		self.outlinefilename 	= jdict['outlinefilename']
		self.outfilename 		= jdict['outfilename']
		self.tilename 			= jdict['tilename']
		self.sortorder 			= jdict['sortorder']
		self.rect 				= jdict['rect']

###############################################################################
class TILE:
	def __init__(self, area, name, row, col, xorigin, yorigin, rectangle = None, path= ""):
		self.area = area
		self.name = name
		self.row 	= row
		self.col 	= col
		self.xorigin = xorigin
		self.yorigin = yorigin
		self.rectangle = rectangle
		self.path = path
		self.candidates = []
###############################################################################
class TILEFILE:
	'''load Petrobras tile conventions '''
	###########################################################################
	def __init__(self, filename):
		self.tiles = []

		# log("Loading tile file %s" % (filename))
		self.filename = filename
		filename = os.path.join(os.path.dirname(__file__), filename)
		if not os.path.exists( filename):
			log("oops, could not find tile files, quitting", True)
			return

		with open(filename, "r", ) as ins:
			for line in ins:
				words = line.split(",")

				if len(words) < 4:
					continue

				name = words[0]
				area = words[1]
				col = int(words[2])
				row = int(words[3])
				x = float(words[4])
				y = float(words[5])
				xorigin = x - 250
				yorigin = y - 250
				p1 = POINT(x - 250, y - 250)
				p2 = POINT(x + 250, y + 250)
				rectangle = RECT(p1, p2)
				self.tiles.append(TILE(area, name, row, col, xorigin, yorigin, rectangle))

###############################################################################
class TILEFILE2:
	'''load GG tile conventions '''
	###########################################################################
	def __init__(self, filename):
		self.tiles = []

		# log("Loading tile file %s" % (filename))
		self.filename = filename
		filename = os.path.join(os.path.dirname(__file__), filename)
		if not os.path.exists( filename):
			log("oops, could not find tile files, quitting", True)
			return

		with open(filename, "r", ) as ins:
			for line in ins:
				line = line.strip()
				words = line.split(",")

				if len(words) < 4:
					continue
				
				if "max" in line.lower():
					continue
				row = 0
				col = 0
				xmax = float(words[0])
				xmin = float(words[1])
				ymax = float(words[2])
				ymin = float(words[3])
				name = words[4]
				p1 = POINT(xmin, ymin)
				p2 = POINT(xmax, ymax)
				rectangle = RECT(p1, p2)
				self.tiles.append(TILE(area, name, row, col, xmin, ymin, rectangle))

###############################################################################
class TILEFILE:
	'''load Petrobras tile conventions '''
	###########################################################################
	def __init__(self, filename):
		self.tiles = []

		# log("Loading tile file %s" % (filename))
		self.filename = filename
		filename = os.path.join(os.path.dirname(__file__), filename)
		if not os.path.exists( filename):
			log("oops, could not find tile files, quitting", True)
			return

		with open(filename, "r", ) as ins:
			for line in ins:
				words = line.split(",")

				if len(words) < 4:
					continue

				name = words[0]
				area = words[1]
				col = int(words[2])
				row = int(words[3])
				x = float(words[4])
				y = float(words[5])
				xorigin = x - 250
				yorigin = y - 250
				p1 = POINT(x - 250, y - 250)
				p2 = POINT(x + 250, y + 250)
				rectangle = RECT(p1, p2)
				self.tiles.append(TILE(area, name, row, col, xorigin, yorigin, rectangle))

	###########################################################################
	def setarea(self, areatouse):
		subset = [] 
		for tile in self.tiles:
			if areatouse.lower() in tile.area.lower():
				subset.append(tile)
		self.tiles = subset

	###########################################################################
	def tileexists(self, tilename):
		for tile in self.tiles:
			if tilename.lower() in tile.name.lower():
				return True
		return False

###############################################################################
def overlap(r1, r2):
	'''Overlapping rectangles overlap both horizontally & vertically
	'''
	return range_overlap(r1.left, r1.right, r2.left, r2.right) and range_overlap(r1.bottom, r1.top, r2.bottom, r2.top)

###############################################################################
def range_overlap(a_min, a_max, b_min, b_max):
	'''Neither range is completely greater than the other
	'''
	overlapping = True
	if (a_min > b_max) or (a_max < b_min):
		overlapping = False
	return overlapping

###############################################################################
# def	savefilebboxrecord(args, config, filecontrol, filetype, prefix=""):

# 	#see if we have this already serialized.  if not then record it.
# 	fc = loadfilebboxrecord(args, config, filecontrol.filename, filetype, prefix)
# 	if fc is None:
# 		# Writing a JSON file
# 		config	= loadconfig(args)
# 		jsonfilename = os.path.join(config['root'], "log", "mission_" + config['{divenumber}'] + "_" + prefix + "_" + filetype + "_control.json")
# 		#jsonfilename = os.path.join(config['root'], "log", "mission_" + prefix + "_" + config['{divenumber}'] + "_" + filetype + "_control.json")
# 		#jsonfilename = os.path.join(config['root'], "log", "mission_" + config['{divenumber}'] + "_" + filetype + "_control.json")
# 		makedirs(os.path.dirname(jsonfilename.upper()))
# 		fp = open(jsonfilename, 'a')
# 		json_string = json.dumps(filecontrol.__dict__)
# 		fp.write(json_string + "\n")
# 		fp.close()

# ###############################################################################
# def	loadfilebboxrecord(args, config, filename, filetype, prefix=""):
# 	'''read a file control object from disc'''
# 	fc = None
# 	jsonfilename = os.path.join(config['root'], "log", "mission_" + config['{divenumber}'] + "_" + prefix + "_" + filetype + "_control.json")
# 	if os.path.exists(jsonfilename):
# 		fp = open(jsonfilename, 'r')
# 		for line in fp:
# 			record = json.loads(line)
# 			fc = RASTERSOURCE(filename)
# 			fc.decode(record)
# 			if fc.filename.lower() in filename.lower():
# 				fp.close()
# 				return fc
# 	return None

###############################################################################
###############################################################################
# class PHOTOINDEX(object):
# 	def __init__(self, filename):
# 		self.filename = filename
# 		self.conn = None
# 		self.cursor = None

# 		if os.path.exists(filename):
# 			self.conn = sqlite3.connect(filename)
# 		else:
# 			self.conn = sqlite3.connect(filename)
# 			c = conn.cursor()
# 			c.execute("CREATE TABLE IF NOT EXISTS photo(unix REAL, filename TEXT, left REAL, right REAL, bottom REAL, top REAL)")
# 			c.commit()
# 			c.close()

# 	###############################################################################
# 	def addphoto(self, filename, rect):
# 		c = self.conn.cursor()
# 		c.execute("INSERT INTO photo VALUES(filename, rect.left, rect.right, rect.bottom, rect.top)")
# 		c.commit()
# 		c.close()

# 	###############################################################################
# 	def findphoto(self, filename):
# 		c = self.conn.cursor()
# 		c.execute('SELECT * FROM photo WHERE filename = filename')
# 		data = c.fetchall()
# 		print(data)
# 		for row in data:
# 			print(row)
# 		c.close()

# 		return data
###############################################################################
###############################################################################
class POINT(object):
	def __init__(self, x, y, z=0.0):
		self.x = x
		self.y = y
		self.z = z

###############################################################################
class RECT(object):
	def __init__(self, p1, p2):
		'''Store the top, bottom, left and right values for points 
				p1 and p2 are the (corners) in either order
		'''
		self.left   = min(p1.x, p2.x)
		self.right  = max(p1.x, p2.x)
		self.bottom = min(p1.y, p2.y)
		self.top    = max(p1.y, p2.y)
		self.minz   = min(p1.z, p2.z)
		self.maxz   = max(p1.z, p2.z)

	def is_intersect(self, other):
			if self.left > other.right or self.right < other.left:
				return False
			if self.bottom > other.top or self.top < other.bottom:
				return False
			return True

	def __and__(self, other):
		if not self.is_intersect(other):
			return RECT([self.left, self.bottom], [self.right,self.top])
		left = min(self.left, other.left)
		right = max(self.right, other.right)
		bottom = min(self.bottom, other.bottom)
		top = max(self.top, other.top)
		return RECT([left, bottom], [right,top])

	intersect = __and__

	def __or__(self, other):
		left = min(self.left, other.left)
		right = max(self.right, other.right)
		bottom = min(self.bottom, other.bottom)
		top = max(self.top, other.top)
		return RECT([left, bottom], [right,top])

	union = __or__

	def __str__(self):
		return 'Rectangle({self.left},{self.right},{self.bottom},{self.top})'.format(self=self)

		@property
		def area(self):
			return (self.right - self.left) * (self.top - self.bottom)
###############################################################################
class POLYGON(object):
	def __init__(self, name, area, geometry):
		self.name= name
		self.area = area
		self.geometry = geometry
###############################################################################

###############################################################################
class CPROGRESS(object):
	'''thread safe class to display progress in command window when in multiprocess mode'''
	# procprogress = CPROGRESS(1000)
	# for i in range(1000):
	# 	time.sleep(0.01)
	# 	procprogress.increment_progress("test", i)

	###########################################################################
	def __init__(self, maxcount=100):
		self.length = 20 # modify this to change the length
		self.maxcount = max(maxcount,1)
		# self.progress = 0
		self.stime = datetime.now()
		self.value = 0
		self.msg = "Progress:"

	###########################################################################
	def setmaximum(self, value, current=0):
		self.maxcount = value
		self.value = current
		self.stime = datetime.now()

	###########################################################################
	def increment_progress(self, msg="", value=0):
		
		if len(str(msg)) > 0:
			self.msg = msg

		if value == 0:
			self.value = self.value + 1
		else:
			self.value = value

		progress = self.value/self.maxcount
		
		# print(value)
		secondsconsumed = (datetime.now() - self.stime).total_seconds()
		secondsperitem = secondsconsumed / max(self.value,1)
		secondsremaining = int((self.maxcount - self.value) * secondsperitem)
		timeremaining = str(timedelta(seconds=secondsremaining))
		block = int(round(self.length*progress))
		msg = "\r{0}: [{1}] {2:2.2f}% Remaining: {3}".format(self.msg, "#"*block + "-"*(self.length-block), round(progress*100, 2), timeremaining)
		if progress >= 1: msg += " DONE\r\n"
		sys.stdout.write(msg)
		sys.stdout.flush()

	###########################################################################
	def complete(self, msg):
		length = 20 # modify this to change the length
		progress = 1
		block = int(round(length*progress))
		msg = "\r{0}: [{1}] {2}%".format(msg, "#"*block + "-"*(length-block), round(progress*100, 2))
		if progress >= 1: msg += " DONE\r\n"
		sys.stdout.write(msg)
		sys.stdout.flush()

class MEMORYSTATUSEX(ctypes.Structure):
	_fields_ = [
		("dwLength", ctypes.c_ulong),
		("dwMemoryLoad", ctypes.c_ulong),
		("ullTotalPhys", ctypes.c_ulonglong),
		("ullAvailPhys", ctypes.c_ulonglong),
		("ullTotalPageFile", ctypes.c_ulonglong),
		("ullAvailPageFile", ctypes.c_ulonglong),
		("ullTotalVirtual", ctypes.c_ulonglong),
		("ullAvailVirtual", ctypes.c_ulonglong),
		("sullAvailExtendedVirtual", ctypes.c_ulonglong),
	]

	def __init__(self):
		# have to initialize this to the size of MEMORYSTATUSEX
		self.dwLength = ctypes.sizeof(self)
		super(MEMORYSTATUSEX, self).__init__()

###################################################################################################
if __name__ == "__main__":
	procprogress = CPROGRESS(100)
	for i in range(100):
		time.sleep(0.01)
		procprogress.increment_progress("test", i)

		# main()
