#name:			survey2mcap.py
#created:		June 2032
#by:			paul.kennedy@guardiangeomatics.com
#description:	python module to scan a folder structure, convert raw MBES files to Foxglove MCAP file format so we can easily QC and analyse data using Foxgloveextract the ellipsoidal heigths from the RAW files and create spatial datasets for creating a hydroid / geoid for a survey area
######################
# DONE
# basic script structure
# clean out SSDM
# simplify command line args
# process entire folder

######################
######################
# 2DO
# add core py module for creating mcap files
# add SPO position so fox can make trackplot
# add SKM motion data so we can make time series plots
# add MRZ as point cloud
# add MRZ ping header so we can QC
######################
######################

import sys
import time
import os
import tempfile
import ctypes
import fnmatch
import math
import json
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from datetime import datetime
from datetime import timedelta
import multiprocessing as mp
import multiprocesshelper

# local from the shared area...
# sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'shared'))
import fileutils
import kmall
import ggmbes
# from py7k import s7kreader
# from pygsf import GSFREADER

##############################################################################
def main():

	parser = ArgumentParser(description='Read any Survey folder and create a OGC compliant GEOPackage in the SSDM Schema summarising the survey. This is a distillation process extracting important spatial attributes from the survey in an automated and rigorous manner.',
			epilog='Example: \n To process all files under a root folder use -i c:/foldername \n', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', 		action='store', 		default="",		dest='inputfolder', 	help='input folder to process.')
	# parser.add_argument('-o', 		action='store', 		default="",		dest='outputFilename', 	help='output GEOPACKAGE filename.')
	parser.add_argument('-s', 		action='store', 		default="1",	dest='step', 			help='decimate the data to reduce the output size. [Default: 1]')
	# parser.add_argument('-odir', 	action='store', 		default="",	dest='odir', 			help='Specify a relative output folder e.g. -odir GIS')
	# parser.add_argument('-opath', 	action='store', 		default="",	dest='opath', 			help='Specify an output path e.g. -opath c:/temp')
	# parser.add_argument('-odix', 	action='store', 		default="",	dest='odix', 			help='Specify an output filename appendage e.g. -odix _coverage')
	# parser.add_argument('-epsg', 	action='store', 		default="4326",	dest='epsg', 			help='Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326')
	# parser.add_argument('-all', 	action='store_true', 	default=True, 	dest='all', 			help='extract all supported forms of data (ie do everything).')
	parser.add_argument('-reprocess', 	action='store_true', 	default=False, 	dest='reprocess', 			help='reprocess the survey folders by re-reading input files and creating new GIS features, ignoring the cache files. (ie do everything).')
	parser.add_argument('-cpu', 		dest='cpu', 			action='store', 		default='0', 	help='number of cpu processes to use in parallel. [Default: 0, all cpu]')

	args = parser.parse_args()
	# if len(sys.argv)==1:
	# 	parser.print_help()
	# 	sys.exit(1)

	if len(args.inputfolder) == 0:
		args.inputfolder = os.getcwd() + args.inputfolder

	if args.inputfolder == '.':
		args.inputfolder = os.getcwd() + args.inputfolder

	process(args)

###############################################################################
def process(args):
	if not os.path.isdir(args.inputfolder):
		print ("oops, input is not a folder.  Please specify a survey folder.")
		return

	surveyname = os.path.basename(args.inputfolder) #this folder should be the Survey NAME
	# if args.opath == "":
	# 	args.opath = os.path.join(args.inputfolder, "FOXGLOVE")

	# process any and all kmall files
	mp_processKMALL(args)

###############################################################################
def update_progress(job_title, progress):
	'''progress value should be a value between 0 and 1'''
	length = 20 # modify this to change the length
	block = int(round(length*progress))
	msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
	if progress >= 1: msg += " DONE\r\n"
	sys.stdout.write(msg)
	sys.stdout.flush()

################################################################################
def processKMALL(filename, outfilename, step):

	lastimestamp = 0
	step = float(step)

	#now read the kmall file and return the mcap filename if success

	# create the mcap file so we can start streaming results to it

	# setup the SPO headers
	# setup the MRZ headers
	# setup the XYZ headers
	# setup the SKM attitude headers
	
	# print("Loading KMALL Navigation...")
	r = kmall.kmallreader(filename)

	while r.moreData():
		try:
			# print(self.fileptr.tell())
			typeofdatagram, datagram = r.readDatagram()
			if (typeofdatagram == 'CORRUPT'):
				#we have seen corrupt kmall files when sis crashes.
				r.rewind()
				return outfilename
			if (r.recordTime - lastimestamp) < step:
				# skip...  performance increase
				continue

			if (typeofdatagram == '#SPO'):
				if (r.recordTime - lastimestamp) < step:
					# skip...  performance increase
					continue
				datagram.read()
				if positioniscorrupt(datagram):
					continue
				# navigation.append([to_timestamp(datagram.date), datagram.longitude, datagram.latitude, 0.0, datagram.heading])
				lastimestamp = r.recordTime

				#add to mcap file...

			if (typeofdatagram == '#MRZ'):
				datagram.read(True)
				if positioniscorrupt(datagram):
					continue

				# navigation.append([to_timestamp(datagram.date), datagram.latitude, datagram.longitude, datagram.ellipsoidHeightReRefPoint_m, datagram.heading, datagram.txTransducerDepth_m, datagram.z_waterLevelReRefPoint_m, 0.0])
				lastimestamp = r.recordTime

				#add ping header to mcap file...

				ph = ggmbes.GGPING()
				ph.timestamp 			= to_timestamp(datagram.date)
				ph.longitude 			= datagram.longitude
				ph.latitude 			= datagram.latitude
				ph.ellipsoidalheight 	= datagram.ellipsoidHeightReRefPoint_m
				ph.heading		 		= datagram.heading
				ph.pitch			 	= 0
				ph.roll			 		= 0
				ph.heave			 	= 0
				ph.tidecorrector	 	= 0 #datagram.txTransducerDepth_m # or is it this one? datagram.z_waterLevelReRefPoint_m
				ph.waterLevelReRefPoint_m = datagram.z_waterLevelReRefPoint_m
				ph.txTransducerDepth_m = datagram.txTransducerDepth_m
				# pingdata.append(ph)
				print (ph)
			if (typeofdatagram == '#SKM'):
				datagram.read()
				for sample in datagram.data:
					timestamp = (sample[3] + sample[4]/1000000000)
					#time, x, y, z, roll, pitch, heading, heave
					# attitude.append([timestamp, sample[6], sample[7], sample[8], sample[9], sample[10], sample[11], sample[12]])

		except:
			e = sys.exc_info()[0]
			print("Error: %s.  Please check file.  it seems to be corrupt: %s" % (e, r.fileName))
	r.rewind()
	return outfilename

################################################################################
def positioniscorrupt(datagram):
	# trap bad values
	if datagram.latitude < -90:
		return True
	
	if datagram.latitude > 90:
		return True
		
	if datagram.longitude < -180:
		return True
		
	if datagram.longitude > 180:
		return True

	return False	
################################################################################
def mp_processKMALL(args):
	''' decode the kmall files using multiple CPU and extract the ping and navigation data'''
	
	boundarytasks 		= []
	results 			= []
	outfilenames 		= []

	matches = fileutils.findFiles2(True, args.inputfolder, "*.kmall")

	for filename in matches:
		root = os.path.splitext(filename)[0]
		root = os.path.basename(filename)
		outfilename = os.path.join(filename + ".mcap").replace('\\','/')
		outfilenames.append(outfilename)
		if args.reprocess:
			if os.path.exists(outfilename):
				os.unlink(outfilename)

		if not os.path.exists(outfilename):
			boundarytasks.append([filename, outfilename])

	if args.cpu == '1':
		update_progress("Converting to MCAP", 0)
		for idx, filename in enumerate(matches):
			root = os.path.splitext(filename)[0]
			root = os.path.basename(filename)
			outfilename = os.path.join(filename + ".mcap").replace('\\','/')

			# the files exist so skip
			if os.path.exists(outfilename):
				continue
			result = processKMALL(filename, outfilename, args.step)
			results.append([filename, result])			
			update_progress("Converting to MCAP", (idx+1)/len(matches))
	else:
		multiprocesshelper.log("New kmall Files to Import: %d" %(len(boundarytasks)))		
		cpu = multiprocesshelper.getcpucount(args.cpu)
		multiprocesshelper.log("Extracting KMALL Navigation with %d CPU's" %(cpu))
		pool = mp.Pool(cpu)
		multiprocesshelper.g_procprogress.setmaximum(len(boundarytasks))
		poolresults = [pool.apply_async(processKMALL, (task[0], task[1], args.step), callback=multiprocesshelper.mpresult) for task in boundarytasks]
		pool.close()
		pool.join()
		for idx, result in enumerate (poolresults):
			results.append([boundarytasks[idx][0], result._value[0]])

	# now we can read the results files and create the geometry into the SSDM table
	multiprocesshelper.log("Files converted to MCAP: %d" %(len(results)))		

	multiprocesshelper.g_procprogress.setmaximum(len(results))

	return outfilenames

###############################################################################
def from_timestamp(unixtime):
	return datetime(1970, 1 ,1) + timedelta(seconds=unixtime)

###############################################################################
def to_timestamp(recordDate):
	return (recordDate - datetime(1970, 1, 1)).total_seconds()

###############################################################################
###############################################################################
###############################################################################
if __name__ == "__main__":
	main()
