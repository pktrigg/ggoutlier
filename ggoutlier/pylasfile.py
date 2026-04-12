#name:		  pylasfile
#created:	   July 2017
#by:			p.kennedy@guardiangeomatics.com
#description:   python module to read and write a ASPRS LAS file natively
#notes:		 See main at end of script for example how to use this
#based on ASPRS LAS 1.4-R13 15 July 2013

# See readme.md for more details

# LAS FORMAT DEFINITION:
# The format contains binary data consisting of
# A public header block,
# Any number of (optional) Variable Length Records (VLRs)
# The Point Data Records
# Any number of (optional) Extended Variable Length Records (EVLRs).
#
# All data are in little-endian format.
# The public header block contains generic data such as point numbers and point data bounds.

import os.path
import struct
import pprint
import time
import datetime
import math
import random

############################################################################
############################################################################
def main():

	testreader("C:/development/python/samplev1.2.las")
	# testreader("C:/development/python/samplev1.4.las")
	# testwriter()

############################################################################
def testwriter():
	'''
	sample write script so we can see how to use the code
	'''
	outFileName = os.path.join(os.path.dirname(os.path.abspath("c:/development/python/laswriter.las")), "laswriter.las")
	outFileName = createOutputFileName(outFileName)
	print("outputfile %s" % outFileName)
	# writer = laswriter(outFileName, 1.2)
	writer = laswriter(outFileName, 1.4)

	# write out a WGS variable length record so users know the coordinate reference system
	writer.writeVLR_WGS84()

	# now write some points
	writer.hdr.PointDataRecordFormat = 10
	for _ in range(100):
		# now add some random point data to the point lists
		writer.x.append(round(random.uniform(1, 100000),6))
		writer.y.append(round(random.uniform(1, 100000),6))
		writer.z.append(round(random.uniform(1, 1000),6))

	start_time = time.time() # time the process so we can keep it quick

	# before we write any points, we need to compute the bounding box, scale and offsets
	writer.computebbox_offsets()
	writer.writepoints()

	# we need to write the header after writing records so we can update the bounding box, point format etc 
	writer.writeHeader()
	writer.close()

	print("Write duration %.3fs" % (time.time() - start_time )) # time the process

	# testreader(outFileName)

###############################################################################
class laswriter:
	############################################################################
	def __init__(self, filename, lasformat=1.4):
		self.fileName = filename
		self.fileptr = open(filename, 'wb+')
		self.hdr = lashdr(lasformat)

		# the lists of all the data we will populate, then write into whatever format the user desires.  
		# these could be numpy arrays, but that introduces a dependency, so we will leave them as lists
		self.x 							= []
		self.y 							= []
		self.z 							= []
		self.intensity 					= []
		self.returnnumber 				= []
		self.numberreturns 				= []
		self.scandirectionflag 			= []
		self.edgeflightline 			= []
		self.classification 			= []
		self.scananglerank 				= []
		self.userdata 					= []
		self.pointsourceid 				= []
		self.gpstime 					= []
		self.red 						= []
		self.green 						= []
		self.blue 						= []
		self.wavepacketdescriptorindex 	= []
		self.byteoffsettowaveformdata 	= []
		self.waveformpacketsize 		= []
		self.returnpointwaveformlocation= []
		self.wavex 						= []
		self.wavey 						= []
		self.wavez 						= []
		self.nir 						= []

		self.classificationflags 		= []
		self.scannerchannel 			= []
		self.userdata 					= []
		self.scanangle 					= []

		self.supportedformats = self.hdr.getsuportedpointformats()

	############################################################################
	def writepointlist(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		self.computebbox_offsets()
		self.writepoints()

		# we need to write the header after writing records so we can update the bounding box, point format etc 
		self.writeHeader()

	############################################################################
	def writeVLR_WGS84(self):
		'''
		compose and write a standard variable length record for the WKTY of WGS84 CRS
		'''

		# before we write, we need to set the file pointer to the end of the VLR section , which is directly after the header block
		vlrl = self.getVLRTotalLength()
		self.fileptr.seek(self.hdr.HeaderSize + vlrl, 0)

		# now write out the vlr record
		vlrReserved				   = 0
		vlrUserid					 = b'LASF_Projection'
		vlrrecordid				   = 2112
		byte_str = 'WKT OGC COORDINATE SYSTEM'.encode('utf-8')
		byte_str = byte_str[:32].decode('utf-8', 'ignore').encode('utf-8')
		vlrDescription				= byte_str 
		vlrdata = b'GEOGCS["WGS 84",DATUM["World_Geodetic_System_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"],AXIS["Longitude",EAST],AXIS["Latitude",NORTH]]\x00'
		# vlrdata = b'PROJCS["WGS 84 / UTM zone 55S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AUTHORITY["EPSG","32755"]]\x00'
		vlrRecordLengthAfterHeader	= len(vlrdata)

		# now we have set the file pointer to the correct spot, write the record to disc
		record_struct = struct.Struct(self.hdr.vlrhdr14fmt)
		self.fileptr.write(record_struct.pack(vlrReserved, vlrUserid, vlrrecordid, vlrRecordLengthAfterHeader, vlrDescription))
		self.fileptr.write(vlrdata)

		self.hdr.NumberofVariableLengthRecords += 1
	############################################################################
	def writeVLR_WKT(self, wkt):
		'''
		compose and write a standard variable length record for the WKT from pyproj or rasterio
		'''

		# before we write, we need to set the file pointer to the end of the VLR section , which is directly after the header block
		vlrl = self.getVLRTotalLength()
		self.fileptr.seek(self.hdr.HeaderSize + vlrl, 0)

		# now write out the vlr record
		vlrReserved				   	= 0
		vlrUserid					= b'LASF_Projection'
		vlrrecordid				   	= 2112
		byte_str 					= 'WKT OGC COORDINATE SYSTEM'.encode('utf-8')
		byte_str 					= byte_str[:32].decode('utf-8', 'ignore').encode('utf-8')
		vlrDescription				= byte_str 
		vlrdata 					= wkt.encode('utf-8')
		vlrRecordLengthAfterHeader	= len(vlrdata)

		# now we have set the file pointer to the correct spot, write the record to disc
		record_struct = struct.Struct(self.hdr.vlrhdr14fmt)
		self.fileptr.write(record_struct.pack(vlrReserved, vlrUserid, vlrrecordid, vlrRecordLengthAfterHeader, vlrDescription))
		self.fileptr.write(vlrdata)

		self.hdr.NumberofVariableLengthRecords += 1
		
	############################################################################
	def getVLRTotalLength(self):
		'''get the variable length record size'''
		VLRTotalLength = 0
		# seek to the start of the vlrdata
		currentPosition = self.fileptr.tell()
		self.fileptr.seek(self.hdr.hdr14len, 0)
		for i in range (self.hdr.NumberofVariableLengthRecords):
			data = self.fileptr.read(self.hdr.vlrhdr14len)
			s = struct.unpack(self.hdr.vlrhdr14fmt, data)
			vlrRecordLengthAfterHeader = s[3]
			# now read the variable data
			vlrdata = self.fileptr.read(vlrRecordLengthAfterHeader)
			VLRTotalLength = VLRTotalLength + 54 + vlrRecordLengthAfterHeader
		self.fileptr.seek(currentPosition)
		return VLRTotalLength

	############################################################################
	def fit(self, s, l):
		u = s.encode("utf8")
		while True:
			if len(s) <= l:
				return s + "\0" * (l - len(s))
			u = u[:-1]
			s = u.encode("utf8")
		return s

	############################################################################
	def round_up(self, n, decimals=0):
		multiplier = 10 ** decimals
		return math.ceil(n * multiplier) / multiplier

	############################################################################
	def round_down(self, n, decimals=0):

		multiplier = 10 ** decimals
		return math.floor(n * multiplier) / multiplier
	
	############################################################################
	def computebbox_offsets(self):
		'''
		compute the bounding box of all records in the list
		'''
		rounding = 3 # was 9
		zrounding = 3
		self.hdr.MaxX = self.round_up(max(self.x), rounding)
		self.hdr.MinX = self.round_down(min(self.x), rounding)

		self.hdr.MaxY = self.round_up(max(self.y), rounding)
		self.hdr.MinY = self.round_down(min(self.y), rounding)

		self.hdr.MaxZ = self.round_up(max(self.z), rounding)
		self.hdr.MinZ = self.round_down(min(self.z), rounding)

		# self.hdr.MaxY = math.ceil(max(self.y))
		# self.hdr.MinY = math.floor(min(self.y))

		# self.hdr.MaxZ = math.ceil(max(self.z)) 
		# self.hdr.MinZ = math.floor(min(self.z))

		# xbuffer = (max(self.x) - min(self.x)) + (max(self.x) - min(self.x))/20
		# ybuffer = (max(self.y) - min(self.y)) + (max(self.y) - min(self.y))/20
		# zbuffer = (max(self.z) - min(self.z)) + (max(self.z) - min(self.z))/20

		# self.hdr.MaxX = max(self.x) + xbuffer
		# self.hdr.MinX = min(self.x) - xbuffer

		# self.hdr.MaxY = max(self.y) + ybuffer
		# self.hdr.MinY = min(self.y) - ybuffer

		# self.hdr.MaxZ = max(self.z) 
		# self.hdr.MinZ = min(self.z) 

		self.hdr.Xoffset = self.hdr.MinX 
		self.hdr.Yoffset = self.hdr.MinY
		self.hdr.Zoffset = self.hdr.MinZ

		digit2, afterDP2 = self.precision_and_scale(self.hdr.MaxX - self.hdr.MinX)
		self.hdr.Xscalefactor = 10**-(8-digit2)
		self.hdr.Xscalefactor = 10**-(rounding)

		digit2, afterDP2 = self.precision_and_scale(self.hdr.MaxY - self.hdr.MinY)
		self.hdr.Yscalefactor = 10**-(8-digit2)
		self.hdr.Yscalefactor = 10**-(rounding)

		digit2, afterDP2 = self.precision_and_scale(self.hdr.MaxZ - self.hdr.MinZ)
		self.hdr.Zscalefactor = 10**-(8-digit2)
		self.hdr.Zscalefactor = 10**-(zrounding)

		# self.hdr.Xscalefactor = 0.000000001
		# self.hdr.Yscalefactor = 0.000000001
		# self.hdr.Zscalefactor = 0.01

	def precision_and_scale(self, x):
		'''
		compute the number of digits beafore and after the decimal place so we can accurately scale a float into an integer
		'''
		max_digits = 14
		int_part = int(abs(x))
		magnitude = 1 if int_part == 0 else int(math.log10(int_part)) + 1
		if magnitude >= max_digits:
			return (magnitude, 0)
		frac_part = abs(x) - int_part
		if frac_part == 0.0:
			frac_part = 0.01
		multiplier = 10 ** (max_digits - magnitude)
		frac_digits = multiplier + int(multiplier * frac_part + 0.5)
		while frac_digits % 10 == 0:
			frac_digits /= 10
		scale = int(math.log10(frac_digits))
		return (magnitude+scale, scale)
		# return (scale) #return the number of digits after the decimal point only

	############################################################################
	def zerolistmaker(self, n):
		listofzeros = [0] * n
		return listofzeros

	############################################################################
	def onelistmaker(self, n):
		listofones = [1] * n
		return listofones

	############################################################################
	def fixemptylists(self):
		if len(self.intensity) == 0: 
			self.intensity = self.zerolistmaker(len(self.x))
		if len(self.returnnumber) == 0:
			self.returnnumber = self.onelistmaker(len(self.x))
		if len(self.numberreturns) == 0: 
			self.numberreturns = self.onelistmaker(len(self.x))
		if len(self.scandirectionflag) == 0: 
			self.scandirectionflag = self.zerolistmaker(len(self.x))
		if len(self.edgeflightline) == 0: 
			self.edgeflightline = self.zerolistmaker(len(self.x))
		if len(self.classification) == 0: 
			self.classification = self.zerolistmaker(len(self.x))
		if len(self.scananglerank) == 0: 
			self.scananglerank = self.zerolistmaker(len(self.x))
		if len(self.userdata) == 0: 
			self.userdata = self.zerolistmaker(len(self.x))
		if len(self.pointsourceid) == 0: 
			self.pointsourceid = self.zerolistmaker(len(self.x))
		if len(self.gpstime) == 0: 
			self.gpstime = self.zerolistmaker(len(self.x))
		if len(self.red) == 0: 
			self.red = self.zerolistmaker(len(self.x))
		if len(self.green) == 0: 
			self.green = self.zerolistmaker(len(self.x))
		if len(self.blue) == 0: 
			self.blue = self.zerolistmaker(len(self.x))
		if len(self.wavepacketdescriptorindex) == 0: 
			self.wavepacketdescriptorindex = self.zerolistmaker(len(self.x))
		if len(self.byteoffsettowaveformdata) == 0: 
			self.byteoffsettowaveformdata = self.zerolistmaker(len(self.x))
		if len(self.waveformpacketsize) == 0: 
			self.waveformpacketsize = self.zerolistmaker(len(self.x))
		if len(self.returnpointwaveformlocation) == 0: 
			self.returnpointwaveformlocation = self.zerolistmaker(len(self.x))
		if len(self.wavex) == 0: 
			self.wavex = self.zerolistmaker(len(self.x))
		if len(self.wavey) == 0: 
			self.wavey = self.zerolistmaker(len(self.x))
		if len(self.wavez) == 0: 
			self.wavez = self.zerolistmaker(len(self.x))
		if len(self.nir) == 0: 
			self.nir = self.zerolistmaker(len(self.x))

		if len(self.classificationflags) == 0:
			self.classificationflags = self.zerolistmaker(len(self.x))
		if len(self.scannerchannel) == 0:
			self.scannerchannel = self.zerolistmaker(len(self.x))
		if len(self.userdata) == 0:
			self.userdata = self.zerolistmaker(len(self.x))
		if len(self.scanangle) == 0:
			self.scanangle = self.zerolistmaker(len(self.x))

	############################################################################
	def writepoints(self):
		'''
		write points in the ASPRS version 1.4 format
		'''
		xs = self.hdr.Xscalefactor
		ys = self.hdr.Yscalefactor
		zs = self.hdr.Zscalefactor

		xo = self.hdr.Xoffset
		yo = self.hdr.Yoffset
		zo = self.hdr.Zoffset
	
		self.fixemptylists()

		self.hdr.LegacyNumberofpointrecords 	+= len(self.x)
		self.hdr.LegacyNumberofpointsbyreturn1 	+= len(self.x)
		self.hdr.Numberofpointrecords 			+= len(self.x)
		self.hdr.Numberofpointsbyreturn1 		+= len(self.x)
		if self.hdr.lasformat == 1.2:
			self.hdr.Offsettopointdata = self.hdr.hdr12len + self.getVLRTotalLength()
		if self.hdr.lasformat == 1.4:
			self.hdr.Offsettopointdata = self.hdr.hdr14len + self.getVLRTotalLength()

		if self.hdr.PointDataRecordFormat == 0:
			for i in range(len(self.x)):
				flags = self.setpointflags(self.returnnumber[i], self.numberreturns[i], self.scandirectionflag[i], self.edgeflightline[i])
				n = (int((self.x[i] - xo) / xs),
					int((self.y [i] - yo) / ys),
					int((self.z [i] - zo) / zs),
					int(self.intensity[i]),
					flags,
					self.classification[i],
					self.scanangle[i],
					self.userdata[i],
					self.pointsourceid[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return
		if self.hdr.PointDataRecordFormat == 1:
			format = self.supportedformats[self.hdr.PointDataRecordFormat][0]
			for i in range(len(self.x)):
				flags = self.setpointflags(self.returnnumber[i], self.numberreturns[i], self.scandirectionflag[i], self.edgeflightline[i])
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flags,
					self.classification[i],
					self.scanangle[i],
					self.userdata[i],
					self.pointsourceid[i],
					self.gpstime[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(format)
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 2:
			for i in range(len(self.x)):
				flags = self.setpointflags(self.returnnumber[i], self.numberreturns[i], self.scandirectionflag[i], self.edgeflightline[i])
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flags,
					self.classification[i],
					self.scanangle[i],
					self.userdata[i],
					self.pointsourceid[i],
					self.red[i],
					self.green[i],
					self.blue[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 3:
			for i in range(len(self.x)):
				flags = self.setpointflags(self.returnnumber[i], self.numberreturns[i], self.scandirectionflag[i], self.edgeflightline[i])
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flags,
					self.classification[i],
					self.scanangle[i],
					self.userdata[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.red[i],
					self.green[i],
					self.blue[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 4:
			for i in range(len(self.x)):
				flags = self.setpointflags(self.returnnumber[i], self.numberreturns[i], self.scandirectionflag[i], self.edgeflightline[i])
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flags,
					self.classification[i],
					self.scanangle[i],
					self.userdata[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.wavepacketdescriptorindex[i],
					self.byteoffsettowaveformdata[i],
					self.waveformpacketsize[i],
					self.returnpointwaveformlocation[i],
					self.wavex[i],
					self.wavey[i],
					self.wavez[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 5:
			for i in range(len(self.x)):
				flags = self.setpointflags(self.returnnumber[i], self.numberreturns[i], self.scandirectionflag[i], self.edgeflightline[i])
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flags,
					self.classification[i],
					self.scanangle[i],
					self.userdata[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.red[i],
					self.green[i],
					self.blue[i],
					self.wavepacketdescriptorindex[i],
					self.byteoffsettowaveformdata[i],
					self.waveformpacketsize[i],
					self.returnpointwaveformlocation[i],
					self.wavex[i],
					self.wavey[i],
					self.wavez[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 6:
			for i in range(len(self.x)):
				flag1 = self.setpointflag1_6_10(self.returnnumber[i], self.numberreturns[i])
				flag2 = self.setpointflag2_6_10(self.classificationflags[i], self.scannerchannel[i], self.scandirectionflag[i], self.edgeflightline[i] )
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flag1,
					flag2,
					self.classification[i],
					self.userdata[i],
					self.scanangle[i],
					self.pointsourceid[i],
					self.gpstime[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 7:
			for i in range(len(self.x)):
				flag1 = self.setpointflag1_6_10(self.returnnumber[i], self.numberreturns[i])
				flag2 = self.setpointflag2_6_10(self.classificationflags[i], self.scannerchannel[i], self.scandirectionflag[i], self.edgeflightline[i] )
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flag1,
					flag2,
					self.classification[i],
					self.userdata[i],
					self.scanangle[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.red[i],
					self.green[i],
					self.blue[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 8:
			for i in range(len(self.x)):
				flag1 = self.setpointflag1_6_10(self.returnnumber[i], self.numberreturns[i])
				flag2 = self.setpointflag2_6_10(self.classificationflags[i], self.scannerchannel[i], self.scandirectionflag[i], self.edgeflightline[i] )
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flag1,
					flag2,
					self.classification[i],
					self.userdata[i],
					self.scanangle[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.red[i],
					self.green[i],
					self.blue[i],
					self.nir[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 9:
			for i in range(len(self.x)):
				flag1 = self.setpointflag1_6_10(self.returnnumber[i], self.numberreturns[i])
				flag2 = self.setpointflag2_6_10(self.classificationflags[i], self.scannerchannel[i], self.scandirectionflag[i], self.edgeflightline[i] )
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flag1,
					flag2,
					self.classification[i],
					self.userdata[i],
					self.scanangle[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.wavepacketdescriptorindex[i],
					self.byteoffsettowaveformdata[i],
					self.waveformpacketsize[i],
					self.returnpointwaveformlocation[i],
					self.wavex[i],
					self.wavey[i],
					self.wavez[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

		if self.hdr.PointDataRecordFormat == 10:
			for i in range(len(self.x)):
				flag1 = self.setpointflag1_6_10(self.returnnumber[i], self.numberreturns[i])
				flag2 = self.setpointflag2_6_10(self.classificationflags[i], self.scannerchannel[i], self.scandirectionflag[i], self.edgeflightline[i] )
				n = (int((self.x[i] - xo) / xs),
					int((self.y[i] - yo) / ys),
					int((self.z[i] - zo) / zs),
					int(self.intensity[i]),
					flag1,
					flag2,
					self.classification[i],
					self.userdata[i],
					self.scanangle[i],
					self.pointsourceid[i],
					self.gpstime[i],
					self.red[i],
					self.green[i],
					self.blue[i],
					self.nir[i],
					self.wavepacketdescriptorindex[i],
					self.byteoffsettowaveformdata[i],
					self.waveformpacketsize[i],
					self.returnpointwaveformlocation[i],
					self.wavex[i],
					self.wavey[i],
					self.wavez[i]
					)
				# now write the record to disc
				record_struct = struct.Struct(self.supportedformats[self.hdr.PointDataRecordFormat][0])
				self.fileptr.write(record_struct.pack(*n))
			return

	############################################################################
	def close(self):
		self.fileptr.close()
		
	############################################################################
	def rewind(self):
		# go back to start of file
		self.fileptr.seek(0, 0)				

	############################################################################
	def seekPointRecordStart(self):
		# set the file pointer to the start of the points block
		self.fileptr.seek(self.hdr.Offsettopointdata, 0)				

	############################################################################
	def seekPointRecordEnd(self):
		# set the file pointer to the start of the points block
		self.fileptr.seek(self.hdr.Offsettopointdata + (self.hdr.Numberofpointrecords * self.hdr.PointDataRecordLength), 0)

	############################################################################
	def writeHeader(self):
		'''
		convert the header variables into a list, then conver the list into a tuple so we can pack it
		'''
		values = self.hdr.hdr2tuple()
		if self.hdr.lasformat == 1.2:
			s = struct.Struct(self.hdr.hdr12fmt)
		if self.hdr.lasformat == 1.4:
			s = struct.Struct(self.hdr.hdr14fmt)
		data = s.pack(*values)
		self.fileptr.seek(0, 0)				
		self.fileptr.write(data)

	############################################################################
	def setpointflags(self, returnnumber, numberreturns, scandirectionflag, edgeflightline ):
		flags = 0
		flags = self.setBitsFor_returnNo(flags, returnnumber)
		flags = self.setBitsFor_numberreturns(flags, numberreturns)
		flags = self.setBitsFor_scandirectionflag(flags, scandirectionflag)
		flags = self.setBitsFor_edgeflightline(flags, edgeflightline)
		return flags

	############################################################################
	def setpointflag1_6_10(self, returnnumber, numberreturns):
		flags = 0
		flags = self.setBitsFor_returnNo6_10(flags, returnnumber)
		flags = self.setBitsFor_numberreturns6_10(flags, numberreturns)
		return flags

	############################################################################
	def setpointflag2_6_10(self, classificationflags, scannerchannel, scandirectionflag, edgeflightline ):
		flags = 0
		flags = self.setBitsFor_classificationflags6_10(flags, classificationflags)
		flags = self.setBitsFor_scannerchannel6_10(flags, scannerchannel)
		flags = self.setBitsFor_scandirectionflag(flags, scandirectionflag)
		flags = self.setBitsFor_edgeflightline(flags, edgeflightline)
		return flags

	############################################################################
	def isBitSet(self, int_type, offset):
		'''testBit() returns a nonzero result, 2**offset, if the bit at 'offset' is one.'''
		mask = 1 << offset
		return (int_type & (1 << offset)) != 0

	############################################################################
	############################################################################
	def bitSet(self, v, offset):
		'''
		Set the index:th bit of v to 1 and return the new value.
		'''
		mask = 1 << offset   # Compute mask, an integer with just bit 'index' set.
		v |= mask		 
		return v

	############################################################################
	def setBitsFor_edgeflightline(self, int_type, edgeflightline):
		'''
		set the bit if this is the edge of a scan
		'''
		if edgeflightline: 
			int_type = self.bitSet(int_type, 7)
		return int_type

	############################################################################
	def setBitsFor_scandirectionflag(self, int_type, scandirectionflag):
		if scandirectionflag: #positive direction
			int_type = self.bitSet(int_type, 6)
		return int_type

	############################################################################
	def setBitsFor_numberreturns(self, int_type, numberreturns):
		'''
		packs the number of returns into the byte at the correct offset
		'''
		if numberreturns == 0:
			return int_type
		if numberreturns == 1:
			int_type = self.bitSet(int_type, 3)
			return int_type
		if numberreturns == 2:
			int_type = self.bitSet(int_type, 4)
			return int_type
		if numberreturns == 3:
			int_type = self.bitSet(int_type, 3)
			int_type = self.bitSet(int_type, 4)
			return int_type
		if numberreturns == 4:
			int_type = self.bitSet(int_type, 5)
			return int_type
		if numberreturns == 5:
			int_type = self.bitSet(int_type, 3)
			int_type = self.bitSet(int_type, 5)
			return int_type
		return int_type

	############################################################################
	def setBitsFor_returnNo(self, int_type, returnNo):
		'''
		packs the return number into the byte at the correct offset
		'''
		if returnNo == 0:
			return int_type
		if returnNo == 1:
			int_type = self.bitSet(int_type, 0)
			return int_type
		if returnNo == 2:
			int_type = self.bitSet(int_type, 1)
			return int_type
		if returnNo == 3:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			return int_type
		if returnNo == 4:
			int_type = self.bitSet(int_type, 2)
			return int_type
		if returnNo == 5:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 2)
			return int_type
		return int_type

	############################################################################
	def setBitsFor_returnNo6_10(self, int_type, returnNo):
		'''
		packs the return number into the byte at the correct offset for the las v1.4
		bits 0-3
		'''

		if returnNo == 0:
			return int_type
		if returnNo == 1:
			int_type = self.bitSet(int_type, 0)
			return int_type
		if returnNo == 2:
			int_type = self.bitSet(int_type, 1)
		if returnNo == 3:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			return int_type
		if returnNo == 4:
			int_type = self.bitSet(int_type, 2)
			return int_type
		if returnNo == 5:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 2)
			return int_type
		if returnNo == 6:
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			return int_type
		if returnNo == 7:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			return int_type
		if returnNo == 8:
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 9:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 10:
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 11:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 12:
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 13:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 14:
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if returnNo == 15:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		return int_type

	############################################################################
	def setBitsFor_numberreturns6_10(self, int_type, numberreturns):
		'''
		packs the number of return into the byte at the correct offset for the las v1.4
		# bits 4-7
		'''
		if numberreturns == 0:
			return int_type
		if numberreturns == 1:
			int_type = self.bitSet(int_type, 4)
			return int_type
		if numberreturns == 2:
			int_type = self.bitSet(int_type, 5)
			return int_type
		if numberreturns == 3:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 5)
			return int_type
		if numberreturns == 4:
			int_type = self.bitSet(int_type, 6)
			return int_type
		if numberreturns == 5:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 6)
			return int_type
		if numberreturns == 6:
			int_type = self.bitSet(int_type, 5)
			int_type = self.bitSet(int_type, 6)
			return int_type
		if numberreturns == 7:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 5)
			int_type = self.bitSet(int_type, 6)
			return int_type
		if numberreturns == 8:
			int_type = self.bitSet(int_type, 8)
			return int_type
		if numberreturns == 9:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 8)
			return int_type
		if numberreturns == 10:
			int_type = self.bitSet(int_type, 5)
			int_type = self.bitSet(int_type, 7)
			return int_type
		if numberreturns == 11:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 6)
			int_type = self.bitSet(int_type, 7)
			return int_type
		if numberreturns == 12:
			int_type = self.bitSet(int_type, 6)
			int_type = self.bitSet(int_type, 7)
			return int_type
		if numberreturns == 13:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 6)
			int_type = self.bitSet(int_type, 7)
			return int_type
		if numberreturns == 14:
			int_type = self.bitSet(int_type, 5)
			int_type = self.bitSet(int_type, 6)
			int_type = self.bitSet(int_type, 7)
			return int_type
		if numberreturns == 15:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 5)
			int_type = self.bitSet(int_type, 6)
			int_type = self.bitSet(int_type, 7)
			return int_type
		return int_type

	############################################################################
	def setBitsFor_classificationflags6_10(self, int_type, classificationflags):
		'''
		packs the classification flag at the correct offset for the las v1.4
		# bits 0-3
		'''
		if classificationflags == 0:
			return int_type
		if classificationflags == 1:
			int_type = self.bitSet(int_type, 0)
			return int_type
		if classificationflags == 2:
			int_type = self.bitSet(int_type, 1)
		if classificationflags == 3:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			return int_type
		if classificationflags == 4:
			int_type = self.bitSet(int_type, 2)
			return int_type
		if classificationflags == 5:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 2)
			return int_type
		if classificationflags == 6:
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			return int_type
		if classificationflags == 7:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			return int_type
		if classificationflags == 8:
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 9:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 10:
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 11:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 12:
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 13:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 14:
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		if classificationflags == 15:
			int_type = self.bitSet(int_type, 0)
			int_type = self.bitSet(int_type, 1)
			int_type = self.bitSet(int_type, 2)
			int_type = self.bitSet(int_type, 3)
			return int_type
		return int_type

	############################################################################
	def setBitsFor_scannerchannel6_10(self, int_type, scannerchannel):
		'''
		packs the scanner channel at the correct offset for the las v1.4
		# bits 4 & 5
		'''
		if scannerchannel == 0:
			return int_type
		if scannerchannel == 1:
			int_type = self.bitSet(int_type, 4)
			return int_type
		if scannerchannel == 2:
			int_type = self.bitSet(int_type, 5)
			return int_type
		if scannerchannel == 3:
			int_type = self.bitSet(int_type, 4)
			int_type = self.bitSet(int_type, 5)
			return int_type
		return int_type

###############################################################################
class lashdr:
	def __init__(self, lasformat=1.4):
		# version 1.2 header format
		self.hdr12fmt = "<4sHHLHH8sBB32s32sHHHLLBHL5L12d"
		self.hdr12len = struct.calcsize(self.hdr12fmt)

		# version 1.4 header format
		self.hdr14fmt = "<4sHHLHH8sBB32s32sHHHLLBHL5LddddddddddddQQLQ15Q"
		self.hdr14len = struct.calcsize(self.hdr14fmt)

		# variable length record same for v1.2 and v1.4
		self.vlrhdr14fmt = "<H16sHH32s"
		self.vlrhdr14len = struct.calcsize(self.vlrhdr14fmt)

		self.lasformat = lasformat # default to version 1.4.

		# create a default template for a V1.4 header.  We use this for writing purposes
		self.FileSignature 							= b'LASF'
		self.FileSourceID  							= 0
		self.GlobalEncoding 						= 17
		self.ProjectIDGUIDdata1 					= 0
		self.ProjectIDGUIDdata2 					= 0
		self.ProjectIDGUIDdata3 					= 0
		self.ProjectIDGUIDdata4 					= b"0"
		self.VersionMajor 							= 1
		if self.lasformat 							== 1.2:
			self.VersionMinor 						= 2
		else:
			self.VersionMinor 						= 4			

		self.SystemIdentifier 						= b'pylasfile'
		self.GeneratingSoftware 					= b'pylasfile'
		self.FileCreationDayofYear 					= datetime.datetime.now().timetuple().tm_yday
		self.FileCreationYear 						= datetime.datetime.now().year
		self.HeaderSize 							= 375
		self.Offsettopointdata 						= 0
		self.NumberofVariableLengthRecords 			= 0
		self.PointDataRecordFormat 					= 1
		self.PointDataRecordLength 					= 28

		self.LegacyNumberofpointrecords 			= 0
		self.LegacyNumberofpointsbyreturn1 			= 0
		self.LegacyNumberofpointsbyreturn2 			= 0
		self.LegacyNumberofpointsbyreturn3 			= 0
		self.LegacyNumberofpointsbyreturn4 			= 0
		self.LegacyNumberofpointsbyreturn5 			= 0
		self.Xscalefactor 							= 1
		self.Yscalefactor 							= 1
		self.Zscalefactor 							= 1

		self.Xoffset								= 0
		self.Yoffset 								= 0
		self.Zoffset 								= 0
		self.MaxX 									= 0
		self.MinX 									= 0
		self.MaxY 									= 0
		self.MinY 									= 0
		self.MaxZ 									= 0
		self.MinZ 									= 0

		self.StartofWaveformDataPacketRecord 		= 0
		self.StartoffirstExtendedVariableLengthRecord =	0
		self.NumberofExtendedVariableLengthRecords	= 0
		self.Numberofpointrecords 					= 0
		self.Numberofpointsbyreturn1 				= 0  
		self.Numberofpointsbyreturn2 				= 0  
		self.Numberofpointsbyreturn3 				= 0  
		self.Numberofpointsbyreturn4 				= 0  
		self.Numberofpointsbyreturn5 				= 0  

		self.Numberofpointsbyreturn6 				= 0  
		self.Numberofpointsbyreturn7 				= 0  
		self.Numberofpointsbyreturn8 				= 0  
		self.Numberofpointsbyreturn9 				= 0  
		self.Numberofpointsbyreturn10 				= 0  
		self.Numberofpointsbyreturn11 				= 0  
		self.Numberofpointsbyreturn12 				= 0  
		self.Numberofpointsbyreturn13 				= 0  
		self.Numberofpointsbyreturn14 				= 0  
		self.Numberofpointsbyreturn15 				= 0  

	############################################################################
	def __str__(self):
		'''
		pretty print this class
		'''
		return pprint.pformat(vars(self))

	############################################################################
	def getsuportedpointformats(self):
		'''
		returns a list of supported point file formats.
		'''
		s = []
		# format 0, v1.2,v1.4
		fmt = "<lllHBBbBH"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 1, v1.2,v1.4
		fmt = "<lllH BB b BH d"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 2, v1.2,v1.4
		fmt = "<lllH B BbBH HHH"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 3, v1.2,v1.4
		fmt = "<lllH B BbBH d HHH"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 4, v1.4
		fmt = "<lllH BBbBH d BQLffff"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])
	
		# format 5, v1.4
		fmt = "<lllH BBbB HdHH H BQLffff"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 6, v1.4
		fmt = "<lllH BBBB hHd"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])
	
		# format 7, v1.4
		fmt = "<lllHBBBBhHdHHH"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 8, v1.4
		fmt = "<lllHBBBBhHdHHHH"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 9, v1.4
		fmt = "<lllH BBBB hH d BQLffff"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		# format 10, v1.4
		fmt = "<lllH BB BB hHd HHHH BQLffff"
		fmtlen = struct.calcsize(fmt)
		s.append([fmt,fmtlen])

		return s

	############################################################################
	def hdr2tuple(self):
		'''
		convert the header properties into a tuple so we can easily write it to disc using struct
		'''	
		if self.lasformat == 1.2:
				return (
				self.FileSignature,
				self.FileSourceID,
				self.GlobalEncoding,
				self.ProjectIDGUIDdata1,
				self.ProjectIDGUIDdata2, 
				self.ProjectIDGUIDdata3, 
				self.ProjectIDGUIDdata4,
				self.VersionMajor,
				self.VersionMinor,
				self.SystemIdentifier,
				self.GeneratingSoftware,
				self.FileCreationDayofYear,
				self.FileCreationYear,
				self.HeaderSize,
				self.Offsettopointdata,
				self.NumberofVariableLengthRecords,
				self.PointDataRecordFormat,
				self.PointDataRecordLength,
				self.LegacyNumberofpointrecords,
				self.LegacyNumberofpointsbyreturn1,
				self.LegacyNumberofpointsbyreturn2,
				self.LegacyNumberofpointsbyreturn3,
				self.LegacyNumberofpointsbyreturn4,
				self.LegacyNumberofpointsbyreturn5,
				self.Xscalefactor,
				self.Yscalefactor,
				self.Zscalefactor,
				self.Xoffset,
				self.Yoffset,
				self.Zoffset,
				self.MaxX,
				self.MinX,
				self.MaxY,
				self.MinY,
				self.MaxZ,
				self.MinZ,
				)
		if self.lasformat == 1.4:
				return (
				self.FileSignature, 
				self.FileSourceID,
				self.GlobalEncoding,
				self.ProjectIDGUIDdata1,
				self.ProjectIDGUIDdata2, 
				self.ProjectIDGUIDdata3, 
				self.ProjectIDGUIDdata4,
				self.VersionMajor,
				self.VersionMinor,
				self.SystemIdentifier,
				self.GeneratingSoftware,
				self.FileCreationDayofYear,
				self.FileCreationYear,
				self.HeaderSize,
				self.Offsettopointdata,
				self.NumberofVariableLengthRecords,
				self.PointDataRecordFormat,
				self.PointDataRecordLength,
				self.LegacyNumberofpointrecords,
				self.LegacyNumberofpointsbyreturn1,
				self.LegacyNumberofpointsbyreturn2,
				self.LegacyNumberofpointsbyreturn3,
				self.LegacyNumberofpointsbyreturn4,
				self.LegacyNumberofpointsbyreturn5,
				self.Xscalefactor,
				self.Yscalefactor,
				self.Zscalefactor,

				self.Xoffset,
				self.Yoffset,
				self.Zoffset,
				self.MaxX,
				self.MinX,
				self.MaxY,
				self.MinY,
				self.MaxZ,
				self.MinZ,

				self.StartofWaveformDataPacketRecord,
				self.StartoffirstExtendedVariableLengthRecord,
				self.NumberofExtendedVariableLengthRecords,
				self.Numberofpointrecords,
				self.Numberofpointsbyreturn1,  
				self.Numberofpointsbyreturn2,  
				self.Numberofpointsbyreturn3,  
				self.Numberofpointsbyreturn4,  
				self.Numberofpointsbyreturn5,  

				self.Numberofpointsbyreturn6,  
				self.Numberofpointsbyreturn7,  
				self.Numberofpointsbyreturn8,  
				self.Numberofpointsbyreturn9,  
				self.Numberofpointsbyreturn10,  
				self.Numberofpointsbyreturn11,  
				self.Numberofpointsbyreturn12,  
				self.Numberofpointsbyreturn13,  
				self.Numberofpointsbyreturn14,  
				self.Numberofpointsbyreturn15,  
				)

	############################################################################
	def decodehdr(self, data):
		'''
		decode a header from a bytearray
		'''
		if self.lasformat == 1.2:
		
			s = struct.unpack(self.hdr12fmt, data)

			self.FileSignature 						= s[0]
			self.FileSourceID 						= s[1]
			self.GlobalEncoding 					= s[2]
			self.ProjectIDGUIDdata1 				= s[3]
			self.ProjectIDGUIDdata2 				= s[4]
			self.ProjectIDGUIDdata3 				= s[5]
			self.ProjectIDGUIDdata4 				= s[6]
			self.VersionMajor 						= s[7]
			self.VersionMinor 						= s[8]

			self.SystemIdentifier 					= s[9]
			self.GeneratingSoftware 				= s[10]
			self.FileCreationDayofYear 				= s[11]
			self.FileCreationYear 					= s[12]
			self.HeaderSize 						= s[13]
			self.Offsettopointdata 					= s[14]
			self.NumberofVariableLengthRecords 		= s[15]
			self.PointDataRecordFormat 				= s[16]
			self.PointDataRecordLength 				= s[17]

			self.LegacyNumberofpointrecords 		= s[18]
			self.LegacyNumberofpointsbyreturn1 		= s[19]
			self.LegacyNumberofpointsbyreturn2 		= s[20]
			self.LegacyNumberofpointsbyreturn3 		= s[21]
			self.LegacyNumberofpointsbyreturn4 		= s[22]
			self.LegacyNumberofpointsbyreturn5 		= s[23]
			self.Xscalefactor 						= s[24]
			self.Yscalefactor 						= s[25]
			self.Zscalefactor 						= s[26]
			self.Xoffset 							= s[27]

			self.Yoffset 							= s[28]
			self.Zoffset 							= s[29]
			self.MaxX 								= s[30]
			self.MinX 								= s[31]
			self.MaxY 								= s[32]
			self.MinY 								= s[33]
			self.MaxZ 								= s[34]
			self.MinZ 								= s[35]

		if self.lasformat == 1.4:
		
			s = struct.unpack(self.hdr14fmt, data)

			self.FileSignature 						= s[0]
			self.FileSourceID 						= s[1]
			self.GlobalEncoding 					= s[2]
			self.ProjectIDGUIDdata1 				= s[3]
			self.ProjectIDGUIDdata2 				= s[4]
			self.ProjectIDGUIDdata3 				= s[5]
			self.ProjectIDGUIDdata4 				= s[6]
			self.VersionMajor 						= s[7]
			self.VersionMinor 						= s[8]

			self.SystemIdentifier 					= s[9]
			self.GeneratingSoftware 				= s[10]
			self.FileCreationDayofYear 				= s[11]
			self.FileCreationYear 					= s[12]
			self.HeaderSize 						= s[13]
			self.Offsettopointdata 					= s[14]
			self.NumberofVariableLengthRecords 		= s[15]
			self.PointDataRecordFormat 				= s[16]
			self.PointDataRecordLength 				= s[17]

			self.LegacyNumberofpointrecords 		= s[18]
			self.LegacyNumberofpointsbyreturn1 		= s[19]
			self.LegacyNumberofpointsbyreturn2 		= s[20]
			self.LegacyNumberofpointsbyreturn3 		= s[21]
			self.LegacyNumberofpointsbyreturn4 		= s[22]
			self.LegacyNumberofpointsbyreturn5 		= s[23]
			self.Xscalefactor 						= s[24]
			self.Yscalefactor 						= s[25]
			self.Zscalefactor 						= s[26]
			self.Xoffset 							= s[27]

			self.Yoffset 							= s[28]
			self.Zoffset 							= s[29]
			self.MaxX 								= s[30]
			self.MinX 								= s[31]
			self.MaxY 								= s[32]
			self.MinY 								= s[33]
			self.MaxZ 								= s[34]
			self.MinZ 								= s[35]

			self.StartofWaveformDataPacketRecord 	= s[36]
			self.StartoffirstExtendedVariableLengthRecord =	s[37]
			self.NumberofExtendedVariableLengthRecords = s[38]
			self.Numberofpointrecords 				= s[39]
			self.Numberofpointsbyreturn1 			= s[40]
			self.Numberofpointsbyreturn2 			= s[41]
			self.Numberofpointsbyreturn3 			= s[42]
			self.Numberofpointsbyreturn4 			= s[43]
			self.Numberofpointsbyreturn5 			= s[44]

			self.Numberofpointsbyreturn6 			= s[45]
			self.Numberofpointsbyreturn7 			= s[46]
			self.Numberofpointsbyreturn8 			= s[47]
			self.Numberofpointsbyreturn9 			= s[48]
			self.Numberofpointsbyreturn10 			= s[49]
			self.Numberofpointsbyreturn11 			= s[50]
			self.Numberofpointsbyreturn12 			= s[51]
			self.Numberofpointsbyreturn13 			= s[52]
			self.Numberofpointsbyreturn14 			= s[53]
			self.Numberofpointsbyreturn15 			= s[54]

	############################################################################
	def get_PointDataRecordFormat(self):
		return self._PointDataRecordFormat

	############################################################################
	def set_PointDataRecordFormat(self, value):
		self._PointDataRecordFormat = value
		formats = self.getsuportedpointformats()
		self.PointDataRecordLength = formats[value][1]

	PointDataRecordFormat = property(get_PointDataRecordFormat,set_PointDataRecordFormat)

###############################################################################
class lasreader:
	def __init__(self, filename):
		if not os.path.isfile(filename):
			print ("file not found:", filename)
		self.fileName = filename
		self.fileptr = open(filename, 'rb')		
		self.fileSize = os.path.getsize(filename)
		self.hdr = lashdr()
		self.supportedformats = self.hdr.getsuportedpointformats()

		# the lists of all the data we will populate, then write into whatever format the user desires.  
		# these could be numpy arrays, but that introduces a dependency, so we will leave them as lists
		self.x 							= []
		self.y 							= []
		self.z 							= []
		self.intensity 					= []
		self.returnnumber 				= []
		self.numberreturns 				= []
		self.scandirectionflag 			= []
		self.edgeflightline 			= []
		self.classification 			= []
		self.scananglerank 				= []
		self.userdata 					= []
		self.pointsourceid 				= []
		self.gpstime 					= []
		self.red 						= []
		self.green 						= []
		self.blue 						= []
		self.wavepacketdescriptorindex 	= []
		self.byteoffsettowaveformdata 	= []
		self.waveformpacketsize			= []
		self.returnpointwaveformlocation= []
		self.wavex 						= []
		self.wavey 						= []
		self.wavez 						= []
		self.nir 						= []

		self.classificationflags 		= []
		self.scannerchannel 			= []
		self.userdata 					= []
		self.scanangle 					= []

	############################################################################
	def close(self):
		'''
		close the file
		'''
		self.fileptr.close()
		
	############################################################################
	def rewind(self):
		'''
		go back to start of file
		'''
		self.fileptr.seek(0, 0)				

	############################################################################
	def seekPointRecordStart(self):
		'''
		set the file pointer to the START of the points block so we can write some records
		'''
		self.fileptr.seek(self.hdr.Offsettopointdata, 0)				

	############################################################################
	def seekPointRecordEnd(self):
		'''
		set the file pointer to the END of the points block so we can add new records
		'''
		self.fileptr.seek(self.hdr.Offsettopointdata + (self.hdr.Numberofpointrecords*self.hdr.PointDataRecordLength), 0)

	############################################################################
	def __str__(self):
		'''
		pretty print this class
		'''
		return pprint.pformat(vars(self))

	############################################################################
	def getformatVersion(self):
		'''
		read the first few parameters from the header to see if this is a las file and what version
		'''
		curr = self.fileptr.tell()

		snifffmt = "<4sHHLHH8sBB"
		data = self.fileptr.read(struct.calcsize(snifffmt))
		s = struct.unpack(snifffmt, data)

		FileSignature 				= s[0].decode('utf-8').rstrip('\x00')
		FileSourceID 				= s[1]
		GlobalEncoding 				= s[2]
		ProjectIDGUIDdata1 			= s[3]
		ProjectIDGUIDdata2 			= s[4]
		ProjectIDGUIDdata3 			= s[5]
		ProjectIDGUIDdata4 			= s[6]
		VersionMajor 				= s[7]
		VersionMinor 				= s[8]

		self.fileptr.seek(curr, 0)

		return (FileSignature == "LASF", VersionMajor+VersionMinor/10)

	############################################################################
	def readhdr(self):
		'''
		read the las file header from disc
		'''
		islas, self.hdr.lasformat = self.getformatVersion()
		if self.hdr.lasformat == 1.2:
			data = self.fileptr.read(self.hdr.hdr12len)
			self.hdr.decodehdr(data)
		if self.hdr.lasformat == 1.4:
			data = self.fileptr.read(self.hdr.hdr14len)
			self.hdr.decodehdr(data)

	############################################################################
	def unpackpoints(self, records):
		'''
		the points read into the list need unpacking into the real world useful data
		'''
		for r in records:
			self.x.append((r[0] * self.hdr.Xscalefactor) + self.hdr.Xoffset)
			self.y.append((r[1] * self.hdr.Yscalefactor) + self.hdr.Yoffset)
			self.z.append((r[2] * self.hdr.Zscalefactor) + self.hdr.Zoffset)

	############################################################################
	def readpointrecords(self, recordsToRead=1):
		'''
		read the required number of records from the file
		'''
		data = self.fileptr.read(self.supportedformats[self.hdr.PointDataRecordFormat][1] * recordsToRead)
		result = []
		i = 0
		for r in range(recordsToRead):
			j = i + self.supportedformats[self.hdr.PointDataRecordFormat][1]
			result.append(struct.unpack(self.supportedformats[self.hdr.PointDataRecordFormat][0], data[i:j]))
			i = j
		
		return result

	############################################################################
	def readvariablelengthrecord(self):
		'''
		read a variable length record from the file
		'''
		vlrhdr14fmt = "<H16sHH32s"
		vlrhdr14len = struct.calcsize(vlrhdr14fmt)
		data = self.fileptr.read(vlrhdr14len)
		s = struct.unpack(vlrhdr14fmt, data)

		self.vlrReserved				 	= s[0]
		self.vlrUserid					 	= s[1]
		self.vlrrecordid				   	= s[2]
		self.vlrRecordLengthAfterHeader		= s[3]
		self.vlrDescription					= s[4]

		# now read the variable data
		self.vlrdata = self.fileptr.read(self.vlrRecordLengthAfterHeader)
		print (self.vlrdata)


###############################################################################
def createOutputFileName(path):
	'''Create a valid output filename. if the name of the file already exists the file name is auto-incremented.'''
	path	  = os.path.expanduser(path)

	if not os.path.exists(os.path.dirname(path)):
		os.makedirs(os.path.dirname(path))

	if not os.path.exists(path):
		return path

	root, ext = os.path.splitext(os.path.expanduser(path))
	dir	   = os.path.dirname(root)
	fname	 = os.path.basename(root)
	candidate = fname+ext
	index	 = 1
	ls		= set(os.listdir(dir))
	while candidate in ls:
			candidate = "{}_{}{}".format(fname,index,ext)
			index	+= 1

	return os.path.join(dir, candidate)


	############################################################################
def testreader(filename):
	'''
	sample read script so we can see how to use the code
	'''
	start_time = time.time() # time the process so we can keep it quick

	# filename = "C:/development/python/sample.las"
	# filename = "C:/development/python/version1.4_format0.las"
	# create a lasreader class and pass the filename
	r = lasreader(filename)

	# sniff the header so we can determine the file format
	islas, format = r.getformatVersion()
	if not islas:
		print ("this is not a las file, quitting...")
		exit()

	# now we know what version las file we are dealing with, we can read teh correct format.

	# read the header
	r.readhdr()

	# print some metadata about the reader
	print (r.hdr)

	# read the variable records
	for i in range(r.hdr.NumberofVariableLengthRecords):
		r.readvariablelengthrecord()

	# now find the start point for the point records
	r.seekPointRecordStart()
	# read the point data
	points = r.readpointrecords(64)
	# points = r.readpointrecords(r.hdr.Numberofpointrecords)
	
	# unpack from the native formmat into lists so we can do something with them
	r.unpackpoints(points)

	for i in range(len(r.x)):
		print ("%.3f, %.3f %.3f" % (r.x[i], r.y[i], r.z[i]))
		
	for p in points:
		print ("%.3f, %.3f %.3f" % ((p[0] * r.hdr.Xscalefactor) + r.hdr.Xoffset, (p[1] * r.hdr.Yscalefactor) + r.hdr.Yoffset, (p[2] * r.hdr.Zscalefactor) + r.hdr.Zoffset))

	print("Duration %.3fs" % (time.time() - start_time )) # time the process

	return

	############################################################################
def isBitSet(int_type, offset):
	'''testBit() returns a nonzero result, 2**offset, if the bit at 'offset' is one.'''
	mask = 1 << offset
	return (int_type & (1 << offset)) != 0


###############################################################################
if __name__ == "__main__":
		main()