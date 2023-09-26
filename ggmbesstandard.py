#name:			ggmbesstandard
#created:		July 2017
#by:			p.kennedy@guardiangeomatics.com
#description:	python module to represent MBES data STANDARDS
 
import math
import pprint
import rasterio
import numpy as np

###############################################################################
class sp44:
	'''used to hold the metadata associated with an IHO MBES standard.'''
	def __init__(self):
		self.name			= ""
		self.longitude 			= 0
		self.standards = []

		self.standards.append(standard("order2", 1.0, 0.023))
		self.standards.append(standard("order1b", 0.5, 0.013))
		self.standards.append(standard("order1a", 0.5, 0.013))
		self.standards.append(standard("specialorder", 0.25, 0.0075))
		self.standards.append(standard("exclusiveorder", 0.15, 0.0075))

		self.standards.append(standard("hipp1", 0.25, 0.0075))
		self.standards.append(standard("hipp2", 0.5, 0.013))
		self.standards.append(standard("hippassage", 1.0, 0.023))

	###############################################################################
	def __str__(self):
		return pprint.pformat(vars(self))

###############################################################################
	def getordernames(self):
		msg = []
		for rec in self.standards:
			msg.append(rec.name)
		return msg
	
###############################################################################
	def loadstandard(self, namerequired):
		for rec in self.standards:
			if namerequired in rec.name:
				return rec

###############################################################################
class standard:
	'''used to hold the metadata associated with an IHO MBES standard.'''
	def __init__(self, name, depthtvu_a, depthtvu_b ):
		self.name					= name
		self.depthtvu_a 			= depthtvu_a
		self.depthtvu_b 			= depthtvu_b

	###############################################################################
	def gettvuat(self, depth):
		'''TVU(d) = sqrt((a*a) + ( b * d)^2)'''
		tvud = math.sqrt((self.depthtvu_a * self.depthtvu_a) + (self.depthtvu_b * depth)**2)
		return tvud
	###############################################################################
	def details(self):
		msg = "Name:" + self.name + ",a=" + str(self.depthtvu_a) + ",b=" + str(self.depthtvu_b) + ",TVU(d)=sqrt((a*a)+(b*d)^2)"
		return msg

	###############################################################################
	def computeTVUSurface(self, filename, outfilename):
		'''compute the TVU for a surface array'''
		with rasterio.open(filename) as src:
			array = src.read(1)
			profile = src.profile
			NODATA = src.nodatavals[0]

		#now compute the TVU for the entire surface using numpy array mathmatics so its fast
		#preserve the NODATA value
		# array[array==NODATA] = 0
		arrayTVU = np.multiply (array, self.depthtvu_b)
		arrayTVU = np.square (arrayTVU, arrayTVU)
		arrayTVU = np.add (arrayTVU, (self.depthtvu_a*self.depthtvu_a))
		arrayTVU = np.sqrt(arrayTVU)
		# arrayTVU[arrayTVU== self.depthtvu_a] = NODATA

		#reset the nodata value...
		arrayTVU[arrayTVU> NODATA] = NODATA

		# Write to tif, using the same profile as the source
		with rasterio.open(outfilename, 'w', **profile) as dst:
			dst.write_band(1, arrayTVU)

		return outfilename

	###############################################################################
	def computeTVUBarometer(self, allowabletvufilename, uncertaintyfilename, outfilename):
		'''compute the TVU barometric pressure. A low pressure represents where the TVU for a survey point is well within specificaiton.  As high pressure is where the TVU is almost using all the allowable TVU'''
		with rasterio.open(allowabletvufilename) as allowedsrc:
			allowedarray = allowedsrc.read(1)
			allowedprofile = allowedsrc.profile
			allowedNODATA = allowedsrc.nodatavals[0]
			allowedarray[allowedarray==allowedNODATA] = 9999

		allowedsrc.close()
		with rasterio.open(uncertaintyfilename) as uncertaintysrc:
			uncertaintyarray = uncertaintysrc.read(1)
			uncertaintyprofile = uncertaintysrc.profile
			uncertaintyNODATA = uncertaintysrc.nodatavals[0]
			uncertaintyarray[uncertaintyarray==uncertaintyNODATA] = 0
		uncertaintysrc.close()
	
		#now compute the TVU barometric pressure for the entire surface using numpy array mathmatics so its fast
		# the TVUBAROMETER is the percentage of the allowable uncertainty compared to the actual uncertainty as computed by CUBE (or other software)
		# eg if the allowable uncertainty is 0.5m and the actual uncertainty is 0.25m then the TVUBAROMETER is 50%
		# eg if the allowable uncertainty is 0.5m and the actual uncertainty is 0.75m then the TVUBAROMETER is 150%
		# eg if the allowable uncertainty is 0.5m and the actual uncertainty is 0.5m then the TVUBAROMETER is 100%
		tvubarometerarray = np.divide (uncertaintyarray, allowedarray)
		tvubarometerarray = np.multiply (tvubarometerarray, 100)

		# Write to tif, using the same profile as the source
		with rasterio.open(outfilename, 'w', **allowedprofile) as dst:
			dst.write_band(1, tvubarometerarray)

		return outfilename

