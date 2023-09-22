#name:			ggmbesstandard
#created:		July 2017
#by:			p.kennedy@guardiangeomatics.com
#description:	python module to represent MBES data STANDARDS
 
import math
import pprint

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
	def gettvatdept(self, depth):
		'''TVU(d) = sqrt((a*a) + ( b * d)^2)'''
		tvud = math.sqrt((self.depthtvu_a * self.depthtvu_a) + (self.depthtvu_b * depth)**2)
		return tvud
	###############################################################################
	def __str__(self):
		return pprint.pformat(vars(self))
