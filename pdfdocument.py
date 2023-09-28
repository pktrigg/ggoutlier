###################################################################################################
#name:			pdfdocument.py
#created:		January 2019
#by:			paul.kennedy@guardiangeomatics.com
#description:	python module to create a standard report in native PDF format

#####################################################
#done
#27/2/20230 initial version
#the first table in the report needs to have a summary of:
# inputs and outputs  - done
# an image of the results
# the important metrics from th log file

#2DO
#####################################################
# nothing

import sys
import os
from argparse import ArgumentParser
from datetime import datetime
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import mm, cm
from reportlab.platypus import BaseDocTemplate, SimpleDocTemplate, PageTemplate, Paragraph, Spacer, Frame, Table, Image, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.lib import utils
from functools import partial
import PIL
import rasterio as rio
from rasterio.plot import show
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'shared'))
import fileutils
import cloud2tif

####################################################################################################
####################################################################################################
def bathyqcreport(logfilename, resultfolder):
	'''create an bathyqc report into PDF'''
	# resultfolder should be 5_grid
	# logfilename should be args.inputfolder\bathyqc.log

	if not os.path.exists(resultfolder):
		return

	outfilename = os.path.join(resultfolder, "BathyQCReport.pdf")
	outfilename = fileutils.createOutputFileName(outfilename)
	myreport = REPORT("BathyQC Report", outfilename)
	log(filename, "BathyQCReport.pdf")
	#parse the bathyqc log file and make a summary table
	if os.path.exists(logfilename):
		bathyqcreportsummary(myreport, logfilename )


####################################################################################################
def collectinformation(line, msgid, username, metrics):
	if msgid in line:
		line = line.replace(msgid,"")
		line = line.strip()
		metrics.append([username, line])
	return line

####################################################################################################
def bathyqcreportsummary(myreport, logfilename):

	if not os.path.exists(logfilename):
		return

	myreport.addtitle("BathyQC Principles")
	myreport.addspace()
	myreport.addspace()

	myreport.addparagraph("BathyQC is a tool developed by Guardian Geomatics to integrate many an unlimited quantity of of multibeam bathymetry survey lines, and generate a series of data products.  These products can be used to quality control in a consistent manner and when acceptable, deliver to the client.  The tool is rigorous and highly automated, producing consistent deliverable products.")
	myreport.addspace()

####################################################################################################
####################################################################################################
####################################################################################################
def reportsummary(myreport, GGOutlierlogfilename):

	if not os.path.exists(GGOutlierlogfilename):
		return

	#process the log file
	surveylines = []
	GGOutlierduration = ""
	surveyline = None
	status = []
	metrics = []
	metrics.append (["Inputs_Summary_Log", GGOutlierlogfilename])

	with open(GGOutlierlogfilename) as fp:
		for line in fp:
			line = line.replace("  ", " ")
			line = line.replace("  ", " ")
			line = line.replace("  ", " ")
			line = line.strip()
			line = line.lstrip()
			line = line.rstrip()

			collectinformation(line, "INFO:root:Username:", "Username", metrics)
			collectinformation(line, "INFO:root:Computer:", "Computer", metrics)
			collectinformation(line, "INFO:root:GGOutlier Version:", "QC_Version", metrics)
			collectinformation(line, "INFO:root:Writing outliers to:", "Outlier_Shape_Filename", metrics)
			collectinformation(line, "INFO:root:Created LAZ file of outliers:", "Outlier_LAZ_Filename", metrics)
			collectinformation(line, "INFO:root:Processing file:", "Input_Filename", metrics)
			collectinformation(line, "INFO:root:QC Duration:", "GGOutlier_Duration", metrics)
			collectinformation(line, "INFO:root:Depths loaded for quality control:", "Depths_Loaded_for_QC", metrics)
			collectinformation(line, "INFO:root:Points tagged for further evaluation:", "Machine_Learning_Candidates", metrics)	
			collectinformation(line, "INFO:root:Points outside specification:", "**Final_Outliers_Exceeding_Standard**", metrics)
			collectinformation(line, "INFO:root:QC to Survey Standard:", "Required_Survey_Standard", metrics)
			collectinformation(line, "INFO:root:Survey_Standard:", "Survey_Standard_Details", metrics)
			collectinformation(line, "INFO:root:EPSGCode for geodetic conversions:", "EPSG_Code", metrics)
			collectinformation(line, "INFO:root:Percentage outside specification:", "Percentage_Outside_Specification", metrics)
			collectinformation(line, "INFO:root:Points checked:", "Points_Checked", metrics)

			msg = "INFO:root:Created REGIONAL TIF file for IHO validation:"
			if msg in line:
				line = line.replace(msg,"")
				line = line.strip()
				metrics.append(["Regional_TIF", line])
				regionalfilename = line
		
			msg = "INFO:root:Created TXT file of outliers:"
			if msg in line:
				line = line.replace(msg,"")
				line = line.strip()
				metrics.append(["Outlier_TXT_Filename", line])
				outliertxtfilename = line

	totalpoints = 0

	#write out the per line stats...
	reportfilename = GGOutlierlogfilename + "_adjustment.txt"
	f = open(reportfilename, 'w')
	f.write("Item Value\n")
	for rec in metrics:
		f.write("%s: %s\n" % (rec[0], rec[1]))
	f.close()

	myreport.addspace()
	myreport.addtitle("GGOutlier Summary of Results")
	myreport.addspace()
	myreport.addtable(reportfilename)

	myreport.addtitle("What is an Outlier?")
	myreport.addspace()
	myreport.addparagraph("In bathymetry, an outlier typically refers to an isolated or anomalous depth measurement or feature on a seafloor depth map or chart. These outliers can be depths that are significantly different from the surrounding seafloor topography. Outliers might be caused by various factors such as errors in data collection, equipment malfunction, or unique geological features like wrecks, obstructions, seamounts or underwater volcanoes that stand out from the surrounding seabed. The scale of an outlier can be considerable. Identifying and understanding outliers in bathymetric data is important for accurate navigation, scientific research, and ocean exploration. Separating a real feature from noise is a complex issue. The final decision comes down to the skill and experience of the Surveyor In Charge. GGOutlier efficiently analyse and highlight outliers for validation.")

	myreport.addtitle("GGOutlier Principles")
	myreport.addspace()

	myreport.addparagraph("GGOutlier is a tool developed by Guardian Geomatics to Quality Control processed a multibeam bathymetry surface, and validate a processed depth surface against a standard such as those published by IHO SP44 or HIPP. The principle is similar in methodology to a traditional review by a surveyor-in-charge (SIC) process in which the SIC would review the depth surface by identifying outliers relative to its nearest neighours, determine if the outlier is significant and if so flag it forinvestigation.")
	myreport.addparagraph("GGOutlier primary purpose is to positively, rigorously identify each and every depth which is considered an outlier relative to the required total vertical uncertainty at that depth. This provides the SIC and client with full confidence that the quality of the depth surface meets the required specification and any remaining outliers are known, have been investigated and are considered features rather than noise.")
	myreport.addspace()
	myreport.addparagraph("Inliers are points which do meet the required specification for allowable total vertical uncertainty.")
	myreport.addparagraph("Outliers are points which do NOT meet the required specification for allowable Total Vertical Uncertainty (TVU).")
	myreport.addparagraph("Inputs are very simple. A depth surface (a floating point TIF file) and a IHO SP44 specification such as 'order1a', 'specialorder'.")
	myreport.addparagraph("A 'Regional Surface' is created using median depths of the nearest neighbours to each depth.")
	myreport.addparagraph("A 'TVU Surface' is created. This is the allowable TVU for the depth of each and every pixel.")
	myreport.addparagraph("A 'DeltaZ Surface is created. This is the difference between the regional surface and the depth surface")
	myreport.addparagraph("The DeltaZ values are then assessed against the TVU for that depth and either flagged as an outlier or accepted as within specification. The flagged depths are called 'outliers'.")
	myreport.addparagraph("Outliers are saved to a point cloud file and a shape file. The shape file contains the processed depth, the Regional Depth, the AllowableTVU, the DeltaZ (difference between Regional Depth and processed depth) and a field for Review/Approval by SIC.")
	myreport.addparagraph("GGOutlier does NOT modify the input file in any way. It is a read-only process.")
	myreport.addparagraph("GGOutlier will generate a QC report (this document) in order to enable rapid assessment of results.")
	myreport.addspace()
	myreport.addtitle("Role of Surveyor In Charge")
	myreport.addparagraph("The SIC role is essential. The resulting shape file identifies all depths which do NOT meet the IHO specificaiton.  These are either outliers missed in processing or depths on stepp slopes which inherently will not meet TVU specification due to gridding resolution limitations. it is the role of the SIC to review these flagged outliers and confirm they are valid or need to be passed back to the data processors for additional cleaning.")
	myreport.addparagraph("The SIC should edit the shape file attribute field to document each outlier has been reveiwed and approved.")
	myreport.addspace()
	myreport.addtitle("Role of Data Processor")
	myreport.addparagraph("The shape file can be loaded into the processing software (ag CARIS, Qimera) and used to guide the data processor to revisit the ungridded raw data points and re-evaluate underlying data and edit if required.")
	myreport.addparagraph("If additional edits are required, the DP shall regenerate the depth surface and rerun GGOutlier.")
	myreport.addspace()
	myreport.addtitle("Role of Client Representative")
	myreport.addparagraph("The shape file and Regional Depth will be delivered as part of a survey report. This can be used by the client to gain confidence all outliers have been reviewed by the SIC and there are no additional outliers in the depth surface.")
	myreport.addspace()
	myreport.addparagraph("Below is an example of how to consume the results from GGOutlier using GIS to analyse outliers which do not meet specification.")
	myreport.addspace()
	
	image = os.path.join(os.path.dirname(__file__), "GGOutliergis.png")
	myreport.addimage(image, 450)

	plt.ioff()
	dtm_dataset = rio.open(regionalfilename)
	NODATA = dtm_dataset.nodatavals[0]
	dtm_data = dtm_dataset.read(1)
	dtm_data[dtm_data > 10000] = 0
	dtm_data[dtm_data <= -999] = 0

	# fig, ax = plt.subplots(1, 1, figsize=(6,6))
	# dtm_map = show(dtm_dataset,title='Digital Terrain Model',ax=ax);
	# show(dtm_dataset,contour=True, ax=ax); #overlay the contours
	# im = dtm_map.get_images()[0]
	# fig.colorbar(im, label = 'Elevation (m)', ax=ax) # add a colorbar
	# ax.ticklabel_format(useOffset=False, style='plain') # turn off scientific notation
	
	#Use hillshade function on the DTM data array
	hs_data = cloud2tif.hillshade(dtm_data,315,10)

	# fig, ax = plt.subplots(1, 1, figsize=(6,6))
	ext = [dtm_dataset.bounds.left, dtm_dataset.bounds.right, dtm_dataset.bounds.bottom, dtm_dataset.bounds.top]
	# plt.imshow(hs_data,extent=ext)
	# plt.colorbar(); 
	# plt.set_cmap('RdYlGn')
	# plt.title('TEAK Hillshade')
	# ax=plt.gca(); ax.ticklabel_format(useOffset=False, style='plain') #do not use scientific notation 
	# rotatexlabels = plt.setp(ax.get_xticklabels(),rotation=90) #rotate x tick labels 90 degrees

	#Overlay transparent hillshade on DTM:
	SMALL_SIZE = 8
	MEDIUM_SIZE = 10
	BIGGER_SIZE = 12

	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
	# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

	plt.figure().set_figwidth(10)
	plt.figure().set_figheight(20)
	plt.rcParams['figure.figsize'] = [8, 8]
	fig, ax = plt.subplots(1, 1)
	ax = plt.gca()
	ax.set_aspect('equal')
	plt.gca().set_aspect('equal', adjustable='box')

	# plt.rcParams.update({'font.size': 6})
	# im1 = plt.plot(dtm_data,cmap='terrain'); 
	im1 = plt.imshow(dtm_data,cmap='terrain',extent=ext); 
	# cbar = plt.colorbar(); cbar.set_label('Elevation, m',rotation=270,labelpad=20)
	im2 = plt.imshow(hs_data,cmap='Greys',alpha=0.8,extent=ext); 
	# im2 = plt.plot(hs_data,cmap='Greys',alpha=0.5); 
	# plt.colorbar()
	# ax.ticklabel_format(useOffset=False, style='plain') #do not use scientific notation 
	# rotatexlabels = plt.setp(ax.get_xticklabels(),rotation=90) #rotate x tick labels 90 degrees
	# ax.set_aspect('equal', adjustable='box')
	plt.axis('off')
	# plt.grid('on'); 
	# plt.colorbar();

	f = open(outliertxtfilename, 'r')
	Lines = f.readlines()
	count = 0
	points = []
	# Strips the newline character
	for line in Lines:
		count += 1
		point = line.strip().split(',')
		points.append(point)
		plt.plot(float(point[0]), float(point[1]), marker="x", markersize=2.5, markeredgewidth=0.25, markeredgecolor="red", markerfacecolor="green")

	plt.grid()
	# plt.show()
	plt.title('Depth Surface With Outliers')
	# plt.axis('on')
	# plt.show()
	overviewimagefilename = regionalfilename + "_hillshade.png"
	plt.savefig(overviewimagefilename, bbox_inches='tight', dpi=640)
	
	myreport.addspace()
	myreport.addspace()
	myreport.addspace()
	myreport.addspace()
	myreport.addimage(overviewimagefilename, width=640)
	myreport.addparagraph("The image above is a screenshot of the outliers relative to a hillshade of the input file. This provides a visual summary of the outliers. The outliers are plotted as red crosses.")
	myreport.addparagraph("Outliers will most frequently be located at rocky outcrops so will be clustered. This is exactly as it should be. Stray, lonely outliers are the ones to closely analyse in GIS and CARIS to confirm if they are valid. Lonely outliers on flat terrain are items for concern.")
	myreport.addparagraph("END OF REPORT.")

	return

####################################################################################################
def GGOutlierreport(GGOutlierlogfilename, resultfolder):
	'''create an infinitpos QC report into PDF'''
	if not os.path.exists(resultfolder):
		return

	outfilename = os.path.join(resultfolder, "GGOutlierQCReport.pdf")
	outfilename = fileutils.createOutputFileName(outfilename)
	myreport = REPORT("GGOutlier QCReport", outfilename)

	#parse the GGOutlier log file and make a summary table
	if os.path.exists(GGOutlierlogfilename):
		reportsummary(myreport, GGOutlierlogfilename )

	myreport.save()
	myreport.viewpdf()

# ###################################################################################################
def addQCImage(myreport, f, fragment, notes):
		if fragment in os.path.basename(f).lower():
			requiredwidth = 512
			cmapfilename = findcmap(os.path.dirname(myreport.filename), fragment+"_cmap")
			outfilename = os.path.join(f + "_QC.png")
			image = myreport.compositeimage(f, requiredwidth, cmapfilename, 15 ,outfilename)
			
			myreport.addtitle("File: %s" % (os.path.basename(f)))
			for note in notes:
				myreport.addparagraph(note)
			myreport.addimage(image, requiredwidth/3)

# ###################################################################################################
def findcmap(folder, text):
	cmapname = ""
	matches = fileutils.findFiles2(False, folder, "*"+text+"*")
	for f in matches:
		cmapname = f
		return cmapname
	return cmapname

# ###################################################################################################
def main():

	parser = ArgumentParser(description='\n * generate a PDF report from GGOutlier Process one or many mission folders using GGOutlier.')
	parser.add_argument('-i', 			dest='inputfolder', action='store', 		default='.',			help='the root folder to find one more more mission folders. Pease refer to procedure for the mission folder layout - e.g. c:/mysurveyarea')
	
	args = parser.parse_args()

	if args.inputfolder == '.':
		args.inputfolder = os.getcwd()

	resultfolder = os.path.join(args.inputfolder, "8_cor").replace('\\','/')
	GGOutlierlogfilename = os.path.join(os.path.dirname(args.inputfolder), "GGOutlier.log").replace('\\','/')
	GGOutlierreport(GGOutlierlogfilename, resultfolder)

	# outfilename = "c:/temp/myfile.pdf"
	# outfilename = fileutils.createOutputFileName(outfilename)
	
	# myreport = REPORT("reportname", outfilename)
	# myreport.addheader("Report %s" % (outfilename))
	# myreport.addtitle("The quick brown fox jumped over the lazy dog.")
	# myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")
	# myreport.addparagraph("X")

	# myreport.addimage("guardian.png", 50, "the guardian logo")
	# myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")

	# myreport.newpage()
	# myreport.addheader("page2")

	# myreport.addimage("guardian.png", 50, "the guardian logo")
	# myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")

	# myreport.newpage()
	# myreport.addheader("Page 3")

	# myreport.save()
	# myreport.viewpdf()

# 	# move the origin up and to the left
# 	c.translate(mm,mm)
# 	# define a large font
# 	c.setFont("Helvetica", 10)
# 	# choose some colors
# 	c.setStrokeColorRGB(0.2,0.5,0.3)
# 	c.setFillColorRGB(1,0,1)
# 	# draw some lines
# 	c.line(0,0,0,10*mm)
# 	c.line(0,0,10*mm,0)
# 	# draw a rectangle
# 	c.rect(2*mm,2*mm,10*mm,15*mm, fill=1)
# 	# make text go straight up
# 	c.rotate(45)
# 	# change color
# 	# say hello (note after rotate the y coord needs to be negative!)
# 	c.drawString(100*mm,100*mm,"Hello World",)
# 	c.setFillColorRGB(0,0,0.77)
# 	c.rotate(-45)

# 	image = "guardian.png"
# 	x = 100
# 	y = 100
# 	c.drawImage(image, x,y, preserveAspectRatio=True, width=10*mm,height=10,mask=None)

# ###################################################################################################
# def header(c, title):
# 	headerheight = 25*mm
# 	top = c._pagesize[1]

# 	# move the origin up and to the left
# 	# c.translate(mm,mm)
# 	# define a font
# 	# choose some colors
# 	c.setStrokeColorRGB(0,0,0)
# 	c.setFillColorRGB(0,0,0)

# 	# set the title
# 	x = 10*mm
# 	y = top - 15*mm
# 	c.setFont("Helvetica", 18)
# 	c.drawString(x, y, title)

# 	# set the date
# 	x = 10*mm
# 	y = top - 22*mm
# 	c.setFont("Helvetica", 8)
# 	str = "Report Date: %s" % (datetime.now().strftime("%Y%m%d%H%M%S"))
# 	c.drawString(x, y, str)

# 	image = "guardian.png"
# 	x = c._pagesize[0]-(20*mm)
# 	y = top - 22*mm
# 	c.drawImage(image, x,y, preserveAspectRatio=True, width=15*mm, height=15*mm, mask=None)

# 	# draw header line
# 	c.line(5*mm, top-headerheight, c._pagesize[0]-(5*mm), c._pagesize[1]-headerheight)

##############################################################################
class REPORTSURVEYLINE:
	'''class to hold a group for reporting'''
	##############################################################################
	def __init__(self, name):
		self.group = ""
		words = name.split(" ")
		if len(words) > 0:
			self.name = words[0]

##############################################################################
class REPORT:
	'''class to create a PDF report'''
	##############################################################################
	def __init__(self, title, filename):
		self.filename = filename
		self.title = title
		styles = getSampleStyleSheet()
		self.normal = styles["Normal"]
		self.heading1 = styles['Heading1']
		# self.doc = SimpleDocTemplate(filename, pagesize=A4)
		# self.canvas = canvas.Canvas(filename, pagesize=A4)
		# self.cursor = self.canvas._pagesize[1]
		# self.leftmargin = self.canvas.leftMargin
		# self.rightmargin = self.canvas.rightMargin
		self.story = []

		self.doc = BaseDocTemplate(self.filename, pagesize=A4)
		frame = Frame(self.doc.leftMargin, self.doc.bottomMargin, self.doc.width, self.doc.height-0.5*cm, id='normal')
		# header_content = Paragraph("This is a multi-line header.  It goes on every page.  " * 8, self.normal)
		template = PageTemplate(id='header', frames=frame, onPage=partial(self.addheader, title=title))
		self.doc.addPageTemplates([template])

###################################################################################################
	def get_image(self, path, width=100*mm):
		img = utils.ImageReader(path)
		iw, ih = img.getSize()
		aspect = ih / float(iw)
		# scalefactor=min((self.doc.height/ih), (self.doc.width/iw))
		width=min(width, self.doc.width)
		
		height=min((width * aspect), self.doc.height - self.doc.topMargin - self.doc.bottomMargin)
		# width = self.doc.width * scalefactor
		# height= self.doc.height * scalefactor
		# return Image(path, width=iw, height=ih)
		return Image(path, width=width, height=height)

		# lowable <Image at 0x2baa2bbd390 frame=normal filename=GGOutliergis.png>(1918 x 1032) too large on page 3 in frame 'normal'(439.27559055118115 x 671.716535433071*) of template 'header'

###################################################################################################
	def compositeimage(self, filename, requiredwidth, legendfilename, legendwidth, outfilename):

		#the image we need to add a legend to
		PIL.Image.MAX_IMAGE_PIXELS = None
		img = PIL.Image.open(filename, 'r')
		img_w, img_h = img.size

		MAXWIDTH = requiredwidth
		ratio = MAXWIDTH/img.size[0]
		newimg = img.resize((int(img.size[0]*ratio), int(img.size[1]*ratio)), PIL.Image.ANTIALIAS)
		newimg_w, newimg_h = newimg.size

		#the new image for compositing into
		background = PIL.Image.new('RGBA', (newimg_w, newimg_h), (255, 255, 255, 255))
		bg_w, bg_h = background.size
	
		offset = ((bg_w - newimg_w) // 2, (bg_h - newimg_h) // 2)
		background.paste(newimg, offset)

		#the legend
		if os.path.exists(legendfilename):
			imglegend = PIL.Image.open(legendfilename, 'r')
			# imglegend_w, imglegend_h = imglegend.size
			MAXWIDTH = legendwidth
			ratio = MAXWIDTH/imglegend.size[0]
			newlegend = imglegend.resize((int(imglegend.size[0]*ratio), int(imglegend.size[1]*ratio)), PIL.Image.ANTIALIAS)
			background.paste(newlegend, (0,0))

		background.save(outfilename)
		return outfilename

###################################################################################################
	def addimagetable(self, filename, width, legendfilename, legendheight):
		'''write an image'''

		self.story.append(self.get_image(filename, width=40*mm))

###################################################################################################
	def addimage(self, filename, width=320, height=640):
		'''write an image'''

		self.story.append(self.get_image(filename, width=width))
		# image = Image(filename, width=height, height=height)
		# image = Image(filename, width=height, height=height, preserveAspectRatio=True)
		# self.story.append(image)

		# x = ((self.docrightmargin - self.leftmargin) / 2 ) + self.leftmargin
		# y = self.cursor - (height*1.5)
		# self.canvas.drawImage(filename, x,y, preserveAspectRatio=True, width=height*mm, height=height*mm, mask=None, anchorAtXY=True, anchor='c')

		# self.setcursor(height*3.2)

		# if len(label) == 0:
		# 	return

		# self.canvas.drawCentredString(x, self.cursor, label)
		# self.setcursor(10)


	###################################################################################################
	def addheader(self,canvas, doc, title):
		canvas.saveState()

		top = doc.height + doc.topMargin + (5 * mm)
		#do the title
		header_content = Paragraph(title, self.heading1)
		w, h = header_content.wrap(doc.width, doc.topMargin)
		header_content.drawOn(canvas, doc.leftMargin, top)
		
		
		# # set the date
		canvas.setFont("Helvetica", 8)
		str = "Creation Date: %s" % (datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
		canvas.drawString(doc.leftMargin, top - 5*mm, str)

		# set the user
		canvas.setFont("Helvetica", 6)
		str = "User : %s" % (os.getlogin())
		canvas.drawString(doc.leftMargin + 130*mm, top - 5*mm, str)

		image = os.path.join(os.path.dirname(__file__), "guardian.png")
		x = doc.width + (15 * mm)
		y = top - 6*mm
		canvas.drawImage(image, x,y, preserveAspectRatio=True, width=15*mm, height=15*mm, mask=None)

		# # draw header line
		y = top - 7*mm
		canvas.line(doc.leftMargin, y, doc.width + (30*mm), y)

		canvas.restoreState()

	###################################################################################################
	# def setcursor(self, dy):
	# 	self.cursor -= dy

	# ###################################################################################################
	# def resetcursor(self):
	# 	self.cursor = self.canvas._pagesize[1]

	###################################################################################################
	def viewpdf(self):
		os.startfile(self.filename, 'open')

	###################################################################################################
	def save(self):
		self.doc.build(self.story)
		# self.doc.save()

	# ###################################################################################################
	def newpage(self):
		self.story.append(PageBreak())
	# 	self.canvas.showPage()
	# 	self.resetcursor()

	###################################################################################################
	def addparagraph(self, text):
		'''write a paragraph of text '''

		# text.append(Paragraph("This is line %d." % i, styleN))

		styles = getSampleStyleSheet()
		style = styles["Normal"]
		# story = [Spacer(1,2*mm)]
		p = Paragraph(text, style)
		self.story.append(p)
		# Story.append(Spacer(1,0.2*mm))
		# self.canvas.build(Story)

		# width = self.rightmargin - self.leftmargin
		# # pixperchar = style.fontSize
		# charsperline = width / style.fontSize * .15
		# height = max(35, len(text) / charsperline)
		# f = Frame(self.leftmargin, self.cursor-height, width, height, leftPadding=0, bottomPadding=1, rightPadding=0, topPadding=0, showBoundary=1 )
		# f.addFromList(story, self.canvas)

		# # self.canvas.setFont("Helvetica", 12)
		# # self.canvas.drawString(self.leftmargin, self.cursor, text)

		# self.setcursor(height)

###################################################################################################
	def addspace(self, height=1):
		self.story.append(Spacer(1,float(height)*mm))

###################################################################################################
	def addtable(self, filename):


		# data= [['00', '01', '02', '03', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04', '04'],
		# ['10', '11', '12', '13', '14'],
		# ['20', '21', '22', '23', '24'],
		# ['30', '31', '32', '33', '34']]
		# t=Table(data)
		# t.setStyle(TableStyle([('BACKGROUND',(1,1),(-2,-2),colors.green),
		# ('TEXTCOLOR',(0,0),(1,-1),colors.red)]))
		# self.story.append(t)
		# return


		data = []
		with open(filename) as fp:
			for line in fp:
				line = line.strip()
				line = line.lstrip()
				line = line.rstrip()
				line = line.replace("  ", " ")
				line = line.replace("  ", " ")
				line = line.replace("  ", " ")
				line = line.replace("  ", " ")
				line = line.replace("  ", " ")
				line = line.replace("  ", " ")
				line = line.replace("  ", " ")
				words = line.split(" ")
				data.append(words)

		t=Table(data, rowHeights=None, style=None, splitByRow=1, repeatRows=0, repeatCols=0, rowSplitRange=None, spaceBefore=None, spaceAfter=None)
		t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,0), "Helvetica-Bold")]))
		t.setStyle(TableStyle([('FONTSIZE',(0,0),(-1,-1), 6)]))
		t.setStyle(TableStyle([('LEFTPADDING',(0,0),(-1,-1), 2)]))
		t.setStyle(TableStyle([('RIGHTPADDING',(0,0),(-1,-1), 2)]))
		t.setStyle(TableStyle([('TOPPADDING',(0,0),(-1,-1), 2)]))
		t.setStyle(TableStyle([('BOTTOMPADDING',(0,0),(-1,-1), 2)]))
		t.setStyle(TableStyle([('INNERGRID',(0,0),(-1,-1), 0.25, colors.black)]))
		t.setStyle(TableStyle([('BOX',(0,0),(-1,-1), 0.25, colors.black)]))
		t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),colors.lightblue), ('TEXTCOLOR',(0,0),(0,0),colors.black)]))
		# t.setStyle(TableStyle([('BACKGROUND',(1,1),(-2,-2),colors.green),
		# ('TEXTCOLOR',(0,0),(1,-1),colors.red)]))
		self.story.append(t)

###################################################################################################
	def addtitle(self, text):
		'''write a title of text '''

		styles = getSampleStyleSheet()
		style = styles["Heading2"]
		p = Paragraph(text, style)
		self.story.append(p)
		# Story.append(Spacer(1,0.2*mm))
		# self.canvas.build(Story)

		# width = self.rightmargin - self.leftmargin
		# # pixperchar = style.fontSize
		# charsperline = width / style.fontSize * .15
		# height = max(50, len(text) / charsperline)
		# f = Frame(self.leftmargin, self.cursor-height, width, height, leftPadding=0, bottomPadding=1, rightPadding=0, topPadding=0, showBoundary=0 )
		# f.addFromList(story, self.canvas)

		# # self.canvas.setFont("Helvetica", 12)
		# # self.canvas.drawString(self.leftmargin, self.cursor, text)

		# self.setcursor(height)

	# ###################################################################################################
	# def oldaddheader(self, title):
	# 	headerheight = 25*mm
	# 	top = self.canvas.height

	# 	# move the origin up and to the left
	# 	# c.translate(mm,mm)
	# 	# define a font
	# 	# choose some colors
	# 	# self.canvas.setStrokeColorRGB(0,0,0)
	# 	# self.canvas.setFillColorRGB(0,0,0)

	# 	# set the title
	# 	x = 10*mm
	# 	y = top - 15*mm
	# 	self.canvas.setFont("Helvetica-Bold", 14)

	# 	self.canvas.drawString(x, y, title)

	# 	# set the date
	# 	x = 10*mm
	# 	y = top - 20*mm
	# 	self.canvas.setFont("Helvetica", 8)
	# 	str = "Creation Date: %s" % (datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
	# 	self.canvas.drawString(x, y, str)

	# 	# set the user
	# 	x = 150*mm
	# 	y = top - 20*mm
	# 	self.canvas.setFont("Helvetica", 8)
	# 	str = "User : %s" % (os.getlogin())
	# 	self.canvas.drawString(x, y, str)

	# 	# set the user
	# 	x = 10*mm
	# 	y = top - 23*mm
	# 	self.canvas.setFont("Helvetica", 8)
	# 	str = "Filename : %s" % (self.filename)
	# 	self.canvas.drawString(x, y, str)

	# 	image = "guardian.png"
	# 	x = self.canvas._pagesize[0]-(20*mm)
	# 	y = top - 22*mm
	# 	self.canvas.drawImage(image, x,y, preserveAspectRatio=True, width=15*mm, height=15*mm, mask=None)

	# 	# draw header line
	# 	self.canvas.line(5*mm, top-headerheight, self.canvas._pagesize[0]-(5*mm), self.canvas._pagesize[1]-headerheight)

	# 	self.setcursor(headerheight)




##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# class REPORTORI:
# 	'''class to create a PDF report'''
# 	##############################################################################
# 	def __init__(self, filename):
# 		self.filename = filename
# 		self.title = ""
# 		self.doc = SimpleDocTemplate(filename, pagesize=A4)
# 		self.canvas = canvas.Canvas(filename, pagesize=A4)
# 		self.cursor = self.canvas._pagesize[1]
# 		self.leftmargin = 10 * mm
# 		self.rightmargin = self.canvas._pagesize[0] - (10 * mm)

# 	###################################################################################################
# 	def setcursor(self, dy):
# 		self.cursor -= dy

# 	###################################################################################################
# 	def resetcursor(self):
# 		self.cursor = self.canvas._pagesize[1]

# 	###################################################################################################
# 	def viewpdf(self):
# 		os.startfile(self.filename, 'open')

# 	###################################################################################################
# 	def save(self):
# 		self.canvas.save()

# 	###################################################################################################
# 	def newpage(self):
# 		self.canvas.showPage()
# 		self.resetcursor()

# 	###################################################################################################
# 	def addparagraph(self, text):
# 		'''write a paragraph of text '''

# 		styles = getSampleStyleSheet()
# 		style = styles["Normal"]
# 		story = [Spacer(1,2*mm)]
# 		p = Paragraph(text, style)
# 		story.append(p)
# 		# Story.append(Spacer(1,0.2*mm))
# 		# self.canvas.build(Story)

# 		width = self.rightmargin - self.leftmargin
# 		# pixperchar = style.fontSize
# 		charsperline = width / style.fontSize * .15
# 		height = max(35, len(text) / charsperline)
# 		f = Frame(self.leftmargin, self.cursor-height, width, height, leftPadding=0, bottomPadding=1, rightPadding=0, topPadding=0, showBoundary=1 )
# 		f.addFromList(story, self.canvas)

# 		# self.canvas.setFont("Helvetica", 12)
# 		# self.canvas.drawString(self.leftmargin, self.cursor, text)

# 		self.setcursor(height)

# ###################################################################################################
# 	def addtable(self, filename):

# 		data = []
# 		with open(filename) as fp:
# 			line = fp.readline()
# 			cnt = 1
# 			while line:
# 				print("Line {}: {}".format(cnt, line.strip()))
# 				line = fp.readline()
# 				words = line.split(" ")
# 				data.append([words])
# 				cnt += 1

# 		styles = getSampleStyleSheet()
# 		style = styles["Normal"]
# 		story = [Spacer(1,2*mm)]
# 		t = Table(data, colWidths=10, rowHeights=None, style=None, splitByRow=1, repeatRows=0, repeatCols=0, rowSplitRange=None, spaceBefore=None, spaceAfter=None)
# 		# t = Table(data)

# 		width = 200
# 		height = 200
# 		# self.canvas.wrapOn(t, width, height)
# 		t.drawOn(self.canvas, self.leftmargin, self.cursor-200)
		
# 		# col_widths = 10
# 		# story.append(Table(data, colWidths=col_widths))

# 		# # self.canvas.append(t)

# 		# height = 200
# 		# t.drawOn(self.canvas, self.leftmargin, self.cursor-height, _sW=0)

# 		# story.append(t)
# 		# width = self.rightmargin - self.leftmargin
# 		# charsperline = width / style.fontSize * .15
# 		# height = max(135, len(data) / charsperline)
# 		# f = Frame(self.leftmargin, self.cursor-height, width, height, leftPadding=0, bottomPadding=1, rightPadding=0, topPadding=0, showBoundary=0 )
# 		# f.addFromList(story, self.canvas)
# 		# f = Frame(self.leftmargin, self.cursor-height, width, height, leftPadding=0, bottomPadding=1, rightPadding=0, topPadding=0, showBoundary=1 )
# 		# f.addFromList(story, self.canvas)

# ###################################################################################################
# 	def addtitle(self, text):
# 		'''write a title of text '''

# 		styles = getSampleStyleSheet()
# 		style = styles["Heading2"]
# 		story = [Spacer(1,2*mm)]
# 		p = Paragraph(text, style)
# 		story.append(p)
# 		# Story.append(Spacer(1,0.2*mm))
# 		# self.canvas.build(Story)

# 		width = self.rightmargin - self.leftmargin
# 		# pixperchar = style.fontSize
# 		charsperline = width / style.fontSize * .15
# 		height = max(50, len(text) / charsperline)
# 		f = Frame(self.leftmargin, self.cursor-height, width, height, leftPadding=0, bottomPadding=1, rightPadding=0, topPadding=0, showBoundary=0 )
# 		f.addFromList(story, self.canvas)

# 		# self.canvas.setFont("Helvetica", 12)
# 		# self.canvas.drawString(self.leftmargin, self.cursor, text)

# 		self.setcursor(height)

# 	###################################################################################################
# 	def addheader(self, title):
# 		headerheight = 25*mm
# 		top = self.canvas._pagesize[1]

# 		# move the origin up and to the left
# 		# c.translate(mm,mm)
# 		# define a font
# 		# choose some colors
# 		self.canvas.setStrokeColorRGB(0,0,0)
# 		self.canvas.setFillColorRGB(0,0,0)

# 		# set the title
# 		x = 10*mm
# 		y = top - 15*mm
# 		self.canvas.setFont("Helvetica-Bold", 14)

# 		self.canvas.drawString(x, y, title)

# 		# set the date
# 		x = 10*mm
# 		y = top - 20*mm
# 		self.canvas.setFont("Helvetica", 8)
# 		str = "Creation Date: %s" % (datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
# 		self.canvas.drawString(x, y, str)

# 		# set the user
# 		x = 150*mm
# 		y = top - 20*mm
# 		self.canvas.setFont("Helvetica", 8)
# 		str = "User : %s" % (os.getlogin())
# 		self.canvas.drawString(x, y, str)

# 		# set the user
# 		x = 10*mm
# 		y = top - 23*mm
# 		self.canvas.setFont("Helvetica", 8)
# 		str = "Filename : %s" % (self.filename)
# 		self.canvas.drawString(x, y, str)

# 		image = "guardian.png"
# 		x = self.canvas._pagesize[0]-(20*mm)
# 		y = top - 22*mm
# 		self.canvas.drawImage(image, x,y, preserveAspectRatio=True, width=15*mm, height=15*mm, mask=None)

# 		# draw header line
# 		self.canvas.line(5*mm, top-headerheight, self.canvas._pagesize[0]-(5*mm), self.canvas._pagesize[1]-headerheight)

# 		self.setcursor(headerheight)

# ###################################################################################################
# 	def addimage(self, filename, height, label=""):
# 		'''write an image'''

# 		x = ((self.rightmargin - self.leftmargin) / 2 ) + self.leftmargin
# 		y = self.cursor - (height*1.5)
# 		self.canvas.drawImage(filename, x,y, preserveAspectRatio=True, width=height*mm, height=height*mm, mask=None, anchorAtXY=True, anchor='c')

# 		self.setcursor(height*3.2)

# 		if len(label) == 0:
# 			return

# 		self.canvas.drawCentredString(x, self.cursor, label)
# 		self.setcursor(10)

# ###################################################################################################
# def main2():
# 	outfilename = "c:/temp/myfile.pdf"
# 	outfilename = fileutils.createOutputFileName(outfilename)
# 	resultfolder = os.path.dirname(outfilename)

# 	styles = getSampleStyleSheet()
# 	styleN = styles['Normal']
# 	# styleH = styles['Heading1']

# 	myreport = REPORT("GGOutlier QC Report %s" % (os.path.dirname(resultfolder)), outfilename)

# 	myreport.addparagraph("hello")
# 	myreport.addparagraph("hello2")
# 	myreport.addimage("guardian.png", 10)
# 	myreport.addparagraph("hello3")
# 	myreport.addparagraph("hello")
# 	myreport.addparagraph("hello2")
# 	myreport.addimage("guardian.png", 10)
# 	myreport.addparagraph("hello3")
# 	myreport.addparagraph("hello")
# 	myreport.addparagraph("hello2")
# 	myreport.addimage("guardian.png", 10)
# 	myreport.addparagraph("hello3")
# 	myreport.addparagraph("hello")
# 	myreport.addparagraph("hello2")
# 	myreport.addimage("guardian.png", 10)
# 	myreport.addparagraph("hello3")
# 	myreport.addparagraph("hello")
# 	myreport.addparagraph("hello2")
# 	myreport.addimage("guardian.png", 10)
# 	myreport.addparagraph("hello3")
# 	myreport.addparagraph("hello")
# 	myreport.addparagraph("hello2")
# 	myreport.addimage("guardian.png", 10)
# 	myreport.addparagraph("hello3")
# 	# story = []
# 	myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")
# 	myreport.addspace(2)

# 	myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")
# 	myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")
# 	myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")
# 	myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")
# 	myreport.doc.build(myreport.story)
# 	# myreport.doc.build(myreport.story)
	
# 	# for i in range(111):
# 	# 	story.append(Paragraph("This is line %d." % i, styleN))
# 	# myreport.doc.build(story)

# 	# myreport.addtitle("The quick brown fox jumped over the lazy dog.")
# 	# myreport.addparagraph("X")

# 	# myreport.addimage("guardian.png", 50, "the guardian logo")
# 	# myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")

# 	# myreport.newpage()
# 	# myreport.addheader("page2")

# 	# myreport.addimage("guardian.png", 50, "the guardian logo")
# 	# myreport.addparagraph("This is a bad world.  we need to take care of The BaseDocTemplate class implements the basic machinery for document formatting. An instance of the class contains a list of one or more PageTemplates that can be used to describe the layout of information on a single page. The build method can be used to process a list of Flowables to produce a PDF document.")

# 	# myreport.newpage()
# 	# myreport.addheader("Page 3")

# 	myreport.save()
# 	myreport.viewpdf()


###################################################################################################
if __name__ == "__main__":

	main()
	# resultfolder = "E:/projects/GGOutlier/A14/IP_Result_20200725013252/8_cor"
	# GGOutlierlogfilename = "E:/projects/GGOutlier/A14/GGOutlier.log"

	# GGOutlierreport(GGOutlierlogfilename, resultfolder)
