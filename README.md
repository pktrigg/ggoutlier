# kmallclean.py
KmallClean uses machine learning to clean a point cloud inside a kmall file and write out a new kmall file which has ALL the original data but with the rejected soundings flagged as rejected. The new file can then be imported into CARIS / Qimera for further processing.  The rejected beams can always be re-accepted in the post processing package.  

*	The filter is a pre-processing filter to be used BEFORE formal processing in a toold such as CARIS/QIMERA
*	The filtering is an automated process which can be overloaded by humans if required.  
*	The algorithm treats each file independent of other files.
*	The results are stored to a NEW KMALL file.  the original KMALL file is not modified in any way
*	The algorithm will iterate through the data to and analyse the level of noise present in the kmall file. By default it will adjust the filters to only reject approximately 1 percent *	of the point cloud data.
*	The NEW KMALL file is created with _C.KMALL extension.  The file will contain identical records, depths, posiitons etc.  the only values changed are the detect quality flags in the MRZ records.  these are decoded by processing tools such as CARIS/QIMERA when importing a file for processing.  
*	The point cloud data is categorised into 'outliers' and 'inliers'  For convenience, an outlier .SHP file is created.  This can be dropped into GIS for review.


***Inputs***  
* KMALL files
* user settings if required. 

***Command Line Parameters***  
*	-epsg	Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326'[Default: 0 == automatic UTM zone calculation]  
*	-i		Input filename/folder to process.  
*	If nothing is specificed the current folder is used.  
*	If a file is specified it will be processed.  
*	If a folder is specified all files in the folder will br processed.  
*	Wildcards are not supported.  
*	-c			clip outer beams each side to this max angle. Set to -1 to disable [Default: -1]' This is not likeyl to be used but is good for testing  
*	-cpu		number of cpu processes to use in parallel. [Default: 0, all cpu]'  
*	-odir		Specify a relative output folder e.g. -odir GIS' [Deafult : ""]  
*	-n			Specify the number of nearest neighbours points to use.  More points means more data will be rejected. ADVANCED ONLY [Default:5]'  
*	-p			Specify the approximate percentage of data to remove.  The engine will analyse the data and learn what filter settings are appropriate for your waterdepth and data quality. This is the MOST IMPORTANT (and only) parameter to consider spherical radius to find the nearest neightbours. [Default:1.0]'  
*	-z			Specify the ZScale to accentuate the depth difference ove the horizontal distance between points. Think of this as how you exxagerate teh vertical scale in a swath editor to more easily spot the outliers. [Default:10]'  
*	-debug		Specify the number of pings to process.  good only for debugging. [Default:-1]'  

***Outputs***
* A <filename>_R.TIF file which contains a quick and dirty TIF file of ALL depths prior to cleaning
* A <filename>_Inlier.TIF file which contains a quick and dirty TIF file of accepted depths following cleaning
* A NEW KMALL file with inliers and outliers flagged in the MRZ beamDetection field.
* A <filename>_outlier.SHP file which contains a multi-point shape file of the detected outliers
* A <filename>_outlier.txt file which contains the rejected points so you can review and plot if needed
* A <filename>_inlier.tif file which contains a quick and dirty TIF file of floating point accepted depths

***How to run the script***  
*	python c:\ggtools\kmallclean\kmallclean.py

***Example Output***  
python C:\ggtools\kmallclean\kmallclean.py -epsg 0 -i C:/sampledata/SI1026_A/Special_Order_Selection/SpecialOrder -cpu 1 -odir Z50_N5_P1 -debug 10000000 -n 5 -p 1 -z 50
Output Folder: C:/sampledata/SI1026_A/Special_Order_Selection/SpecialOrder\Z50_N5_P1  
Loading Point Cloud...  
EPSGCode for geodetic conversions: 32751  
Extracting Point Cloud: [####################] 99.96%  
Depths Loaded for cleaning: 1256448  
Creating tif file... C:/sampledata/SI1026_A/Special_Order_Selection/SpecialOrder\Z50_N5_P1\A_Infill_001_0509_20221009_024512.kmall_R.txt_R.tif  
Percentage rejection 0.00  
Filter level increasing to reject a few more points...  
Percentage rejection 0.00  
Filter level increasing to reject a few more points...  
Percentage rejection 0.07  
Filter level increasing to reject a few more points...  
Percentage rejection 0.67  
Filter level increasing to reject a few more points...  
Percentage rejection 39.44  
Filter level decreasing to reject a few less points...  
Percentage rejection 4.12  
Filter level decreasing to reject a few less points...  
Percentage rejection 1.46  
Filter level decreasing to reject a few less points...  
Percentage rejection 0.95  
Points accepted: 1244474.00  
Points rejected: 11974.00  
Creating tif file...  
C:/sampledata/SI1026_A/Special_Order_Selection/SpecialOrder\Z50_N5_P1\A_Infill_001_0509_20221009_024512.kmall_Inlier.txt.tif  
Writing NEW KMALL file C:/sampledata/SI1026_A/Special_Order_Selection/SpecialOrder\Z50_N5_P1\A_Infill_003_0504_20221009_023541_CLEANED.kmall  
Writing cleaned data: [####################] 99.91%  
