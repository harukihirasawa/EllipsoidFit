# Colloidal Ellipsoid Fitting
# Using Nima Moshtagh's MVEE algorithm
# ASSUMES 3D IMAGES:
# WILL FIX TO HANDLE 2D later

from PIL import Image
import math as m
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import tifffile as tiff
from scipy import ndimage
from scipy import spatial
import MVEE
from skimage import filters
from skimage import morphology
import string
import sys

'''
FUNC: EllipsoidFit
PURPOSE: Primary overhead function
INPUT:	image [str or np.ndarray] : 	
		either a string for the image name or a 3D numpy array representing 
		the image.
	threshold [int or str] : 	
		either a given thresholding value or a string indicating the type of 
		thresholding to be performed. Currently available: global otsu, 
		local otsu, yin, li.
		local otsu also requires a radius so the input is formatted "local otsu #". Where # is the radius
		DEFAULT: global otsu
	voxelscale [list]: 
		List indicating the size of the voxels in !![x,y,z]!!
	outfile [str] :
		String specifying the name of the files that the output data are
		saved to
		DEFAULT: image name (if image is a string) or the current time and 
		date.
	minsize [int] : 
		The minimum size of a object to which a ellipsoid will be fitted.
		DEFAULT: 10
	tol [float] :
		The criterion used to end the iterative ellipsoid fitting.
		DEFAULT: 0.01
	output [bool] :
		Toggles file output
		DEFAULT: True
	display [bool] :
		Toggle graph GUI display
		DEFAULT: True
	hull [bool]:
		Toggles fitting to the convex hull of the objects
		DEFAULT: True
	fittype [string]:
		Key for choosing type of bounding volume
		DEFAULT: ellipsoid
	smooth [str]:
		Specifies smoothing type (median,mean,none) currently always uses
		radius 1 smoothing
		DEFAULT: median
	
OUTPUT:	data [list]: 
		A list of dictionaries containing information on each ellipsoid
	if output is True, program will write a number of files (SEE: 
	graphrods, datafileprint, graphdata)
NOTES: 	WIP
	Plan to add more outputs
	Currently fits all thresholded regions over the minimum size, so rods that 
	are near one another can and often will get fitted together. Finding and 
	separating such objects is WIP.
	Due to the way tifffile reads in the stacks, the pixel coords used are [z,x,y]
'''
def EllipsoidFit(image,threshold = "otsu global",voxelscale=[1,1,1],
	outfile=None,minsize = 30,tol = 0.01,
	outputdata = True,output3D = False, outputhist = True,
	display3D=False, displayhist=False, hull = True,
	fittype = "ellipsoid",smooth = "median"):
	
	starttime = time.time()
	# If outfile is not provided, set automatically
	if outfile is None:
		if isinstance(image,str):
			outfile = image.split(".")[0]
		else:
			outfile = time.asctime(time.localtime(time.time()))
	# Read image into array
	if isinstance(image,str): 
		imgtiff = tiff.TIFFfile(image)
		image = imgtiff.asarray()

	ishape = image.shape #image shape
	if smooth == "median":
		for i in range(ishape[0]):
			image[i] = filters.rank.median(image[i],morphology.disk(1))
	elif smooth == "mean":
		for i in range(ishape[0]):
			image[i] = filters.rank.mean(image[i],morphology.disk(1))
	print image.shape
	# Rescale voxels to be cubic
	if len(set(voxelscale)) > 1:
		image = cubifyvoxels(image,voxelscale)
	print image.shape

	# Binarize images and label objects
	binimg, thresh = thresholding(image,threshold)
	lbimg,objnum = ndimage.label(binimg)
	print "Objects: ",objnum
	print "Fitting Ellipsoids:"
	# Fit ellipsoids
	ellipsoids = get_ellipsoids(lbimg,objnum,minsize=minsize,tol=tol,hull=hull,voxelscale=voxelscale)
	print " Done!"
	print "Calculating Ellipsoid Information:"
	# Calculate information on ellipsoids
	data = elldatacalc(ellipsoids)
	print "Done!"
	endtime = time.time()
	# Output various results of the fitting
	print "Graphing and Writing to files: "
	if outputdata:
		datafileprint(data,outfile)
		summaryfileprint(data,tol,minsize,outfile,threshold,hull,objnum,thresh,endtime-starttime,smooth)
	if display3D or output3D:
		graphrods(data,binimg,outfile,display3D,output3D,voxelscale=voxelscale)
	if displayhist or outputhist:
		if len(data)>1: graphdata(data,outfile,displayhist,outputhist)
	print "Done!"
	return data

'''
FUNC: cubifyvoxels
PURPOSE: Uses given voxel scaling and tricubic interpolation to convert voxel sizes to cubes
INPUT: 	image [np.ndarray]
	voxelscale [list]
OUTPUT:
NOTES:
'''
def cubifyvoxels(image,voxelscale):
	vox = min(voxelscale)
	zoomimg = ndimage.interpolation.zoom(image,[v/vox for v in voxelscale])
	return zoomimg

'''
FUNC: graphdata
PURPOSE: plot various data types calculated from the fitted ellipsoids
INPUT:	data [list] : output of elldatacalc
	outfile [str] : SEE outfile in EllipsoidFit
	display [bool]: SEE display in EllipsoidFit
OUTPUT: saves images containing the plots
NOTES:
'''
def graphdata(data,outfile,displayhist,outputhist):
	# phi (angle in x y plane)
	figphi = plt.figure()
	axphi = figphi.add_subplot(111)
	n,bins,patchs =axphi.hist([d["phi"] for d in data if isinstance(d["phi"],float)],bins=50)
	if outputhist: figphi.savefig(outfile+"_phihist.png")
	if displayhist: figphi.show()
	plt.close(figphi)

	# length
	figlen = plt.figure()
	axlen=figlen.add_subplot(111)
	n,bins,patchs =axlen.hist([d["len"] for d in data if isinstance(d["len"],float)],bins=50,range=(0,50))
	if outputhist: figlen.savefig(outfile+"_lenhist.png")
	if displayhist: figlen.show()
	plt.close(figlen)

	# theta (angle from z-axis)
	figtheta = plt.figure()
	axtheta=figtheta.add_subplot(111)
	n,bins,patchs =axtheta.hist([d["theta"] for d in data if isinstance(d["theta"],float)],bins=50)
	if outputhist: figtheta.savefig(outfile+"_thetahist.png")
	if displayhist: figtheta.show()
	plt.close(figtheta)

	# z-coord
	figz = plt.figure()
	axz=figz.add_subplot(111)
	n,bins,patchs =axz.hist([d["centroid"][0] for d in data if isinstance(d["centroid"][0],float)],bins=50)
	if outputhist: figz.savefig(outfile+"_zhist.png")
	if displayhist: figz.show()
	plt.close(figz)
'''
FUNC: summaryfileprint
PURPOSE: prints a .dat file with average data for the ellipsoids
INPUT:	data [list] : output of elldatacalc
	outfile [str] : SEE outfile in EllipsoidFit
OUTPUT: a .dat file
NOTES:
'''	
def summaryfileprint(data,tol,minsize,outfile,threshold,hull,objnum,thresh,timeel,smooth):
	f = open(outfile+"_metadata.dat","w")
	f.write(outfile+"_metadata.dat ... "+time.asctime(time.localtime(time.time())))
	f.write("\n {0:40s} {1:10s}".format("Time Elapsed: ",str(timeel)))
	f.write("\n {0:40s} {1:10s}".format("Threshold Method: ",str(threshold)))
	f.write("\n {0:40s} {1:10f}".format("Threshold: ",float(thresh)))
	f.write("\n {0:40s} {1:10d}".format("Minimum Object Size: ",minsize))
	f.write("\n {0:40s} {1:10f}".format("MVEE tolerance: ",tol))
	f.write("\n {0:40s} {1:10s}".format("Hull: ",str(hull)))
	f.write("\n {0:40s} {1:10s}".format("Smoothing: ",smooth))
	f.write("\n {0:40s} {1:10d}".format("Object Number: ",objnum))
	f.write("\n {0:40s} {1:10d}".format("Ellipsoid Number: ",len(data)))
	f.write("\n --------------------------------------------------------")
	for i in data[0].keys():
		if not isinstance(data[0][i],np.ndarray):
			f.write("\n {0:40s} {1:10f}".format("Mean "+i+": ",np.mean([d[i] for d in data])))
			f.write("\n {0:40s} {1:10f}".format("Median "+i+": ",np.median([d[i] for d in data])))
			f.write("\n {0:40s} {1:10f}".format("Standard Deviation "+i+": ",np.std([d[i] for d in data])))
			f.write("\n --------------------------------------------------------")
	f.close()
	
'''
FUNC: datafileprint
PURPOSE: prints ellipsoid data to a .csv file
INPUT:	data [list] : output of elldatacalc
	outfile [str] : SEE outfile in EllipsoidFit
OUTPUT: a .csv files containing data for each ellipsoid
NOTES:
'''
def datafileprint(data,outfile):
	# open files
	featdata = open(outfile+"_featuredata.csv","w")

	# Retrieves data types and prints as header
	datakeys = data[0].keys()
	for k in datakeys:
		featdata.write(k+",")
	featdata.write("\n")
	
	# prints the data of each ellipsoid
	for d in data:
		for i in d:
			# Format array outputs if needed
			if isinstance(d[i],np.ndarray):
				dims = d[i].shape
				if len(dims) == 1:
					featdata.write("[")
					for j in d[i]:
						featdata.write(str(j)+" ")
					featdata.write("]")
				elif len(dims) == 2:
					featdata.write("[")
					for j in d[i]:
						featdata.write("[")
						for k in j:
							featdata.write(str(k)+" ")
						featdata.write("]")
					featdata.write("]")
				else: featdata.write("2+ dimen array!")
				featdata.write(",")
			else: featdata.write(str(d[i])+",")
		featdata.write("\n")
	featdata.close()
			
'''
FUNC: graphrods
PURPOSE:	3D Plots and saves/displays the binarized image and the ellipsoids 
		fitted to the objects	
INPUT:	data [list] : output of elldatacalc
	outfile [str] : SEE outfile in EllipsoidFit
	binimg [np.ndarray]: SEE output of thresholding
	display [bool]: SEE display in EllipsoidFit
OUTPUT:	image of the 3D plot
NOTES:
'''
def graphrods(data,binimg,outfile,display3D,output3D,voxelscale=[1,1,1]) :
	# Set up 3D plotting
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	pxinds = np.transpose(np.argwhere(binimg))
	ax.scatter(pxinds[1]*voxelscale[0],pxinds[2]*voxelscale[1],pxinds[0]*voxelscale[2])

	# Plot each ellipsoid
	for i in data:
		Mrot = i["rotation matrix"]
		radii = i["lengths"]
		c = i["centroid"]
		# generate parametric variables
		u = np.linspace(0,2*np.pi,100)
		v = np.linspace(0,np.pi,100)

		# Generate an ellipse centered at the origin with longest axis on
		# the z-axis.
		z = radii[0]*np.outer(np.cos(u),np.sin(v))
		x = radii[1]*np.outer(np.sin(u),np.sin(v))
		y = radii[2]*np.outer(np.ones(np.size(u)),np.cos(v))

		# Transform each point using the rotation matrix and centroid of the
		# ellipsoids.
		for j in range(len(x)):
			for k in range(len(x)):
				[z[j,k],x[j,k],y[j,k]] = np.dot([z[j,k],x[j,k],y[j,k]],Mrot)+c

		maxaxis = max(binimg.shape)*voxelscale[np.argmax(binimg.shape)]
		ax.plot_wireframe(x, y, z,  rstride=4, cstride=4,color="red")
		ax.set_xlim(0,maxaxis)
		ax.set_ylim(0,maxaxis)
		ax.set_zlim(0,maxaxis)
	
	if output3D: fig.savefig(outfile+"_ellipsoidfits.png")
	if display3D:
		fig.show()
	#plt.close(fig)

'''
FUNC: get_ellipsoids
PURPOSE: Takes the binarized and labelled 3D image, isolates the objects of interest
	and sends them into MVEE for fitting.
INPUT:	lbimg [np.ndarray] - 3D image with contiguous patches of pixels labelled
	objnum [int] - number of objects
	minsize [int] - minimum size threshold to fit an object
	tol [float] - the criterion used to determine when to end the MVEE fitting
	hull [bool] - SEE: hull in EllipseFit
OUTPUT: ellipsoids [list] - a list where each list element has a dxd array 
	containing the characteristic matrix of the ellipsoid and the centroid of 
	the ellipsoid.
NOTES:
'''
def get_ellipsoids(lbimg,objnum,minsize = 10,tol = 0.01,hull=True,voxelscale=[1,1,1]):
	ishape = lbimg.shape
	ellipsoids = [] # Stores ellipsoid Matrices and centroids
	goodobj = [] # Store indices of fitted groups
	for i in range(1,objnum+1):
		# Counter for the user
		sys.stdout.write("\r {0} / {1}".format(i,objnum))
		sys.stdout.flush()

		# Create binary image of the object in question only
		objarr = np.zeros(ishape)
		objarr[lbimg == i] = 1
		objindices = np.argwhere(objarr)
		# Check against the minimum size
		if objindices.shape[0] > minsize:
			if hull:
				
				# Attempt to get convex hull for the object
				try:
					hull = spatial.ConvexHull(objindices)
					points= np.transpose(objindices[hull.vertices])
					tempE = MVEE.MVEE(points,tol = tol)
				# If this fails, get the ellipsoid of the object itself
				except:
					print "Convex Hull error"
					tempE = MVEE.MVEE(np.transpose(objindices),tol=tol)
				# MVEE returns None if there is an error
				# Saves the matrices and centroid of the fitted ellipsoid
				if tempE:
					goodobj.append(i)
					tempE.append(objindices.shape[0]*np.cumprod(voxelscale)[2])
					ellipsoids.append(tempE)
			else: 
				tempE = MVEE.MVEE(np.transpose(objindices),tol=tol)
				if tempE:
					goodobj.append(i)
					tempE.append(objindices.shape[0]*np.cumprod(voxelscale)[2])
					ellipsoids.append(tempE)
	return ellipsoids
	
# --------------------------------------------------------------------
'''
FUNC: thresholding
PURPOSE: calculates a suitable threshold and returns a binary 3D image
INPUT:	image [np.ndarray] - the image of interest
	threshold [str or int] - SEE: threshold in EllipseFit
OUTPUT: binimg - binary thresholded image
	threshold - calculated threshold
NOTES: 	Currently uses threshold finding methods available through skimage.
	May add more methods later.
'''
def thresholding(image, threshold):
	if isinstance(threshold,int) or isinstance(threshold,float):
		thresh = threshold
	elif isinstance(threshold,str):
		# Assume its Ok to use the same threshold for each layer.
		parsestr = threshold.split(" ")
		parsestr = [i.lower() for i in parsestr]
		parsestr = set(parsestr)
		if "otsu" in parsestr:
			if "global" in parsestr:
				ind = np.argmax([np.mean(i) for i in image])
				thresh = filters.threshold_otsu(image[ind])
				print thresh, ind
			elif "local" in parsestr:
				radius = int(parsestr[2])
				mask = morphology.disk(radius)
				thresh = filters.rank.otsu(image,mask)
		if "li" in parsestr:
			thresh = filters.threshold_li(image)
		if "yen" in parsestr:
			thresh = filters.threshold_yen(image)

	threshinds = image<thresh
	binimg =  np.ones(image.shape)
	binimg[threshinds] = 0
	return binimg, thresh

# --------------------------------------------------------------------
'''
FUNC: elldatacalc
PURPOSE: calculates a variety of information on the fitted ellipsoids
INPUT: ellipsoids [list] - SEE: ellipsoids output in get_ellipsoids
OUTPUT: data [list] - a list of dictionaries with the calculated values
NOTES:
'''
# Calculating information on the ellipsoids
def elldatacalc(ellipsoids):
	data = []
	for i in ellipsoids:
		datadict = {}
		U,Q,V = np.linalg.svd(i[0])
		lengths = np.array([m.sqrt(1./q) for q in Q])
		eivals,eivecs = np.linalg.eig(i[0])
		theta,phi = spherical_coords(eivecs[:,np.argmin(eivals)])

		datadict["theta"] = theta
		datadict["phi"] = phi
		datadict["ellipsoid matrix"] = i[0]
		datadict["eigenvectors"] = np.transpose(eivecs)
		datadict["eigenvalues"] = eivals
		datadict["centroid"] = i[1]
		datadict["rotation matrix"] = V
		datadict["principal numbers"] = Q
		datadict["phi"] = phi
		datadict["theta"] = theta
		datadict["lengths"] = lengths
		datadict["len"] = max(lengths)
		datadict["ellipsoid vol"] = 4*m.pi*lengths[0]*lengths[1]*lengths[2]/3
		datadict["object vol"] = i[2]
		datadict["ellipsoidiness"] = datadict["object vol"]/datadict["ellipsoid vol"]
		datadict["aspect ratio"] = lengths[2]/lengths[1]
		datadict["width"] = lengths[1]
		data.append(datadict)
	return data

# --------------------------------------------------------------------
'''
FUNC: spherical_coords
PURPOSE: calculated angles in spherical coords of a vector
INPUT: vec [list] - vector 
OUTPUT: theta, phi - angles in spherical. Phi is the angle in the xy-plane,
	theta is the angle from the z-axis.
NOTES: The vector is in the form [z,x,y]
'''
def spherical_coords(vec):
	theta = m.acos(vec[0]/m.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]))
	if vec[1] ==0:
		phi = m.pi/2
	else:
		phi = m.atan(vec[2]/vec[1])
	return theta, phi
	


