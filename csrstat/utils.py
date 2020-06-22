# File to contain all utility functionality required for work with point processes
# in R^2.
#
# Enables user to produce and visualise functional and real-valued summary statistics
#
# =============================================================================

from __future__ import annotations
from typing import Any, List, Type
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sbn
import utm

class Point:
	"""Simple base class for a point in R^2. Can also be assigned a mark/weight.
	Useful for analysis of marked point processes."""

	def __init__(self, x_coord: float, y_coord: float, mark:Any=1):
		self.x = x_coord
		self.y = y_coord
		self.mark = mark

	def magnitude(self)-> float:
		return np.sqrt(self.x**2 + self.y**2)

	def translate(self, s: float, t: float)-> Type[Point]:
		return Point(self.x + s, self.y + t, self.mark)

	def scale(self, lmbda: float)-> Type[Point]:
		return Point(self.x * lmbda, self.y * lmbda, self.mark)
	
	def update_mark(self, new_mark: Any):
		self.mark = new_mark

	def __str__(self)->str:
		return "({}, {})".format(self.x, self.y)


class Point_Process:

	"""Class to represent a general point process in R^2."""

	def __init__(self, points: np.ndarray, window:dict=None, boundary=None):

		# Store both array of Points, and arrays of co-ordinates
		# for convenience
		self.points = points
		self.x = np.array([p.x for p in points])
		self.y = np.array([p.y for p in points])
		self.marks = np.array([p.mark for p in points])
		self.distance_matrix = None

		# If no window argument passed, 'guess' it from the data
		if (window is None):
			self.window = {'x_min': self.x.min(), 'x_max': self.x.max(),
			               'y_min': self.y.min(), 'y_max': self.y.max()}
		else:
			self.window = window

		self.num_points = len(points)

		# TODO - add functionality for non-rectangular windows
		# What format will this be?
		self.boundary = None

	def area(self)->float:
		return (self.window['x_max'] - self.window['x_min'])\
		     * (self.window['y_max'] - self.window['y_min'])
	
	def distance_matrix_init(self):
		d_matrix=[]
		for i in range(0, self.num_points):
			d=[]
			for j in range(0, i):
				d=np.append(d, np.sqrt((self.x[j]-self.x[i])**2+(self.y[j]-self.y[i])**2))
				
			d_matrix.append(d)

		self.distance_matrix = d_matrix

	def display(self, marker_size: float = 10.0):
		sbn.scatterplot(self.x, self.y, size=self.marks, alpha=0.7, s=marker_size)
		plt.xlabel('x')
		plt.ylabel('y')
		plt.show()

	def intensity_map(self, partitions: int = 100, bandwidth: float = 1, display = True)-> np.ndarray:
		inc_x=np.linspace(self.window["x_min"], self.window["x_max"], partitions)
		inc_y=np.linspace(self.window["y_min"], self.window["y_max"], partitions)
		gridspace_area=(self.area()/(partitions*partitions))
		x=[]
		y=[]
		for i in range(0, self.num_points):
			for _ in range(0, self.marks[i]):
				x.append(self.x[i])
				y.append(self.y[i])

		smooth=sp.ndimage.gaussian_filter(np.flip(np.histogram2d(y, x,
		                                  [inc_y, inc_x])[0], 0)/gridspace_area,
		                                  bandwidth)

		if(display):
			sbn.heatmap(smooth, xticklabels=False, yticklabels=False)

		return smooth

	def possion_num_mark(self, mean_num):

		for i in range(0,self.num_points):
			self.marks[i]=np.random.poisson(mean_num)


	
	# Return summary stats for PP
	def summary(self):
		return

	# =======================
	# Statistical Data Tests
	# =======================

	def homogeneity_test(self):
		# Chi-Squared Quadrat test?
		return

	def isotropy_test(self):
		# test or assume?
		return

	def SOIRS_test(self):
		# If SOIRS satisified, use inhomo estimates of K and PC as in BMW-2020
		return


	# ====================
	# Intensity Estimates
	# ====================

	# Simple Homogeneous case
	def homo_intensity_estimate(self)-> float:
		return self.num_points/self.area()

	# Can be useful for visualising the intensity function. Relies heavily on
	# the choice of bandwith.
	def kernel_smooth_intensity(self, bandwidth: float):
		return

	# Fit an Inhomo PP to the data. This is advised to test if the source of
	# inhomogeneity is understood
	def fit_inhom_pp(self, parameters):
		return


# ==========================================
#          Utility Functions
# ==========================================

# Generate a homogenoues point process, option for marks
def hom_poisspp(intensity: float, window: dict, factor_marks: bool = False)->Type[Point_Process]:
	mean_num_events = intensity * (window["x_max"] - window["x_min"]) * \
	                              (window["y_max"] - window["y_min"])
	num_events = np.random.poisson(mean_num_events)

	# Always generates a PP with at least one event. Should the homo PP
	# simulation always simulate the same number of points?
	if num_events == 0:
		num_events += 1

	points=[]

	for _ in range(num_events):
		point=Point(np.random.uniform(window["x_min"], window["x_max"]),\
		            np.random.uniform(window["y_min"], window["y_max"]))
		points=np.append(points, point)
	return Point_Process(points, window)

def thinningpp(PP):
	return

# Generate an inhomogenous point process, option for marks. Use thinning based on accept reject
def inhom_poisspp(intensity_surface: np.ndarray, window: dict, factor_marks: bool = False)->Type[Point_Process]:
	points=[]

	#create a probability surface by dividing all points by the maximum of 
	# intensity sufaces
	P=intensity_surface/intensity_surface.max()
	#create a homogenous point process intensity equal to maximum from intensity
	# surface
	homPP = hom_poisspp(intensity_surface.max(), window) 

	#loop through all points in the point proces
	for p in range(0, homPP.num_points): 
		#convert x & y to grid locations of the intensity or probability surface
		gridx=int(homPP.x[p]*(window["x_min"]+(len(P)/(window["x_max"]-window["x_min"])))) 
		gridy=int(homPP.y[p]*(window["y_max"]-(len(P)/(window["y_max"]-window["y_min"])))) 
		
		#generate random number, if below probability then accept- otherwise reject
		bool_accept = np.random.uniform(0,1)<=P[gridy][gridx] 
		if(bool_accept):
			points=np.append(points, Point(homPP.x[p], homPP.y[p]))
		
	return(Point_Process(points, window))

# Generate Point process from data at filename
def scanpp(filename: str, window: dict, marks: bool=False):
	return

def pandas2points(df: pd.DataFrame, x_label: str, y_label: str, 
                  marks_label: str=None)-> np.ndarray:

	"""Utility function to convert a dataframe of x, y columns into
	an array of Point objects.

	Args:
		df (pandas Dataframe): Input data
		x_label (str): Name of the column in df associated to x co-ords
		y_label (str): Name of the column in df associated to y co-ords
		marks_label (str): Name of the column in df associated
		                   to the mark values of a (x, y) co-ord

	Returns:
		Array of Point objects 
	"""
	df = df.dropna()
	x = df[x_label].to_numpy()
	y = df[y_label].to_numpy()

	if len(x) != len(y):
		print('x and y arrays must have same length.')
		return

	points = []
	if marks_label is not None:
		marks = df[marks_label].to_numpy()
		for i in range(len(x)):
			points.append(Point(x[i], y[i], marks[i]))
	else:
		for i in range(len(x)):
			points.append(Point(x[i], y[i]))

	return points

def convert_long_lat_2_easting_northing(df:pd.DataFrame) -> pd.DataFrame:

	"""Utility function to convert a SCOOT dataframe of long, lat columns into
	one with easting northing co-ordinates

	Args:
		df (pandas Dataframe): Input data  with columns ['lat', 'lon']

	Returns:
		pd.DataFrame with converted co-ordinates. 
	"""

	# Check if lat and lon in df
	assert set(['lon','lat']) <= set(df.columns)

	lambdafunc = lambda x: pd.Series(utm.from_latlon(x['lat'], x['lon'])[:2])
	df[['x', 'y']] = df.apply(lambdafunc, axis=1)
	return df
