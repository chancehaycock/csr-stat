# File to contain all utility functionality required for work with point processes
# in R^2.
#
# Enables user to produce and visualise functional and real-valued summary statistics
#
# =============================================================================

# ===============
#   Formatting
# ===============
# TODO
# - Add Type Hinting
# - Google Docstring format: see example_google_docstring.py

# ===============
#  Testing
# ===============
# - Write tests for functions via PyTest

# =============================================================================
# =============================================================================
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sbn
from typing import Any, List, Type


class Point:
	"""Simple base class for a point in R^2."""

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

	"""Class to represent a general point process in R^2"""

	def __init__(self, points: np.ndarray, window:dict=None, boundary=None):

		# Store both array of Points, and arrays of co-ordinates
		# for convenience
		self.points = points
		self.x = np.array([p.x for p in points])
		self.y = np.array([p.y for p in points])
		self.marks = np.array([p.mark for p in points])

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

	def display(self):
		sbn.scatterplot(self.x, self.y, alpha=0.7, s=15.0)
		plt.xlabel('x')
		plt.ylabel('y')
		plt.show()

	def intensity_map(self, bandwidth: float = 0.1):
		return
	
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
	points=[]

	for _ in range(num_events):
		point=Point(np.random.uniform(window["x_min"], window["x_max"]),\
		            np.random.uniform(window["y_min"], window["y_max"]))
		points=np.append(points, point)
	return Point_Process(points, window)

def thinningpp(PP):
	return

# Generate an inhomogenoues point process, option for marks. Use thinning.
def inhom_poisspp(intensity_function, window: dict, factor_marks: bool = False):
	return

# Generate Point process from data at filename
def scanpp(filename: str, window: dict, marks: bool=False):
	return

def pandas2points(df: pd.DataFrame, x_label: str, y_label: str, marks_label: str=None)-> np.ndarray:

	"""Utility function to convert a dataframe of x, y columns into
	an array of Point objects.

	Args:
		df (pandas Dataframe): Input data
		x_label (str): Name of the column in df associated to x co-ords
		y_label (str): Name of the column in df associated to y co-ords
		marks_label (str): [Inactive] Name of the column in df associated
		                   to the mark values of a (x, y) co-ord

	Returns:
		Array of Point objects 
	"""

	df = df[x_label, y_label, marks_label]
	df = df.drop_na()
	x = df[x_label].to_numpy()
	y = df[y_label].to_numpy()

	if len(x) != len(y):
		print('x and y arrays must have same length.')
		return

	points = []
	for i in range(len(x)):
		points.append(Point(x[i], y[i]))

	if marks_label is not None:
		marks = df[marks_label].to_numpy()
		return points, marks

	return points


