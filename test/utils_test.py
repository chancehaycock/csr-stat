# Contains all tests for functions declared in src/point_process_utils.py
# ======================================================================
from csrstat.utils import Point, Point_Process
import pytest

# Point Class Tests
@pytest.fixture
def point():
	p = Point(3, 4)
	return p

def test_point1(point):
	assert point.x == 3 and point.y == 4

def test_point2(point):
	assert point.magnitude() == 5

def test_point3(point):
	assert point.mark == 1

def test_point4(point):
	point.update_mark(5)
	assert point.mark == 5

def test_point5(point):
	point.update_mark(4)
	new_point = point.translate(1, 2)
	assert new_point.x == 4 and new_point.y == 6 and new_point.mark == point.mark

# Point Process Class Tests
@pytest.fixture
def pp():
	p1 = Point(0, 0)
	p2 = Point(1, 2)
	p3 = Point(4, 3)
	return Point_Process([p1, p2, p3])

def test_pp1(pp):
	x = pp.x
	y = pp.y
	marks = pp.marks
	x_expected = [0, 1, 4]
	y_expected = [0, 2, 3]
	for i in range(len(x)):
		if x_expected[i] != x[i] or y_expected[i] != y[i] or marks[i] != 1:
			assert False
	assert True

def test_pp2(pp):
	assert pp.num_points == 3

def test_pp3(pp):
	assert pp.area() == 12

def test_pp4(pp):
	cond = pp.window['x_min'] == 0 and\
	       pp.window['x_max'] == 4 and\
	       pp.window['y_min'] == 0 and\
	       pp.window['y_max'] == 3
	assert cond

# Example G Summary Statistic Test
from pointpats import PointPattern, G
import numpy as np
from csrstat.summary import G_hat

@pytest.fixture
def pysal_pp():
	tuple_points = np.array([[66.22, 32.54], [22.52, 22.39], [31.01, 81.21],
    	                     [9.47, 31.02],  [30.78, 60.10], [75.21, 58.93],
    	                     [79.26,  7.68], [8.23, 39.93],  [98.73, 77.17],
    	                     [89.78, 42.53], [65.19, 92.08], [54.46, 8.48]])
	return PointPattern(tuple_points)

@pytest.fixture
def native_pp():
	tuple_points = np.array([[66.22, 32.54], [22.52, 22.39], [31.01, 81.21],
    	                     [9.47, 31.02],  [30.78, 60.10], [75.21, 58.93],
    	                     [79.26,  7.68], [8.23, 39.93],  [98.73, 77.17],
    	                     [89.78, 42.53], [65.19, 92.08], [54.46, 8.48]])
	points = []
	for point in tuple_points:
	    points.append(Point(point[0], point[1]))
	return Point_Process(points)


def test_G(pysal_pp, native_pp):
	# Currently Depreciation Warning

	G_pysal = G(pysal_pp, intervals=100)
	domain = G_pysal.d
	G_native = G_hat(native_pp, domain, plot=False)

	pysal_array = G_pysal.G
	native_array = G_native
	print(pysal_array)
	print(native_array)

	assert np.array_equal(pysal_array, native_array)
