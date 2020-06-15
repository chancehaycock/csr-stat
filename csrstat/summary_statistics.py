import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import seaborn as sbn
from csrstat.point_process_utils import Point_Process, Point, hom_poisspp
from typing import Any, List, Type



def min_neighbour_dists(PP: Type[Point_Process]) -> np.ndarray :

	"""Calculates the nearest neighbour distance to each of the N samples PP. 
	Returns nearest neighbour distances as a np.array of floats of length N.
	Args:
		PP: a point process object
	Returns: 
		np array"""

	minimum_distances = []
	for i in range(0, PP.num_points):
		d_min=float('inf')
		for j in range(0, PP.num_points):
			d=np.sqrt((PP.x[j] - PP.x[i])**2 + (PP.y[j] - PP.y[i])**2)
			if ((d<d_min) and (i!=j)):
				d_min=d
                
		minimum_distances=np.append(minimum_distances, d_min)
        
	return minimum_distances

def min_empty_space_dists(PP: Type[Point_Process], num_sampled_points:int)-> np.ndarray:

	"""Calculates the nearest neighbour distance to each of the N samples PP. 
	Returns nearest neighbour distances as a np.array of floats of length N.
	Args:
		PP: a point process object
	Returns: 
		np array"""

	empty_space_distances = []
	for i in range(0, num_sampled_points):
		x = np.random.uniform(PP.window["x_min"], PP.window["x_max"])
		y = np.random.uniform(PP.window["y_min"], PP.window["y_max"])
		d_min=float('inf')
		for j in range(0, PP.num_points):
			d=np.sqrt((PP.x[j] - x)**2 + (PP.y[j] - y)**2)
			if ((d<d_min) and (i!=j)):
				d_min=d
                
		empty_space_distances=np.append(empty_space_distances, d_min)
        
	return empty_space_distances



def G_hat(PP: Type[Point_Process], r:np.ndarray, plot:bool=False)-> np.ndarray:

	"""Calculatesan estimate of the G(r) over a range of distances.
     
	Args:
		PP: a point process object
		r: an np array of input values
        plot: default False, True to see plots of G(r) 
	Returns: 
		G(r) as np.ndarray of floats"""

	minimum_distances=min_neighbour_dists(PP)
	intensity=PP.homo_intensity_estimate()
	expected = 1 - np.exp(-intensity*np.pi*r**2)
	total_num_points = PP.num_points
	G_hat = []
	for r_in in r:
		count = 0
		for dist in minimum_distances:
			if dist <= r_in:
				count += 1
		G_hat.append(count/total_num_points)
    

	if(plot):
		sbn.scatterplot(r, G_hat, alpha=0.75, label='Observed G')
		plt.plot(r, expected, c='k', label='Expected G', alpha=0.6)
		plt.legend()
    
		plt.xlabel('r')
		plt.ylabel('G(r)')

	return np.array(G_hat)



def F_hat(PP: Type[Point_Process], r:np.ndarray, num_sampled_points:int = 1000, plot:bool =False) -> np.ndarray:

	"""Calculates and plots an estimate of the F(r) over a range of distances.
     
	Args:
		PP: a point process object
		r: an np array of input values
        num_sampled_points: default 1000, number of random points used to determine F
        plot: default False, True to see plots of F(r)
	Returns: 
		F(r) as np.ndarray of floats"""

	minimum_distances = min_empty_space_dists(PP, num_sampled_points)
	intensity = PP.homo_intensity_estimate()
	expected = 1 - np.exp(-intensity*np.pi*r**2)
	F_hat = []
	for r_in in r:
		count = 0
		for dist in minimum_distances:
			if dist <= r_in:
				count += 1
		F_hat.append(count/num_sampled_points)
   
	if(plot):
		sbn.scatterplot(r, F_hat, alpha=0.75, label='Observed F')
		plt.plot(r, expected, c='k', label='Expected F', alpha=0.6)
		plt.legend()
    
		plt.xlabel('r')
		plt.ylabel('F(r)')

	return np.array(F_hat)

def K_hat(PP:Type[Point_Process], r:np.ndarray, restrict_domain:bool = True, plot:bool = False) -> np.ndarray:

	"""Calculates and plots an estimate of the K(r) over a range of distances.

	Args:
		PP: a point process object
		r: an np array of input values
		restrict_domain: default True, restricting the domain to counter edge effect
		plot: default False, True to see plots of K(r)
	Returns: 
		K(r) as np.ndarray of floats"""

	intensity=PP.homo_intensity_estimate()
	expected = np.pi*r**2
	K_hat=[]

	distance_matrix=[]
	for i in range(0, PP.num_points):
		d=[]
		for j in range(0, PP.num_points):
			d=np.append(d, np.sqrt((PP.x[j]-PP.x[i])**2+(PP.y[j]-PP.y[i])**2))
			
		distance_matrix.append(d)

	for r_in in r:
		count=0
		n=0
		for p in range(0, PP.num_points):
			
			if(restrict_domain):
				cond = (( PP.x[p]>r_in+PP.window["x_min"]) and (PP.x[p]<PP.window["x_max"]-r_in) 
				and (PP.y[p]>r_in+PP.window["y_min"]) and (PP.y[p]<PP.window["y_max"]-r_in))

			else:
				cond=True
			
			if(cond):
				n+=1
				for k in range(0, PP.num_points):
					
					if(distance_matrix[p][k]<=r_in and distance_matrix[p][k]>0):
						count+=1            
		K_hat.append(count/(n*intensity))

	if(plot):
		sbn.scatterplot(r, K_hat, alpha=0.75, label='Observed K')
		plt.plot(r, expected, c='k', label='Expected K', alpha=0.6)
		plt.legend()

		plt.xlabel('r')
		plt.ylabel('K(r)')

	return np.array(K_hat)


def L_hat(PP: Type[Point_Process], r:np.ndarray, restrict_domain:bool = True, plot:bool = False)-> np.ndarray:

	"""Calculates and plots an estimate of the L(r) over a range of distances.
	Note here, under CSR, L is the zero function.
     
	Args:
		PP: a point process object
		r: an np array of input values
        restrict_domain: default True, restricting the domain to counter edge effect
        plot: default False, True to see plots of L(r)
	Returns: 
		L(r) as np.ndarray of floats"""

	# Call Ripley K
	K = K_hat(PP, r, restrict_domain=restrict_domain, plot=False)

	# Normalise as per definition of L
	L_hat = np.sqrt(K / np.pi) - r

	# Option to Plot
	if(plot):
		sbn.scatterplot(r, L_hat, alpha=0.75, label='Observed L')
		plt.plot(r, np.zeros(len(r)), c='k', label='Expected L', alpha=0.6)
		plt.legend()

		plt.xlabel('r')
		plt.ylabel('L(r)')
	return np.array(L_hat)

def O_hat(PP: Type[Point_Process], r:np.ndarray, bandwidth:float=0.1, restrict_domain:bool=True, kernel:str="BK", plot:bool=False)-> np.ndarray:

	"""Calculates and plots an estimate of the O(r) over a range of distances.
		Options for kernel and edge effects.

	Args:
		PP: a point process object
		r: an np array of input values 
		bandwidth : float determining the bandwidth of ring
		restrict_domain: if False edge effects will be ignored, default True
		kernel: kernel for O_hat. Default is "BK" for box kernel, which is a 
		flat count. Option of "EK" for Epanechnikov kernel, a quadratic count 
		which favour points closer to center of ring.
	Returns: 
		O(r) as np.ndarray of floats"""

	intensity=PP.homo_intensity_estimate()
	O_hat=[]

	distance_matrix=[]
	for i in range(0, PP.num_points):
		d=[]
		for j in range(0, PP.num_points):
			d=np.append(d, np.sqrt((PP.x[j]-PP.x[i])**2+(PP.y[j]-PP.y[i])**2))

		distance_matrix.append(d)

	for r_in in r:
		if r_in < bandwidth: #r values below bandwidth will be excluded as O_ring not well defined
				O_hat.append(np.nan)
				continue
		count=0
		n=0
		for p in range(0, PP.num_points):

			if(restrict_domain):
				cond = (( PP.x[p]>r_in+PP.window["x_min"]) and (PP.x[p]<PP.window["x_max"]-r_in) 
				and (PP.y[p]>r_in+PP.window["y_min"]) and (PP.y[p]<PP.window["y_max"]-r_in))

			else:
				cond=True

			if(cond):
				n+=1
				for k in range(0, PP.num_points):

					if(distance_matrix[p][k]>=r_in-bandwidth and distance_matrix[p][k]<=r_in+bandwidth 
						and distance_matrix[p][k]>0):
						if(kernel=="BK"):
							count+=(1/(2*bandwidth))
						if(kernel=="EK"):
							count+=(3/(4*bandwidth))*(1-((distance_matrix[p][k]-r_in)**2/(bandwidth**2)))
			
		O_hat.append(count/(2*np.pi*r_in*n))


	if(plot):

		sbn.scatterplot(r, O_hat, alpha=0.75, label='Observed O(r)')
		plt.plot(r, intensity*np.ones(len(r)), c='k', label='Expected O(r)', alpha=0.6)
		plt.legend()

		plt.xlabel('r')
		plt.ylabel('O(r)')
	return np.array(O_hat)

def PC_hat(PP: Type[Point_Process], r:np.ndarray, bandwidth:float=0.1, restrict_domain:bool=True, kernel:str="BK", plot:bool=False)-> np.ndarray:

	"""Calculates and plots an estimate of the PCF(r) over a range of distances.
	Options for kernels include  and edge effects
		
	Args:
		PP: a point process object
		r: an np array of input values
		bandwidth : float determining the bandwidth of ring
		restrict_domain: if False edge effects will be ignored, default True
		kernel: kernel for O_hat. Default is "BK" for box kernel, which is a 
		flat count. Option of "EK" for Epanechnikov kernel, a quadratic count 
		which favour points closer to center of ring.
	Returns: 
		PCF(r) as np.ndarray of floats"""

	intensity=PP.homo_intensity_estimate()
	O=O_hat(PP, r, bandwidth=bandwidth, restrict_domain=restrict_domain, kernel=kernel)
	PC_hat=O/intensity

	if(plot):

		sbn.scatterplot(r, PC_hat, alpha=0.75, label='Observed PC(r)')
		plt.plot(r, np.ones(len(r)), c='k', label='Expected PC(r)', alpha=0.6)
		plt.legend()

		plt.xlabel('r')
		plt.ylabel('PC(r)')
	return np.array(PC_hat)

def J_hat(PP: Type[Point_Process], r:np.ndarray, num_F_sampled_points:int=1000, plot:bool=False)-> np.ndarray:

	"""Calculates and plots an estimate of the J(r) over a range of distances
	such that F(r) < 1. 

	                   J(r) := (1 - G(r)) / (1 - F(r))

	Returns J(r) as a np.ndarray of length len(r). As the function is only defined
	for values of r such that F(r) < 1, the returned array consists of

	         [ {J(r) : r satisfies F(r) < 1} + {array of np.nans} ]

	Note here, under CSR, J is the constant function of value one. As the function 
	is comprised of a quotient of F and G, any edge effects should cancel out.
     
	Args:
		PP: a point process object
		r: an np array of input values
		num_F_sampled_points: number of sampled points in the domain used to 
		                      estimate the empty space function, F.
        plot: default False, True to see plots of J(r)
	Returns: 
		J(r) as np.ndarray"""

	F = F_hat(PP, r, num_sampled_points=num_F_sampled_points, plot=False)
	G = G_hat(PP, r, plot=False)

	# Max value of F function
	F_max = np.max(F)
	# If F(r) != 1 anyway, proceed as normal
	if (F_max < 1):
		J_hat = (1 - G) / (1 - F)
		max_index = len(r) - 1
	# Otherwise, only calculate values for J_hat where F_hat < 1
	else:
		max_index = np.min(np.where(F == 1.0)) - 1
		J_hat = (1 - G[:max_index]) / (1 - F[:max_index])

	# Create array of nans to be concatenated with J_hat to ensure
	# the returned J_hat array has the same length as input array r.
	nan_buffer = np.full(len(r) - len(J_hat), np.nan)

	# Option to Plot
	if(plot):
		fig, (ax1, ax2 , ax3) = plt.subplots(3, 1, figsize=(8,6))
		sbn.scatterplot(r, F, alpha=0.75, label='Observed F', ax=ax1)
		sbn.scatterplot(r, G, alpha=0.75, label='Observed G', ax=ax2)
		
		sbn.lineplot(r[:max_index], J_hat[:max_index], alpha=0.75,\
                     label='Observed J', linewidth=1.0, ax=ax3)
		plt.plot(r, np.ones(len(r)), c='k', label='Expected J', alpha=0.6)
		plt.legend()

		plt.xlabel('r')
		ax1.set_ylabel('F(r)')
		ax2.set_ylabel('G(r)')
		ax3.set_ylabel('J(r)')
		fig.tight_layout()
		plt.show()

	return np.concatenate((J_hat, nan_buffer))


def mean_min_neighbour_dist(PP: Type[Point_Process])-> float:

	"""Simple utility function to calculate the mean minimum neighbour distance
	of a point process. Used in CSR hypothesis testing.

	Args: 
		PP (Point_Process): Point_Process to calculate the statistic on.

	Returns:
		Average Minimum Neighbour Distance (float)
	"""

	minimum_distances = min_neighbour_dists(PP)

	return np.mean(minimum_distances)


def CSR_hypothesis_test_dmin(observed_PP, significance_level:float, n_sim:int=500,
                             plot_dist:bool=False)->None:

	"""Assumes CSR as the null hypothesis and tests whether the observed Point
	Process is significantly different from the expected CSR case using the
	mean minimum nearest neighbour distance statistic. Performs a two tail
	test to include cases of clustering and dispersion. 

	Args:
		observed_PP (Point_Process): A realisation of a Point Process object to
		                             be tested.
		significance level (float): Statistical significance level
		n_sim (int): Number of simulations performed to estimate the
		             distribution.
		plot_dist (bool): Option to display the simulated distribution along
		                  with the observed test statistic.

		"""

	window = observed_PP.window

	# Test statistics from the observed point process
	observed_min_dist_mean = mean_min_neighbour_dist(observed_PP)

	# Assuming homogeneous
	intensity = observed_PP.homo_intensity_estimate()

	d_min_means = []
	for i in range(n_sim):

		# Simulate a homogeneous poisson process on the same window
		# as the observed
		sim_PP = hom_poisspp(intensity, window)
		d_min_means.append(mean_min_neighbour_dist(sim_PP))

		# Progress bar for user. Simulation can be computationally intensive.
		if i%50 == 0:
			print('{}% of simulations complete.'.format(i*100/n_sim))

	if (plot_dist):
		sbn.distplot(d_min_means, label='Estimated Distribution')
		plt.axvline(x=observed_min_dist_mean, c='k', label='observed')
		plt.xlabel('Mean Minimum Neighbour Distance')
		plt.ylabel('Count')
		plt.title('Sampled Distribution using {} samples'.format(n_sim))
		plt.legend()
		plt.show()

	# Estimate p-value by calculating the proportion of simulations that
	# have mean min neighbour distance GREATER than the observed.
	tail_count = 0
	for mean in d_min_means:
		if mean >= observed_min_dist_mean:
			tail_count += 1

	# Times by two to get two sided test.
	p_value = 2 * np.minimum(n_sim + 1 - tail_count, tail_count + 1) / (n_sim + 1)

	# Print summary to the user.
	print('p-value obtained from 2-sided Test: {0:.3f}'.format(p_value))
	if p_value < significance_level:
		print('p-value is less than {} and so the CSR hypothesis'\
              .format(significance_level), 'can be rejected.')
	else :
		print('p-value is greater than or equal to {} and so the CSR'\
              .format(significance_level), 'hypothesis can not be rejected.')

	return

def simulate_summary_statistic(PP: Type[Point_Process], r:np.ndarray, summary_func:object, n_sims:int=100, plot:bool=False)->tuple:

	"""Compares metrics between a given point process, and an average of n_sims homogenous 
	poisson point processes. Returns the max and mean difference between the two. 

	Args:
		PP: a point process object
		r: an np array of input values
		summary_func: a defined point process functions
        plot: default False, True to see plots of metric compared to metric for poisson with
			min and max envelopes for n_sims iterations
	Returns: 
		tuple of (mean difference, max difference)"""



	func_stack=[]
	intensity=PP.homo_intensity_estimate()
	for i in range(0, n_sims):
		sim=hom_poisspp(intensity, PP.window)
		func_stack.append(summary_func(sim, r))
		b = "Loading{" + "#"*round(i*10/n_sims)+'{}%'.format(round(i*10/n_sims)*10) + "}"
		print (b, end="\r")

	func_mean = np.mean(func_stack, axis=0)
	func_max = np.max(func_stack, axis=0)
	func_min = np.min(func_stack, axis=0)
	observed = summary_func(PP, r)
    
	if(plot):
		plt.plot(r, func_mean, c='k', label='Mean ' + summary_func.__name__[0] + '(r) CSR Simulation', alpha=0.9)
		sbn.scatterplot(r, observed, label='Observed '+ summary_func.__name__[0] + '(r)', s=25.0)
		plt.fill_between(r, func_min, func_max, alpha=0.3, label='Max/Min Envelope')
		plt.xlabel('r')
		plt.ylabel(summary_func.__name__[0] + '(r)')
		plt.legend()
		plt.show()
    
	max_diff=0
	mean_diff=0
	for i in range(0, len(r)):
		diff=abs(observed[i]-func_mean[i])
		if(diff>max_diff):
			max_diff=diff

		mean_diff+=diff
		
	mean_diff=mean_diff/len(r)
		
	return(mean_diff, max_diff)
