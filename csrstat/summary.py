"""
Module containing summary statistics for a 2D point process.
"""
from typing import Type
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sbn
from csrstat.utils import Point_Process, hom_poisspp
sbn.set()

def min_neighbour_dists(PP: Type[Point_Process]) -> np.ndarray:

    """Calculates the nearest neighbour distance of each of the N events of
    the point process PP. Returns nearest neighbour distances as a np.array
    of floats of length N.
    Args:
        PP: a point process object
    Returns:
        np array"""

    minimum_distances = []
    for i in range(0, PP.num_points):
        d_min = float('inf')
        for j in range(0, PP.num_points):
            d = np.sqrt((PP.x[j] - PP.x[i])**2 + (PP.y[j] - PP.y[i])**2)
            if d < d_min and i != j:
                d_min = d

        minimum_distances.append(d_min)

    return np.array(minimum_distances)


def min_empty_space_dists(PP: Type[Point_Process], num_sampled_points: int)-> np.ndarray:

    """Calculates the nearest neighbour distance to each of the `num_sampled_points`
    points in the domain of the point process, PP. Returns nearest neighbour distances
    as a np.array of floats of length N.
    Args:
        PP: a point process object
    Returns:
        np array"""

    empty_space_distances = []
    for i in range(0, num_sampled_points):
        x = np.random.uniform(PP.window["x_min"], PP.window["x_max"])
        y = np.random.uniform(PP.window["y_min"], PP.window["y_max"])
        d_min = float('inf')
        for j in range(0, PP.num_points):
            d = np.sqrt((PP.x[j] - x)**2 + (PP.y[j] - y)**2)
            if d < d_min and i != j:
                d_min = d

        empty_space_distances.append(d_min)

    return np.array(empty_space_distances)


def G_hat(PP: Type[Point_Process], r: np.ndarray, plot: bool = False)-> pd.DataFrame:

    """Calculates an estimate of the nearest neighbour summary statistic, G(r)
    over a range of distances, r, specified by the user.

    Args:
        PP: a point process object
        r: an np array of input values
        plot: default False, True to see plots of G(r)
    Returns:
        pd.DataFrame of floats with columns ['r', 'G_hat', 'csr']"""

    minimum_distances = min_neighbour_dists(PP)
    intensity = PP.homo_intensity_estimate()
    expected = 1 - np.exp(-intensity*np.pi*r**2)
    total_num_points = PP.num_points
    G = []
    for r_in in r:
        count = 0
        for dist in minimum_distances:
            if dist <= r_in:
                count += 1
        G.append(count/total_num_points)

    df = pd.DataFrame(data=np.column_stack((r, G, expected)),
                      columns=['r', 'G_hat', 'csr'])

    if plot:
        sbn.scatterplot(data=df, x='r', y='G_hat', alpha=0.75, label='Observed G')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR G', color='k')
        plt.xlabel('r')
        plt.ylabel('G(r)')
        plt.show()

    return df


def F_hat(PP: Type[Point_Process], r: np.ndarray, num_sampled_points: int = 1000,
          plot: bool = False) -> pd.DataFrame:

    """Calculates and plots an estimate of the  empty space summary statistics,
    F(r) over a range of distances, r, specified by the user. Uniformly samples
    `num_sampled_points` in the domain to estimate.

    Args:
        PP: a point process object
        r: an np array of input values
        num_sampled_points: default 1000, number of random points used to
                            determine F
        plot: default False, True to see plots of F(r)
    Returns:
        pd.DataFrame of floats with columns ['r', 'F_hat', 'csr']"""

    minimum_distances = min_empty_space_dists(PP, num_sampled_points)
    intensity = PP.homo_intensity_estimate()
    expected = 1 - np.exp(-intensity*np.pi*r**2)
    F = []
    for r_in in r:
        count = 0
        for dist in minimum_distances:
            if dist <= r_in:
                count += 1
        F.append(count/num_sampled_points)

    df = pd.DataFrame(data=np.column_stack((r, F, expected)),
                      columns=['r', 'F_hat', 'csr'])

    if plot:
        sbn.scatterplot(data=df, x='r', y='F_hat', alpha=0.75, label='Observed F')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR F', color='k')
        plt.xlabel('r')
        plt.ylabel('F(r)')
        plt.show()

    return df

def K_hat(PP: Type[Point_Process], r: np.ndarray, restrict_domain: bool = False,
          plot: bool = False) -> pd.DataFrame:

    """Calculates and plots an estimate of the homogeneous Ripley K function, K(r),
    over a range of distances specified by the user.

    Args:
        PP: a point process object
        r: an np array of input values
        restrict_domain: default False, restricting the domain to counter edge effects.
        plot: default False, True to see plots of K(r)
    Returns:
        pd.DataFrame of floats with columns ['r', 'K_hat', 'csr']"""

    intensity = PP.homo_intensity_estimate()
    expected = np.pi*r**2
    K = []

    distance_matrix = []
    for i in range(0, PP.num_points):
        d = []
        for j in range(0, PP.num_points):
            d = np.append(d, np.sqrt((PP.x[j]-PP.x[i])**2+(PP.y[j]-PP.y[i])**2))

        distance_matrix.append(d)

    for r_in in r:
        count = 0
        n = 0
        for p in range(0, PP.num_points):

            if restrict_domain:
                cond = ((PP.x[p] > r_in+PP.window["x_min"]) and (PP.x[p] < PP.window["x_max"]-r_in)
                and (PP.y[p] > r_in+PP.window["y_min"]) and (PP.y[p] < PP.window["y_max"]-r_in))

            else:
                cond = True

            if cond:
                n += 1
                for k in range(0, PP.num_points):

                    if(distance_matrix[p][k] <= r_in and distance_matrix[p][k] > 0):
                        count += 1
        K.append(count/(n*intensity))

    df = pd.DataFrame(data=np.column_stack((r, K, expected)),
                      columns=['r', 'K_hat', 'csr'])

    if plot:
        sbn.scatterplot(data=df, x='r', y='K_hat', alpha=0.75, label='Observed K')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR K', color='k')
        plt.xlabel('r')
        plt.ylabel('K(r)')
        plt.show()

    return df


def K_hat_fast(PP: Type[Point_Process], r: np.ndarray, plot: bool = False) -> pd.DataFrame:

    """Calculates and plots an estimate of the homogeneous Ripley K function, K(r),
    over a range of distances specified by the user. Recommended for simulation
    inference. No accounting for edge effects.

    Args:
        PP: a point process object
        r: an np array of input values
        plot: default False, True to see plots of K(r)
    Returns:
        pd.DataFrame of floats with columns ['r', 'K_hat', 'csr']"""

    # Estimate Intensity of the PP
    # intensity = PP.homo_intensity_estimate()
    expected = np.pi*r**2
    N = PP.num_points
    intensity = PP.homo_intensity_estimate()
    assert N > 0
    # Ripley K scale factor outside of sum computer once
    scale_factor = 1 / (N * intensity)


    counts = [] # count[i] is number of points within radius r_i

    # First collect array of all point to point distances. We only
    # store (N^2 - N) / 2 elements compared to the full N^2 matrix
    dists = []

    # This could be brought in as a separate 'warm-start' function
    # for K, L, O, PC. No need to compute everytime.
    for i in range(0, N):
        for j in range(i+1, N):
            dists.append(np.sqrt((PP.x[j]-PP.x[i])**2 + (PP.y[j]-PP.y[i])**2))

    # Cumalative count of points.
    count = 0
    for r_in in r:
        # If r is large enough s.t. all distances have been
        # deleted, no need to check. Just append the current
        # count.
        if len(dists) == 0:
            counts.append(count)
            continue

        # Find elements of array which represent points
        # within distance r_in. We will delete these elements
        # in order to reduce the number of points to check
        # in the next iteration.
        delete = np.where(dists <= r_in)

        # How many points were found within distance r_in
        num_deletes = len(delete[0])

        # Quick fail statement
        if num_deletes == 0:
            counts.append(count)
            continue

        # Remove the discovered points from the array and accumulate
        # the count. i.e. a point within distance r_in=4 is known to
        # be within a distance r_in <=5.
        dists = np.delete(dists, delete)

        # Add double the amount of points found to the count accounting
        # for each pair of points.
        count += 2 * num_deletes
        counts.append(count)

    # Finally convert to numpy array rather than use np.append, and scale
    # appropiately.
    K = np.array(counts) * scale_factor

    df = pd.DataFrame(data=np.column_stack((r, K, expected)),
                      columns=['r', 'K_hat', 'csr'])

    # Usual option to plot.
    if plot:
        sbn.scatterplot(data=df, x='r', y='K_hat', alpha=0.75, label='Observed K')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR K', color='k')
        plt.xlabel('r')
        plt.ylabel('K(r)')
        plt.show()

    return df


def L_hat(PP: Type[Point_Process], r: np.ndarray, restrict_domain: bool = False,
          plot: bool = False)-> pd.DataFrame:

    """Calculates and plots an estimate of the NORMALISED homogeneous Ripley K
    function, L(r), over a range of distances specified by the user.

    Args:
        PP: a point process object
        r: an np array of input values
        restrict_domain: default False, restricting the domain to counter edge
                         effects.
        plot: default False, True to see plots of L(r)
    Returns:
        pd.DataFrame of floats with columns ['r', 'L_hat', 'csr']"""

    # Call Ripley K
    K = K_hat(PP, r, restrict_domain=restrict_domain, plot=False).K_hat.to_numpy()

    # Normalise as per definition of L
    L = np.sqrt(K / np.pi) - r
    expected = np.zeros(len(r))

    df = pd.DataFrame(data=np.column_stack((r, L, expected)),
                      columns=['r', 'L_hat', 'csr'])

    # Option to plot.
    if plot:
        sbn.scatterplot(data=df, x='r', y='L_hat', alpha=0.75, label='Observed L')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR L', color='k')
        plt.xlabel('r')
        plt.ylabel('L(r)')
        plt.show()

    return df


def L_hat_fast(PP: Type[Point_Process], r: np.ndarray, plot: bool = False) -> pd.DataFrame:

    """Calculates and plots an estimate of the NORMALISED homogeneous Ripley K
    function, L(r), over a range of distances specified by the user. Recommended
    for simulation inference. No accounting for edge effects.

    Args:
        PP: a point process object
        r: an np array of input values
        plot: default False, True to see plots of L(r)
    Returns:
        pd.DataFrame of floats with columns ['r', 'L_hat', 'csr']"""

    # Call Ripley K
    K = K_hat_fast(PP, r, plot=False).K_hat.to_numpy()

    # Normalise as per definition of L
    L = np.sqrt(K / np.pi) - r
    expected = np.zeros(len(r))

    df = pd.DataFrame(data=np.column_stack((r, L, expected)),
                      columns=['r', 'L_hat', 'csr'])

    # Option to plot.
    if plot:
        sbn.scatterplot(data=df, x='r', y='L_hat', alpha=0.75, label='Observed L')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR L', color='k')
        plt.xlabel('r')
        plt.ylabel('L(r)')
        plt.show()

    return df


def O_hat(PP: Type[Point_Process], r: np.ndarray, bandwidth: float = 0.1,
          restrict_domain: bool = False, kernel: str = "BK", plot: bool = False)-> pd.DataFrame:

    """Calculates and plots an estimate of the O Ring function, O(r) over a range
    of distances specified by the user. Includes options for kernel-smoothing
    choice and whether to account for edge effects.

    Args:
        PP: a point process object
        r: an np array of input values
        bandwidth : float determining the kernel bandwidth of ring
        restrict_domain: if False edge effects will be ignored, default False
        kernel: kernel for O_hat. Default is "BK" for box kernel, which is a
        flat count. Option of "EK" for Epanechnikov kernel, a quadratic count
        which favour points closer to center of ring.
    Returns:
        pd.DataFrame of floats with columns ['r', 'O_hat', 'csr']"""

    intensity = PP.homo_intensity_estimate()
    O = []

    distance_matrix = []
    for i in range(0, PP.num_points):
        d = []
        for j in range(0, PP.num_points):
            d = np.append(d, np.sqrt((PP.x[j]-PP.x[i])**2+(PP.y[j]-PP.y[i])**2))

        distance_matrix.append(d)

    for r_in in r:
        if r_in < bandwidth: #r values below bandwidth will be excluded as O_ring not well defined
            O.append(np.nan)
            continue
        count = 0
        n = 0
        for p in range(0, PP.num_points):

            if restrict_domain:
                cond = ((PP.x[p] > r_in+PP.window["x_min"]) and (PP.x[p] < PP.window["x_max"]-r_in)
                and (PP.y[p] > r_in+PP.window["y_min"]) and (PP.y[p] < PP.window["y_max"]-r_in))

            else:
                cond = True

            if cond:
                n += 1
                for k in range(0, PP.num_points):

                    if(distance_matrix[p][k] >= r_in-bandwidth and distance_matrix[p][k] <= r_in+bandwidth
                        and distance_matrix[p][k] > 0):
                        if kernel == "BK":
                            count += (1/(2*bandwidth))
                        if kernel == "EK":
                            count += (3/(4*bandwidth))*(1-((distance_matrix[p][k]-r_in)**2/(bandwidth**2)))

        O.append(count/(2*np.pi*r_in*n))

    expected = intensity * np.ones(len(r))
    df = pd.DataFrame(data=np.column_stack((r, O, expected)),
                      columns=['r', 'O_hat', 'csr'])

    # Option to plot.
    if plot:
        sbn.scatterplot(data=df, x='r', y='O_hat', alpha=0.75, label='Observed O')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR O', color='k')
        plt.xlabel('r')
        plt.ylabel('O(r)')
        plt.show()

    return df


def PC_hat(PP: Type[Point_Process], r: np.ndarray, bandwidth: float = 0.1,
           restrict_domain: bool = False, kernel: str = "BK", plot: bool = False)-> pd.DataFrame:

    """Calculates and plots an estimate of the Pair Correlation function, PC(r)
    over a range of distances specified by the user. Includes options for kernel
    -smoothing choice and whether to account for edge effects.

    Args:
        PP: a point process object
        r: an np array of input values
        bandwidth : float determining the bandwidth of ring
        restrict_domain: if False edge effects will be ignored, default False
        kernel: kernel for O_hat. Default is "BK" for box kernel, which is a
        flat count. Option of "EK" for Epanechnikov kernel, a quadratic count
        which favour points closer to center of ring.
    Returns:
        pd.DataFrame of floats with columns ['r', 'PC_hat', 'csr']"""

    intensity = PP.homo_intensity_estimate()
    O = O_hat(PP, r, bandwidth=bandwidth, restrict_domain=restrict_domain,
            kernel=kernel).O_hat.to_numpy()
    PC = O/intensity

    expected = np.ones(len(r))

    df = pd.DataFrame(data=np.column_stack((r, PC, expected)),
                      columns=['r', 'PC_hat', 'csr'])

    # Option to plot.
    if plot:
        sbn.scatterplot(data=df, x='r', y='PC_hat', alpha=0.75, label='Observed PC')
        sbn.lineplot(data=df, x='r', y='csr', alpha=0.75, label='CSR PC', color='k')
        plt.xlabel('r')
        plt.ylabel('PC(r)')
        plt.show()

    return df


def J_hat(PP: Type[Point_Process], r: np.ndarray, num_F_sampled_points: int = 1000,
          plot: bool = False)-> pd.DataFrame:

    """Calculates and plots an estimate of the J(r) over a range of distances
    such that F(r) < 1.

                       J(r) := (1 - G(r)) / (1 - F(r))

    As the function is only defined for values of r such that F(r) < 1,
    the returned J_hat series consists of

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
        pd.DataFrame of floats with columns ['r', 'J_hat', 'csr']"""

    F = F_hat(PP, r, num_sampled_points=num_F_sampled_points, plot=False).F_hat.to_numpy()
    G = G_hat(PP, r, plot=False).G_hat.to_numpy()

    # Max value of F function
    F_max = np.max(F)
    # If F(r) != 1 anyway, proceed as normal
    if F_max < 1:
        J = (1 - G) / (1 - F)
        max_index = len(r) - 1
    # Otherwise, only calculate values for J_hat where F_hat < 1
    else:
        max_index = np.min(np.where(F == 1.0)) - 1
        J = (1 - G[:max_index]) / (1 - F[:max_index])

    # Create array of nans to be concatenated with J_hat to ensure
    # the returned J_hat array has the same length as input array r.
    expected = np.ones(max_index)
    nan_buffer = np.full(len(r) - len(J_hat), np.nan)
    J = np.concatenate((J, nan_buffer))
    expected = np.concatenate((expected, nan_buffer))

    df = pd.DataFrame(data=np.column_stack((r, J, expected)),
                      columns=['r', 'J_hat', 'csr'])

    # Option to Plot
    if plot:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
        sbn.scatterplot(r, F, alpha=0.75, label='Observed F', ax=ax1)
        sbn.scatterplot(r, G, alpha=0.75, label='Observed G', ax=ax2)

        sbn.lineplot(r[:max_index], J[:max_index], alpha=0.75,\
                     label='Observed J', linewidth=1.0, ax=ax3)
        plt.plot(r, expected, c='k', label='Expected J', alpha=0.6)
        plt.legend()

        plt.xlabel('r')
        ax1.set_ylabel('F(r)')
        ax2.set_ylabel('G(r)')
        ax3.set_ylabel('J(r)')
        fig.tight_layout()
        plt.show()

    return df


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


def CSR_hypothesis_test_dmin(observed_PP, significance_level: float, n_sim: int = 500,
                             plot_dist: bool = False)->None:

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

    if plot_dist:
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
    else:
        print('p-value is greater than or equal to {} and so the CSR'\
              .format(significance_level), 'hypothesis can not be rejected.')

    return


def simulate_summary_statistic(PP: Type[Point_Process], r: np.ndarray, summary_func: object,
                               n_sims: int = 100, plot: bool = False)->pd.DataFrame:

    """Compares metrics between a given point process, and n_sims homogenous
    poisson point processes.

    Args:
        PP: a point process object
        r: an np array of input values
        summary_func: a defined point process function (e.g. F_hat, G_hat, K_hat,
                      K_hat_fast)
        plot: default False, True to see plots of metric compared to metric for poisson with
            min and max envelopes of n_sims simulations.
    Returns:
        pd.DataFrame with columns ['r', 'summary_func', 'min', 'max', 'mean']"""

    name = summary_func.__name__
    if name == 'K_hat_fast':
        vals_column = 'K_hat'
    elif name == 'L_hat_fast':
        vals_column = 'L_hat'
    else:
        vals_column = name

    func_stack = []
    intensity = PP.homo_intensity_estimate()
    for i in range(0, n_sims):
        sim = hom_poisspp(intensity, PP.window)
        func_stack.append(summary_func(sim, r)[vals_column].to_numpy())
        b = "Loading{" + "#"*round(i*10/n_sims)+'{}%'.format(round(i*10/n_sims)*10) + "}"
        print(b, end="\r")

    func_mean = np.mean(func_stack, axis=0)
    func_max = np.max(func_stack, axis=0)
    func_min = np.min(func_stack, axis=0)
    observed = summary_func(PP, r)[vals_column].to_numpy()

    if plot:
        plt.plot(r, func_mean, c='k', label='Mean ' + name[0] + '(r) CSR Simulation', alpha=0.9)
        sbn.scatterplot(r, observed, label='Observed '+ name[0] + '(r)', s=25.0)
        plt.fill_between(r, func_min, func_max, alpha=0.3, label='Max/Min Envelope')
        plt.xlabel('r')
        plt.ylabel(name[0] + '(r)')
        plt.legend()
        plt.show()

    return pd.DataFrame(data=np.column_stack((r, observed, func_min, func_max, func_mean)),\
                        columns=['r', 'observed', 'min', 'max', 'mean'])
