import network_dynamics
import sys
import os
from datetime import datetime
import numpy as np
from pypet import Environment
from pypet import PickleResult
from pypet import pypetconstants
from pypet.utils.explore import cartesian_product
import networkx as nx
from idtxl.estimators_jidt import JidtGaussianTE
from idtxl.estimators_jidt import JidtDiscreteTE
from idtxl.data import Data
import jpype


def compute_bTE_all_pairs(traj):
    nodes_n = traj.par.topology.initial.nodes_n
    # Generate initial network
    G = network_dynamics.generate_network(traj.par.topology.initial)
    # Get adjacency matrix
    adjacency_matrix = np.array(nx.to_numpy_matrix(
        G,
        nodelist=np.array(range(0, nodes_n)),
        dtype=int))
    # Add self-loops
    np.fill_diagonal(adjacency_matrix, 1)
    # Generate initial node coupling
    coupling_matrix = network_dynamics.generate_coupling(
        traj.par.node_coupling.initial,
        adjacency_matrix)
    # Generate delay
    delay_matrices = network_dynamics.generate_delay(traj.par.delay.initial, adjacency_matrix)
    # Generate coefficient matrices
    coefficient_matrices = np.transpose(
        delay_matrices * coupling_matrix,
        (0, 2, 1))
    # Run dynamics
    time_series = network_dynamics.run_dynamics(
        traj.par.node_dynamics,
        coefficient_matrices)
    # initialise an empty data object
    dat = Data()
    # Load time series
    if traj.par.node_dynamics.model == 'boolean_random':
        normalise = False
    elif traj.par.node_dynamics.model == 'AR_gaussian_discrete':
        normalise = True
    dat = Data(
        time_series,
        dim_order='psr',
        normalise=normalise)
    data = dat.data

    # Compute empirical bTE between all pairs
    lag = 1
    history_target = traj.par.estimation.history_target
    settings = {}
    settings['source_target_delay'] = lag
    settings['history_source'] = traj.par.estimation.history_source
    if traj.par.node_dynamics.model == 'boolean_random':
        settings['history_target'] = history_target
        est = JidtDiscreteTE(settings)
    elif traj.par.node_dynamics.model == 'AR_gaussian_discrete':
        settings['history_target'] = history_target
        est = JidtGaussianTE(settings)
    bTE_empirical_matrix = np.full((nodes_n, nodes_n), np.NaN)
    for X in range(nodes_n):
        for Y in range(nodes_n):
            if (adjacency_matrix[X, Y] > 0) and (X != Y):
                bTE_empirical_matrix[X, Y] = est.estimate(data[X, :, 0], data[Y, :, 0])

    # Add results to the trajectory
    # The wildcard character $ will be replaced by the name of the current run,
    # formatted as `run_XXXXXXXX`
    traj.f_add_result(
        '$.topology.initial',
        adjacency_matrix=adjacency_matrix,
        comment='')
    traj.f_add_result(
        '$.node_coupling.initial',
        coupling_matrix=coupling_matrix,
        coefficient_matrices=coefficient_matrices,
        comment='')
    traj.f_add_result(
        '$.delay.initial',
        delay_matrices=delay_matrices,
        comment='')
    # traj.f_add_result(
    #     '$.node_dynamics',
    #     time_series=time_series,
    #     comment='')
    traj.f_add_result(
        PickleResult,
        '$.bTE',
        bTE_empirical_matrix=bTE_empirical_matrix,
        comment='')

    jSystem = jpype.JPackage("java.lang").System
    jSystem.gc()


def main():
    # Get current directory
    traj_dir = os.getcwd()
    # Read output path (if provided)
    if len(sys.argv) > 1:
        # Only use specified folder if it exists
        if os.path.isdir(sys.argv[1]):
            # Get name of directory
            traj_dir = os.path.dirname(sys.argv[1])
            # Convert to full path
            traj_dir = os.path.abspath(traj_dir)
    # Add time stamp (final '' is to make sure there is a trailing slash)
    traj_dir = os.path.join(traj_dir, datetime.now().strftime("%Y_%m_%d_%Hh%Mm%Ss"), '')
    # Create directory with time stamp
    os.makedirs(traj_dir)
    # Change current directory to the one containing the trajectory files
    os.chdir(traj_dir)
    print('Trajectory and results will be stored to: {0}'.format(traj_dir))

    # Create an environment that handles running.
    # Let's enable multiprocessing with scoop:
    env = Environment(
        trajectory='traj',
        comment='',
        add_time=False,
        log_config='DEFAULT',
        log_stdout=True,  # log everything that is printed, will make the log file HUGE
        filename=traj_dir,  # filename or just folder (name will be automatic in this case)
        multiproc=True,
        use_scoop=True,
        wrap_mode=pypetconstants.WRAP_MODE_LOCAL,
        memory_cap=10,
        swap_cap=1
    )
    traj = env.trajectory

    # -------------------------------------------------------------------
    # Add config parameters (those that DO NOT influence the final result of the experiment)
    traj.f_add_config('parallel_target_analysis', False, comment='Analyse targets in parallel')
    # -------------------------------------------------------------------
    # Parameters characterizing the initial topology of the network
    #traj.f_add_parameter('topology.initial.model', 'BA')
    traj.f_add_parameter('topology.initial.model', 'WS')
    traj.f_add_parameter('topology.initial.nodes_n', 5, comment='Number of nodes')
    traj.f_add_parameter('topology.initial.WS_k', 4, comment='Number of neighbours (and mean degree) in the Watts-Strogatz model')
    traj.f_add_parameter('topology.initial.WS_p', 0.0, comment='Rewiring probability in the Watts-Strogatz model')
    #traj.f_add_parameter('topology.initial.BA_m', 1, comment='Number of edges to attach from a new node to existing nodes in the Barabási–Albert model')


    # -------------------------------------------------------------------
    # Parameters characterizing the coupling between the nodes
    traj.f_add_parameter('node_coupling.initial.model', 'linear', comment='Linear coupling model: the input to each target node is the weighted sum of the outputs of its source nodes')
    traj.f_add_parameter('node_coupling.initial.weight_distribution', 'fixed')
    traj.f_add_parameter('node_coupling.initial.fixed_coupling', 0.15)

    # -------------------------------------------------------------------
    # Parameters characterizing the delay
    traj.f_add_parameter('delay.initial.distribution', 'uniform')
    traj.f_add_parameter('delay.initial.delay_links_n_max', 1, comment='Maximum number of delay links')
    traj.f_add_parameter('delay.initial.delay_min', 1, comment='')
    traj.f_add_parameter('delay.initial.delay_max', 1, comment='')
    traj.f_add_parameter('delay.initial.delay_self', 1, comment='')

    # -------------------------------------------------------------------
    # Parameters characterizing the estimator
    traj.f_add_parameter('estimation.history_source', 1, comment='Embedding length for the source')
    traj.f_add_parameter('estimation.history_target', 14, comment='Embedding length for the target')


    # -------------------------------------------------------------------
    # Parameters characterizing the dynamics of the nodes
    traj.f_add_parameter('node_dynamics.model', 'AR_gaussian_discrete')
    #traj.f_add_parameter('node_dynamics.model', 'boolean_random')
    traj.f_add_parameter('node_dynamics.samples_n', 100, comment='Number of samples (observations) to record')
    traj.f_add_parameter('node_dynamics.samples_transient_n', 1000 * traj.topology.initial.nodes_n, comment='Number of initial samples (observations) to skip to leave out the transient')
    traj.f_add_parameter('node_dynamics.replications', 1, comment='Number of replications (trials) to record')
    traj.f_add_parameter('node_dynamics.noise_std', 1, comment='Standard deviation of Gaussian noise')
    #traj.f_add_parameter('node_dynamics.RBN_r', 0.5, comment='Activity level (i.e. probability of state "1") in Boolean dynamics')
    #traj.f_add_parameter('node_dynamics.noise_flip_p', 0.005, comment='Probability of flipping bit in Boolean dynamics')

    # -------------------------------------------------------------------
    # Parameters characterizing the repetitions of the same run
    traj.f_add_parameter('repetition_i', 0, comment='Index of the current repetition') # Normally starts from 0

    # -------------------------------------------------------------------
    # Define parameter combinations to explore (a trajectory in
    # the parameter space)
    # The second argument, the tuple, specifies the order of the cartesian product,
    # The variable on the right most side changes fastest and defines the
    # 'inner for-loop' of the cartesian product
    explore_dict = cartesian_product(
        {
            'node_coupling.initial.weight_distribution': ['fixed'],
            'repetition_i': np.arange(0, 10, 1).tolist(),
            'topology.initial.nodes_n': np.arange(100, 100+1, 30).tolist(),
            'node_dynamics.samples_n': np.array([100000]).tolist(),
            'topology.initial.WS_p': np.around(np.logspace(-2.2, 0, 10), decimals=4).tolist(),
        },
        (
            'node_coupling.initial.weight_distribution',
            'node_dynamics.samples_n',
            'topology.initial.nodes_n',
            'topology.initial.WS_p',
            'repetition_i',
        )
    )
    print(explore_dict)
    traj.f_explore(explore_dict)
    # Run the experiment
    env.run(compute_bTE_all_pairs)
    # Check that all runs are completed
    assert traj.f_is_completed()
    # Finally disable logging and close all log-files
    env.disable_logging()


if __name__ == '__main__':
    # This will execute the main function in case the script is called from the one true
    # main process and not from a child processes spawned by your environment.
    # Necessary for multiprocessing under Windows.
    main()
