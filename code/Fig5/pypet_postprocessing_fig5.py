from pypet import Trajectory
import os
import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from scipy.linalg import solve_discrete_lyapunov

fdr = False
debug_mode = False
save_results = True

# Load the trajectory from the hdf5 file
# Only load parameters, results will be loaded at runtime (auto loading)
traj_dir = 'TE_from_couplings_RBF_WS_sweep_noise0.005_100nodes_100000samples_20rep_history14'
traj_filename = 'traj.hdf5'
traj_fullpath = os.path.join(traj_dir, traj_filename)
traj = Trajectory()
traj.f_load(
    filename=traj_fullpath,
    index=0,
    load_parameters=2,
    load_results=0,
    load_derived_parameters=0,
    force=True)
# Turn on auto loading
traj.v_auto_load = True

# Count number of runs
runs_n = len(traj.f_get_run_names())
print('Number of runs = {0}'.format(runs_n))

# Get list of explored parameters
parameters_explored = [str.split(par, '.').pop() for par in (
    traj.f_get_explored_parameters())]

# Initialise analysis summary table
# (it is important that the columns with the explored parameters
# preceed the ones with the results)
df = pd.DataFrame(
    index=traj.f_get_run_names(),
    columns=parameters_explored + [
        'bTE_empirical_causal_vars',
        'bTE_approx2_causal_vars',
        'bTE_approx4_causal_vars',
        'bTE_motifs_acde_causal_vars',
        'bTE_motifs_ae_causal_vars',
        'bTE_theoretical_causal_vars',
        ],
    dtype=object)
# Add other useful columns to the DataFrame
if 'weight_distribution' not in parameters_explored:
    df['weight_distribution'] = np.NaN

# Loop over runs
for run_name in traj.f_get_run_names():

    # Make trajectory behave like a particular single run:
    # all explored parameter’s values will be set to the corresponding
    # values of one particular run.
    traj.f_set_crun(run_name)

    print('post-processing of {0} in progress...'.format(run_name))

    # Fill in current explored parameter values
    for par in parameters_explored:
        df.loc[run_name, par] = traj.parameters[par]
        if debug_mode:
            print('{0} = {1}'.format(par, traj.parameters[par]))
    if 'weight_distribution' not in parameters_explored:
        current_weight_distribution = traj.parameters[
            'node_coupling.initial.weight_distribution']
        df.loc[run_name, 'weight_distribution'] = current_weight_distribution
        if debug_mode:
            print('weight_distribution added to DataFrame: {0}'.format(
                current_weight_distribution))

    nodes_n = traj.parameters.topology.initial['nodes_n']

    # Load results object
    res = traj.results[run_name].bTE
    # Get real adjacency matrix
    adjacency_matrix_real = np.asarray(
        traj.results[run_name].topology.initial.adjacency_matrix > 0,
        dtype=float).astype(float)
    # Remove self-loops from real and inferred adjacency matrices
    np.fill_diagonal(adjacency_matrix_real, 0)

    # Get real and inferred delay matrices
    history_target = 5
    delay_matrices_real = np.zeros(
        shape=(history_target, nodes_n, nodes_n))
    delay_matrices_real[0, :, :] = (
        traj.results[run_name].delay.initial.delay_matrices
        ).astype(float)[0, :, :]

    # Get max delays
    delay_max_real = np.shape(delay_matrices_real)[0]
    # Load real coupling matrix
    coupling_matrix_real = (
        traj.results[run_name].node_coupling.initial.coupling_matrix)

    # print(np.transpose(coupling_matrix_real))

    # -------------------------------------------------------------------------
    # region TE
    # -------------------------------------------------------------------------

    # Approximate using motifs
    bTE_empirical_matrix = res.bTE_empirical_matrix
    bTE_empirical_matrix[adjacency_matrix_real == 0] = np.NaN
    bTE_approx2_matrix = np.full((nodes_n, nodes_n), np.NaN)
    bTE_approx4_matrix = np.full((nodes_n, nodes_n), np.NaN)
    bTE_motifs_acde_matrix = np.full((nodes_n, nodes_n), np.NaN)
    bTE_motifs_ae_matrix = np.full((nodes_n, nodes_n), np.NaN)

    C = coupling_matrix_real
    CC = C.dot(C)
    CTC = (C.T).dot(C)
    CTCC = (C.T).dot(CC)
    for X in range(nodes_n):
        for Y in range(nodes_n):
            if adjacency_matrix_real[X, Y] > 0:
                CXX = C[X, X]
                CXY = C[X, Y]
                CYX = C[Y, X]
                CYY = C[Y, Y]
                # Compute approximation order 2
                bTE_approx2_matrix[X, Y] = 0.5 * CXY ** 2
                # Compute approximation order 4
                bTE_approx4_matrix[X, Y] = 0.5 * (
                    + 1 * CXY ** 2
                    + 2 * CXY * CTCC[X, Y]
                    - 1 * (CXY ** 2) * CTC[X, X]
                    - 1 * (CXY ** 2) * CTC[Y, Y]
                    - 2 * CXY * CYY * CTC[Y, X]
                    - 2 * CXY * CYX * CC[Y, Y]
                    + 0.5 * (CXY ** 4)
                    + (CXY ** 2) * (CYY ** 2)
                    + (CXY ** 2) * (CYX ** 2)
                    + 2 * CXY * CYX * (CYY ** 2))
                # Compute approximation using only
                # motifs a, c, d, and e
                bTE_motifs_acde_matrix[X, Y] = (
                    # a
                    + 0.5 * (CXY ** 2)
                    - 0.25 * (CXY ** 4)
                    # c
                    + 0.5 * (CXY ** 2) * CTC[X, X]
                    - 0.5 * (CXY ** 2) * (CXX ** 2)
                    - 0.5 * (CXY ** 2) * (CYX ** 2)
                    # d
                    - 0.5 * (CXY ** 2) * CTC[Y, Y]
                    + 0.5 * (CXY ** 2) * (CYY ** 2)
                    + 0.5 * (CXY ** 2) * (CXY ** 2)
                    # e
                    + 0.5 * (CXY ** 2) * (CXX ** 2)
                    )
                # Compute approximation using only
                # motifs a and e
                bTE_motifs_ae_matrix[X, Y] = (
                    # a
                    + 0.5 * (CXY ** 2)
                    - 0.25 * (CXY ** 4)
                    # e
                    + 0.5 * (CXY ** 2) * (CXX ** 2)
                    )
    if debug_mode:
        print('adjacency matrix real:\n {0}'.format(adjacency_matrix_real))
        print('bTE_empirical_matrix:\n {0}'.format(bTE_empirical_matrix))
        print('bTE_approx2_matrix:\n {0}'.format(bTE_approx2_matrix))
        print('bTE_approx4_matrix:\n {0}'.format(bTE_approx4_matrix))
        print('bTE_motifs_acde_matrix:\n {0}'.format(bTE_motifs_acde_matrix))
        print('bTE_motifs_ae_matrix:\n {0}'.format(bTE_motifs_ae_matrix))
    # Add results to DataFrame
    df.loc[run_name, 'bTE_empirical_causal_vars'] = (
        bTE_empirical_matrix)
    df.loc[run_name, 'bTE_approx2_causal_vars'] = (
        bTE_approx2_matrix)
    df.loc[run_name, 'bTE_approx4_causal_vars'] = (
        bTE_approx4_matrix)
    df.loc[run_name, 'bTE_motifs_acde_causal_vars'] = (
        bTE_motifs_acde_matrix)
    df.loc[run_name, 'bTE_motifs_ae_causal_vars'] = (
        bTE_motifs_ae_matrix)

    # Initialise vectors
    bTE_theoretical_causal_vars = np.full(
        (delay_max_real, nodes_n, nodes_n), np.NaN)
    # Compute theoretical info-theoretic measures if VAR process
    if traj.par.node_dynamics.model == 'AR_gaussian_discrete':
        # Recover coefficient matrices
        coefficient_matrices_real = np.transpose(
            delay_matrices_real * coupling_matrix_real,
            (0, 2, 1))
        # Build VAR reduced form
        # See Appendix Faes et al. (PRE, 2015, doi: 10.1103/PhysRevE.91.032904)
        lags = history_target
        VAR_reduced_form = np.zeros((nodes_n * lags, nodes_n * lags))
        VAR_reduced_form[0:nodes_n, 0:nodes_n*delay_max_real] = np.reshape(
            np.transpose(coefficient_matrices_real, (1, 0, 2)),
            [nodes_n, nodes_n * delay_max_real])
        VAR_reduced_form[nodes_n:, 0:-nodes_n] = np.eye(nodes_n * (lags - 1))

        # Use VAR reduced form to compute the spectral radius
        radius = max(np.abs(np.linalg.eigvals(VAR_reduced_form)))
        print('spectral radius = {0}'.format(radius))
        assert radius < 1

        # Recover process noise covariance matrix
        variance = traj.parameters.node_dynamics.noise_std ** 2
        process_noise_cov = variance * np.eye(nodes_n, nodes_n)
        # Pad noise covariance matrix with zeros along both dimensions
        # (to the right and at the bottom)
        VAR_reduced_noise_cov = block_diag(
            process_noise_cov,
            np.zeros((nodes_n * (lags - 1), nodes_n * (lags - 1))))
        # Compute lagged cov matrix by solving discrete Lyapunov equation
        # cov = VAR_reduced_form * cov * VAR_reduced_form.T + noise_cov
        # (See scipy documentation for 'solve_discrete_lyapunov' function)
        VAR_reduced_cov = solve_discrete_lyapunov(
            VAR_reduced_form,
            VAR_reduced_noise_cov,
            method='bilinear')  # use 'bilinear' if 'direct' fails or too slow
        # Check solution
        abs_diff = np.max(np.abs((
            VAR_reduced_cov
            - VAR_reduced_form.dot(VAR_reduced_cov).dot(VAR_reduced_form.T)
            - VAR_reduced_noise_cov)))
        assert abs_diff < 10**-6, "large absolute error = {}".format(abs_diff)

    def partial_variance(variance, cross_cov, cov):
        inverse_cov = np.linalg.inv(cov)
        return variance - np.dot(np.dot(cross_cov, inverse_cov), cross_cov.T)

    def compute_CMI(X_IDs, Y_IDs, Z_IDs):
        # Compute theoretical CMI from covariance matrix of VAR (reduced form)
        if len(X_IDs) > 0 and len(Y_IDs) > 0:
            # Use np._ix to extract submatrix from row and column indices
            Y_cov = VAR_reduced_cov[np.ix_(Y_IDs, Y_IDs)]
            Y_Z_crosscov = VAR_reduced_cov[np.ix_(Y_IDs, Z_IDs)]
            Z_cov = VAR_reduced_cov[np.ix_(Z_IDs, Z_IDs)]
            numerator = partial_variance(Y_cov, Y_Z_crosscov, Z_cov)
            ZX_IDs = Z_IDs + X_IDs
            Y_ZX_crosscov = VAR_reduced_cov[np.ix_(Y_IDs, ZX_IDs)]
            ZX_cov = VAR_reduced_cov[np.ix_(ZX_IDs, ZX_IDs)]
            denominator = partial_variance(Y_cov, Y_ZX_crosscov, ZX_cov)
            CMI = 0.5 * np.log(numerator / denominator)
        else:
            print('Empty source or target set. Will return CMI=NaN')
            CMI = np.NaN
        return CMI

    def compute_bTE(source_IDs, target_IDs, targetPast_IDs):
        # Compute theoretical apparent transfer entropy from covariance
        # matrix of the VAR (reduced form)
        bTE = np.full(shape=(len(source_IDs)), fill_value=np.NaN)
        for i, s in enumerate(source_IDs):
            bTE[i] = compute_CMI(
                [s],
                target_IDs,
                targetPast_IDs)
        return bTE

    for t in range(nodes_n):
        if debug_mode:
            print('\nTarget = {0}'.format(t))

        # Compute theoretical TE if VAR process
        if traj.par.node_dynamics.model == 'AR_gaussian_discrete':
            if debug_mode:
                print('Computing theoretical TE for VAR process...')
            target_IDs = [t]

            # Compute theoretical TE on causal variables
            # get all causal vars as (variable, delay) pairs
            vars_causal = np.rot90(
                np.array(delay_matrices_real[:, :, t].nonzero()),
                axes=(1, 0))  # rotate clockwise
            # Split into source and target variables
            target_self = vars_causal[:, 0] == t
            others = np.invert(target_self)
            targetPast_causal = vars_causal[target_self]
            source_causal = vars_causal[others]
            source_IDs_causal = (
                (source_causal[:, 1] + 1) * nodes_n
                + source_causal[:, 0]).tolist()
            # targetPast_IDs_causal = (
            #     (targetPast_causal[:, 1] + 1) * nodes_n
            #     + targetPast_causal[:, 0]).tolist()
            targetPast_IDs_causal = []
            for lag in range(history_target - 1):
                targetPast_IDs_causal = targetPast_IDs_causal + (
                    (targetPast_causal[:, 1] + 1 + lag) * nodes_n
                    + targetPast_causal[:, 0]).tolist()
            if debug_mode:
                print('Causal variables:')
                print('source_IDs_causal = {0}'.format(source_IDs_causal))
                print('targetPast_IDs_causal = {0}'.format(
                    targetPast_IDs_causal))
            # Compute apparent TE
            bTE_theoretical_causal_vars[
                source_causal[:, 1], source_causal[:, 0], t] = compute_bTE(
                source_IDs_causal,
                target_IDs,
                targetPast_IDs_causal)
            if debug_mode:
                print('bTE_theoretical_causal_vars = {0}'.format(
                    bTE_theoretical_causal_vars[
                        source_causal[:, 1], source_causal[:, 0], t]))
        else:
            bTE_theoretical_causal_vars = np.NaN

    if traj.par.node_dynamics.model == 'AR_gaussian_discrete':
        bTE_theoretical_causal_vars = bTE_theoretical_causal_vars[0]
    # Add results to DataFrame
    df.loc[
        run_name,
        'bTE_theoretical_causal_vars'
        ] = bTE_theoretical_causal_vars

    # -------------------------------------------------------------------------
    # endregion
    # -------------------------------------------------------------------------

# Reset trajectory to the default settings, to release its belief to
# be the last run:
traj.f_restore_default()

if save_results:
    # Save DataFrame
    if fdr:
        df.to_pickle(os.path.join(traj_dir, 'postprocessing_fdr.pkl'))
    else:
        df.to_pickle(os.path.join(traj_dir, 'postprocessing.pkl'))
else:
    print('\nWARNING: Postprocessing DataFrame NOT saved!')
