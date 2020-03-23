import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from cycler import cycler


def get_parameters_explored(first_not_explored, ignore_par):
    # Determine explore parameters
    # WARNING: this method is prone to errors, when the columns are rearranged,
    # however, in such cases the "assertion" sanity checks below will fail, so I
    # will know
    parameters_explored = df.loc[[], :first_not_explored].keys().tolist()[:-1]
    if (len(ignore_par) > 0 and bool(
        set(ignore_par).intersection(parameters_explored))):
        print('\nWARNING: Ignoring the following explored parameters: {0}'.format(
            ignore_par))
        parameters_explored = [
            item for item in parameters_explored if item not in ignore_par]
    print('parameters_explored = : {0}\n'.format(parameters_explored))
    return parameters_explored


def my_subplots(subplots_v, subplots_h, sharex=True, sharey=True):
    fig, axs = plt.subplots(
        subplots_v,
        subplots_h,
        sharex=sharex,
        sharey=sharey)
    # Avoid indexing issues if only one row or one column
    if subplots_v == 1 and subplots_h == 1:
        axs = np.array([axs])
        print('One vertical and one horizontal subplots only')
    else:
        if subplots_v == 1:
            print('One vertical subplot only')
            axs = np.array([axs])
            if subplots_h == 1:
                print('One horizontal subplot only')
                axs = np.array([axs])
                axs = np.transpose(axs)
    # Set figure size
    fig.set_figheight(3 * subplots_v)
    if subplots_h > 1:
        fig.set_figwidth(4 * subplots_h)

    return fig, axs


def save_figures_to_pdf(figures, path, pdf_metadata={}):
    """
    Create a PDF file with several pages.
    Also add PDF file metadata if provided.
    Some possible metadata keys: ['Title', 'Author', 'Subject', 'Keywords',
        'CreationDate', 'ModDate']
    """

    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed
    # properly at the end of the block, even if an Exception occurs.
    with PdfPages(path) as pdf:
        for (fig_i, fig) in enumerate(fig_list):

            # Set figure as current figure
            plt.figure(fig.number)

            # # Set figure size
            # ndim = axs.ndim
            # if ndim == 1:
            #     fig.set_figheight(3 * len(axs))
            # if ndim > 1:
            #     fig.set_figheight(3 * len(axs[:, 0]))
            #     fig.set_figwidth(4 * len(axs[0, :]))
            
            # Set tight layout to avoid overlap between subplots
            fig.tight_layout()

            # Add a pdf note to attach metadata to a page
            # pdf.attach_note('note...')

            # Save figure to PDF page
            pdf.savefig(fig)

            # Also export single figures
            # figure_format = 'PNG'
            # export_path = os.path.join(traj_dir, '{0}.{1}'.format(fig.number, figure_format))
            # plt.savefig(export_path, format=figure_format)

            print('figure {0} saved'.format(fig_i + 1))

        # Set PDF file metadata via the PdfPages object
        d = pdf.infodict()
        for key in pdf_metadata.keys():
            d[key] = pdf_metadata.get(key, '')


def check_remaining_dimensions(df_keys_remaining, parameters_to_average):
    # Ensure that only the desired parameters are aggregated or averaged
    averaging_set = set.intersection(
        set(df_keys_remaining),
        set(parameters_explored))
    assert averaging_set == parameters_to_average, (
        "Attempting to average over {0}, which differs from the "
        "specified set {1}".format(averaging_set, parameters_to_average))


def plot_bTE_empirical_vs_WS_p():
    # Plot bTE vs. Watts-Strogatz rewiring prob
    # Legend: methods for computing bTE
    bTE_methods = ['bTE_empirical_causal_vars',]
    fig, axs = my_subplots(1, 1, sharex=True, sharey=True)
    # Select data of interest
    df_interest = df[parameters_explored + bTE_methods]
    df_interest = df_interest.loc[
        df_interest['nodes_n'] == N_interest].drop('nodes_n', 1)
    df_interest = df_interest.loc[
        df_interest['samples_n'] == T_interest].drop('samples_n', 1)
    df_interest = df_interest.loc[
        df_interest['weight_distribution'] == weight_interest].drop('weight_distribution', 1)
    # Choose which of the explored parameters will be aggregated or averaged
    parameters_to_average = {'repetition_i'}
    for (method_i, bTE_method) in enumerate(bTE_methods):
        df_interest[bTE_method] = df_interest[bTE_method].apply(np.nanmean)
        # Group by WS_p
        df_aggregate = df_interest.groupby('WS_p').agg(
            lambda x: list(x))
        # Ensure that only the desired parameters are aggregated or averaged
        df_keys_remaining = df_aggregate.columns.get_level_values(0)
        check_remaining_dimensions(df_keys_remaining, parameters_to_average)
        # Plot average path lenght inferred
        axs[0].plot(
            WS_p_range,
            [np.nanmean(df_aggregate[bTE_method][WS_p])
             for WS_p in WS_p_range],
            marker=markers_default[method_i],
            label=bTE_methods_legends[bTE_method])
        axs[0].scatter(
            df_interest['WS_p'].values.astype(float).tolist(),
            df_interest[bTE_method].values,
            alpha=0.3,
            marker=markers_default[method_i])
    # Set axes properties
    axs[0].set_xlabel(r'$\text{Rewiring probability }(\gamma)$')
    axs[0].set_ylabel(r'$\text{Transfer Entropy}$')
    axs[0].set_xscale('log')
    axs[0].yaxis.set_ticks_position('both')
    axs[0].set_xlim(left=0.002)

    fig.set_figheight(3)
    fig.set_figwidth(4)

    return fig


# Choose folder
traj_dir = 'TE_from_couplings_RBF_WS_sweep_noise0.005_100nodes_100000samples_20rep_history14'

# Set up plot style
# use latex based font rendering
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command
# Use "matplotlib.rcdefaults()" to restore the default plot style"
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['lines.markersize'] = 4
mpl.rcParams['errorbar.capsize'] = 3
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.labelspacing'] = 0.5  # 0.5 default
mpl.rcParams['legend.handletextpad'] = 0.8  # 0.8 default

# Colours
# Other colorblind palette proposed in https://github.com/matplotlib/matplotlib/issues/9460
colors_petroff = [
    '#7B85D4',
    '#f37738',
    '#83c995',
    '#d7369e',
    '#c4c9d8',
    '#859795',
    '#e9d043',
    '#ad5b50']

# Markers
markers_default = ['x', 'o', '^', 's', 'D', 'v', 'h', '*']
cycler_default = cycler(
    color=colors_petroff,
    marker=markers_default)

# Load DataFrame
df = pd.read_pickle(os.path.join(traj_dir, 'postprocessing.pkl'))

# Initialise empty figure and axes lists
fig_list = []

# Select value of interest (for those plots where only one value is used)
N_interest = 100
T_interest = 100000
weight_interest = 'fixed'
first_not_explored = 'bTE_empirical_causal_vars'
ignore_par = {
    'jidt_threads_n',
    'n_perm_max_stat',
    'n_perm_min_stat',
    'n_perm_max_seq',
    'n_perm_omnibus',
    }
parameters_explored = get_parameters_explored(first_not_explored, ignore_par)
# Get parameter ranges
nodes_n_range = np.unique(df['nodes_n']).astype(int)
samples_n_range = np.unique(df['samples_n']).astype(int)
weight_distributions = np.unique(df['weight_distribution'])

bTE_methods_legends = {
    'bTE_empirical_causal_vars': 'Empirical',
    'bTE_theoretical_causal_vars': 'Theoretical',
    'bTE_approx2_causal_vars': 'Order 2',
    'bTE_approx4_causal_vars': r'All motifs up to $\mathcal{O}(\|C\|^4)$',
    'bTE_motifs_acde_causal_vars': 'Motifs a + c + d + e',
    'bTE_motifs_ae_causal_vars': 'Motifs a + e',
}

WS_p_range = np.unique(df['WS_p']).astype(float)
fig_list.append(plot_bTE_empirical_vs_WS_p())

# Save figures to PDF file
pdf_metadata = {}
pdf_path = os.path.join(traj_dir, 'figures.pdf')
save_figures_to_pdf(fig_list, pdf_path)
