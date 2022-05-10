# PCMCI Python script
# Heili Lowman
# April 11, 2022

# Overview of package available at:
# https://github.com/jakobrunge/tigramite

# PCMCI assumes "Causal stationarity, no contemporaneous causal links, no hidden variables"
# PCMCI outputs "Directed lagged links, undirected contemporaneous links (for tau_min=0)"

# Using the tutorial found at:
# github.com/jakobrunge/tigramite/blob/master/tutorials/tigramite_tutorial_basics.ipynb

# But instead of simulating data as they do,
# the following uses an existing USGS dataset.

# Imports
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
# %matplotlib inline
## use '%matplotlib notebook' for interactive figures
# plt.style.use('ggplot')
import sklearn

# Run 'pip install tigramite' in console to install the package
# you may need to restart python to use it though.
import tigramite
from tigramite import data_processing as pp
from tigramite.toymodels import structural_causal_processes as toys

from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb

# Save the dataset in as a variable.
ts_df = pd.read_csv("pcmci_test_dataset.csv")

# Trim dataset for easier use and viewable plotting.
ts_df_ed = ts_df[["GPP_raw", "ER_raw", "temp.water", "discharge", "PAR_sum"]]

# Make it an array? Yes, this is necessary.
ts_array = ts_df_ed.to_numpy()

# Save the dataframe in the appropriate tigramite format.
variable_names = ts_df_ed.columns
ts_dataframe = pp.DataFrame(ts_array,
                            datatime = np.arange(len(ts_array)),
                            var_names = variable_names)

# Plot the time series
tp.plot_timeseries(ts_dataframe); plt.show()
# So, if trying to mimic the structure of the tutorial:
    # X0 = GPP
    # X1 = ER
    # X2 = Temperature
    # X3 = Discharge
    # X4 = Light (PAR)

# Check to make sure all the parts are there.
ts_dataframe.var_names # it's not empty
var_names = ts_dataframe.var_names # create for labeled plotting below
ts_dataframe.values # it has variables to call
ts_dataframe.datatime # and it is indexed

# Implement linear partial correlation.
parcorr = ParCorr(significance = 'analytic')

# Initialize the PCMCI method.
pcmci = PCMCI(
    dataframe=ts_dataframe,
    cond_ind_test=parcorr,
    verbosity=1)

# Plot lagged correlations
# tau_max refers to the maximal time lag
# Changed tau to 5, because I don't think it'll be much more than 2
correlations = pcmci.get_lagged_dependencies(tau_max=5, val_only=True)['val_matrix']
lag_func_matrix = tp.plot_lagfuncs(val_matrix=correlations, setup_args={'var_names':var_names}); plt.show()

# Use scatterplots to check and see if dependencies are truly linear.
# scatter_lags sets the lag to use for every pair of variables to be their maximal absolute value from calculations above.
scatter_lags = np.argmax(np.abs(correlations), axis=2)
tp.plot_scatterplots(dataframe=ts_dataframe, add_scatterplot_args={'scatter_lags':scatter_lags}); plt.show()

# Just going to keep tau_max at 5 for now, but this coudl certainly be changed.
# Set pc_alpha to None so the PCMCI will choose the significance level in the condition-selection step itself.
# and set alpha_level to 0.05 to be the threshold for significance in plotting (less conservative than theirs).
pcmci.verbosity = 1
results = pcmci.run_pcmci(tau_max=5, pc_alpha=None, alpha_level=0.05)

# ##
# ## Step 1: PC1 algorithm with lagged conditions
# ##

# Parameters:
# independence test = par_corr
# tau_min = 1
# tau_max = 5
# pc_alpha = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
# max_conds_dim = None
# max_combinations = 1



# ## Resulting lagged parent (super)sets:

#     Variable GPP_raw has 8 link(s):
#     [pc_alpha = 0.3]
#         (GPP_raw -1): max_pval = 0.00000, min_val =  0.166
#         (GPP_raw -2): max_pval = 0.00000, min_val =  0.130
#         (GPP_raw -5): max_pval = 0.00001, min_val =  0.117
#         (GPP_raw -3): max_pval = 0.00005, min_val =  0.109
#         (discharge -1): max_pval = 0.00037, min_val = -0.096
#         (GPP_raw -4): max_pval = 0.00658, min_val =  0.073
#         (PAR_sum -3): max_pval = 0.29255, min_val =  0.028
#         (discharge -5): max_pval = 0.29370, min_val =  0.028

#     Variable ER_raw has 9 link(s):
#     [pc_alpha = 0.4]
#         (ER_raw -1): max_pval = 0.00000, min_val =  0.371
#         (ER_raw -2): max_pval = 0.00000, min_val =  0.196
#         (discharge -1): max_pval = 0.00000, min_val =  0.180
#         (ER_raw -3): max_pval = 0.00000, min_val =  0.173
#         (GPP_raw -1): max_pval = 0.00008, min_val = -0.106
#         (discharge -4): max_pval = 0.09758, min_val =  0.045
#         (ER_raw -4): max_pval = 0.10550, min_val =  0.044
#         (GPP_raw -5): max_pval = 0.28087, min_val =  0.029
#         (GPP_raw -2): max_pval = 0.30627, min_val = -0.028

#     Variable temp.water has 8 link(s):
#     [pc_alpha = 0.1]
#         (temp.water -1): max_pval = 0.00000, min_val =  0.761
#         (temp.water -2): max_pval = 0.00000, min_val = -0.213
#         (PAR_sum -1): max_pval = 0.00000, min_val =  0.170
#         (temp.water -4): max_pval = 0.00258, min_val =  0.081
#         (GPP_raw -1): max_pval = 0.00748, min_val =  0.072
#         (temp.water -5): max_pval = 0.01883, min_val =  0.064
#         (PAR_sum -2): max_pval = 0.02949, min_val =  0.059
#         (discharge -1): max_pval = 0.07246, min_val = -0.049

#     Variable discharge has 10 link(s):
#     [pc_alpha = 0.4]
#         (discharge -1): max_pval = 0.00000, min_val =  0.648
#         (GPP_raw -1): max_pval = 0.00326, min_val = -0.079
#         (discharge -2): max_pval = 0.01072, min_val = -0.069
#         (temp.water -1): max_pval = 0.06032, min_val =  0.051
#         (discharge -4): max_pval = 0.11624, min_val =  0.042
#         (temp.water -4): max_pval = 0.23025, min_val = -0.032
#         (ER_raw -2): max_pval = 0.25461, min_val =  0.031
#         (temp.water -2): max_pval = 0.29485, min_val = -0.028
#         (GPP_raw -2): max_pval = 0.31318, min_val =  0.027
#         (PAR_sum -5): max_pval = 0.39081, min_val =  0.023

#     Variable PAR_sum has 6 link(s):
#     [pc_alpha = 0.2]
#         (PAR_sum -1): max_pval = 0.00000, min_val =  0.315
#         (PAR_sum -4): max_pval = 0.00000, min_val =  0.171
#         (PAR_sum -3): max_pval = 0.00000, min_val =  0.146
#         (PAR_sum -5): max_pval = 0.00000, min_val =  0.130
#         (PAR_sum -2): max_pval = 0.00008, min_val =  0.106
#         (discharge -1): max_pval = 0.18160, min_val =  0.036

# ##
# ## Step 2: MCI algorithm
# ##

# Parameters:

# independence test = par_corr
# tau_min = 0
# tau_max = 5
# max_conds_py = None
# max_conds_px = None

# ## Significant links at alpha = 0.05:

#     Variable GPP_raw has 9 link(s):
#         (ER_raw  0): pval = 0.00000 | val =  0.249 | unoriented link
#         (discharge  0): pval = 0.00000 | val = -0.204 | unoriented link
#         (GPP_raw -1): pval = 0.00000 | val =  0.156
#         (PAR_sum  0): pval = 0.00000 | val =  0.144 | unoriented link
#         (GPP_raw -2): pval = 0.00001 | val =  0.120
#         (GPP_raw -5): pval = 0.00009 | val =  0.106
#         (GPP_raw -3): pval = 0.00016 | val =  0.102
#         (discharge -1): pval = 0.00216 | val = -0.083
#         (GPP_raw -4): pval = 0.00612 | val =  0.074

#     Variable ER_raw has 10 link(s):
#         (ER_raw -1): pval = 0.00000 | val =  0.375
#         (GPP_raw  0): pval = 0.00000 | val =  0.249 | unoriented link
#         (ER_raw -2): pval = 0.00000 | val =  0.195
#         (discharge -1): pval = 0.00000 | val =  0.184
#         (ER_raw -3): pval = 0.00000 | val =  0.165
#         (GPP_raw -1): pval = 0.00000 | val = -0.126
#         (discharge -4): pval = 0.00031 | val = -0.098
#         (discharge  0): pval = 0.00132 | val = -0.087 | unoriented link
#         (discharge -2): pval = 0.01276 | val = -0.068
#         (GPP_raw -2): pval = 0.01419 | val = -0.067

#     Variable temp.water has 9 link(s):
#         (temp.water -1): pval = 0.00000 | val =  0.758
#         (temp.water -2): pval = 0.00000 | val = -0.316
#         (PAR_sum -1): pval = 0.00000 | val =  0.151
#         (temp.water -3): pval = 0.00000 | val =  0.144
#         (temp.water -4): pval = 0.00001 | val =  0.118
#         (PAR_sum  0): pval = 0.00186 | val =  0.084 | unoriented link
#         (GPP_raw -1): pval = 0.00392 | val =  0.078
#         (PAR_sum -2): pval = 0.00906 | val =  0.071
#         (temp.water -5): pval = 0.01956 | val = -0.063

#     Variable discharge has 6 link(s):
#         (discharge -1): pval = 0.00000 | val =  0.647
#         (GPP_raw  0): pval = 0.00000 | val = -0.204 | unoriented link
#         (PAR_sum  0): pval = 0.00000 | val = -0.179 | unoriented link
#         (discharge -2): pval = 0.00010 | val = -0.106
#         (ER_raw  0): pval = 0.00132 | val = -0.087 | unoriented link
#         (GPP_raw -1): pval = 0.00359 | val = -0.079

#     Variable PAR_sum has 9 link(s):
#         (PAR_sum -1): pval = 0.00000 | val =  0.298
#         (discharge  0): pval = 0.00000 | val = -0.179 | unoriented link
#         (GPP_raw  0): pval = 0.00000 | val =  0.144 | unoriented link
#         (PAR_sum -4): pval = 0.00000 | val =  0.132
#         (PAR_sum -3): pval = 0.00012 | val =  0.104
#         (temp.water  0): pval = 0.00186 | val =  0.084 | unoriented link
#         (PAR_sum -2): pval = 0.00515 | val =  0.076
#         (PAR_sum -5): pval = 0.00641 | val =  0.074
#         (temp.water -1): pval = 0.01043 | val = -0.069

# So, as demonstrated with the output above, the information regarding
# causality increases exponentially with the addition of added covariates.

# Adjust p_matrix to be more conservative.
q_matrix = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], tau_max=5, fdr_method='fdr_bh')
pcmci.print_significant_links(
        p_matrix = q_matrix,
        val_matrix = results['val_matrix'],
        alpha_level = 0.01)

# ## Significant links at alpha = 0.01:

#     Variable GPP_raw has 7 link(s):
#         (ER_raw  0): pval = 0.00000 | val =  0.249
#         (discharge  0): pval = 0.00000 | val = -0.204
#         (GPP_raw -1): pval = 0.00000 | val =  0.156
#         (PAR_sum  0): pval = 0.00000 | val =  0.144
#         (GPP_raw -2): pval = 0.00008 | val =  0.120
#         (GPP_raw -5): pval = 0.00067 | val =  0.106
#         (GPP_raw -3): pval = 0.00105 | val =  0.102

#     Variable ER_raw has 8 link(s):
#         (ER_raw -1): pval = 0.00000 | val =  0.375
#         (GPP_raw  0): pval = 0.00000 | val =  0.249
#         (ER_raw -2): pval = 0.00000 | val =  0.195
#         (discharge -1): pval = 0.00000 | val =  0.184
#         (ER_raw -3): pval = 0.00000 | val =  0.165
#         (GPP_raw -1): pval = 0.00003 | val = -0.126
#         (discharge -4): pval = 0.00196 | val = -0.098
#         (discharge  0): pval = 0.00132 | val = -0.087

#     Variable temp.water has 6 link(s):
#         (temp.water -1): pval = 0.00000 | val =  0.758
#         (temp.water -2): pval = 0.00000 | val = -0.316
#         (PAR_sum -1): pval = 0.00000 | val =  0.151
#         (temp.water -3): pval = 0.00000 | val =  0.144
#         (temp.water -4): pval = 0.00010 | val =  0.118
#         (PAR_sum  0): pval = 0.00186 | val =  0.084

#     Variable discharge has 5 link(s):
#         (discharge -1): pval = 0.00000 | val =  0.647
#         (GPP_raw  0): pval = 0.00000 | val = -0.204
#         (PAR_sum  0): pval = 0.00000 | val = -0.179
#         (discharge -2): pval = 0.00071 | val = -0.106
#         (ER_raw  0): pval = 0.00132 | val = -0.087

#     Variable PAR_sum has 6 link(s):
#         (PAR_sum -1): pval = 0.00000 | val =  0.298
#         (discharge  0): pval = 0.00000 | val = -0.179
#         (GPP_raw  0): pval = 0.00000 | val =  0.144
#         (PAR_sum -4): pval = 0.00001 | val =  0.132
#         (PAR_sum -3): pval = 0.00086 | val =  0.104
#         (temp.water  0): pval = 0.00186 | val =  0.084

graph = pcmci.get_graph_from_pmatrix(p_matrix=q_matrix, alpha_level=0.01, 
            tau_min=0, tau_max=5, selected_links=None)
results['graph'] = graph

# Process graph
# From the tutorial: "In the process graph, the node color denotes 
# the auto-MCI value and the link colors the cross-MCI value. If 
# links occur at multiple lags between two variables, the link color 
# denotes the strongest one and the label lists all significant lags 
# in order of their strength."

tp.plot_graph(
    val_matrix=results['val_matrix'],
    graph=results['graph'],
    var_names=var_names,
    link_colorbar_label='cross-MCI',
    node_colorbar_label='auto-MCI',
    ); plt.show()

# Time series graph
# Plot time series graph    
tp.plot_time_series_graph(
    figsize=(6, 4),
    val_matrix=results['val_matrix'],
    graph=results['graph'],
    var_names=var_names,
    link_colorbar_label='MCI',
    node_size = 0.01, # removes points
    ); plt.show()

# The online tutorial then goes on to describe how to deal with
# non-linear relationships, but I have chosen to stop here for the moment.

# Additional considerations, including latent variables acting as
# colliders in the process graphs, are discussed in this tutorial:
    # https://github.com/jakobrunge/tigramite/blob/master/tutorials/tigramite_tutorial_assumptions.ipynb

#...............

# additional code used to get the above working

# While I try to get my data working, let's try theirs for tutorial purposes.
np.random.seed(42)     # Fix random seed
links_coeffs = {0: [((0, -1), 0.7), ((1, -1), -0.8)],
                1: [((1, -1), 0.8), ((3, -1), 0.8)],
                2: [((2, -1), 0.5), ((1, -2), 0.5), ((3, -3), 0.6)],
                3: [((3, -1), 0.4)],
                }
T = 1000     # time series length
data, true_parents_neighbors = toys.var_process(links_coeffs, T=T)
T, N = data.shape

# Initialize dataframe object, specify time axis and variable names
var_names = [r'$X^0$', r'$X^1$', r'$X^2$', r'$X^3$']

dataframe_jk = pp.DataFrame(data, 
                         datatime = np.arange(len(data)), 
                         var_names=var_names)

# Plot the time series
tp.plot_timeseries(dataframe_jk); plt.show()
# This works, so make sure data is read in as an array.

#...............

# End of script.
