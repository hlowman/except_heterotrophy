# PCMCI Python script
# Heili Lowman
# April 11, 2022

# Using the tutorial found at:
# github.com/jokobrunge/tigramite/blob/master/tutorials/tigramite_tutorial_basics.ipynb

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
# the following is still throwing an error
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

# Check to make sure all the parts are there.
ts_dataframe.var_names # it's not empty
ts_dataframe.values # it has variables to call
ts_dataframe.datatime # and it is indexed

#...............

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
# Hmmmm, ok so this works. So, make sure data is read in as an array.

#...............

# Implement linear partial correlation.
parcorr = ParCorr(significance = 'analytic')

# Initialize the PCMCI method.
pcmci = PCMCI(
    dataframe=ts_dataframe,
    cond_ind_test=parcorr,
    verbosity=1)

# Plot lagged correlations
# Changed tau to 10, because I don't think it'll be much more than 2
correlations = pcmci.get_lagged_dependencies(tau_max=10, val_only=True)['val_matrix']
lag_func_matrix = tp.plot_lagfuncs(val_matrix=correlations); plt.show()

# Use scatterplots to check and see if dependencies are truly linear.
# scatter_lags sets the lag to use for every pair of variables to be their maximal absolute value from calculations above.
scatter_lags = np.argmax(np.abs(correlations), axis=2)
tp.plot_scatterplots(dataframe=ts_dataframe, add_scatterplot_args={'scatter_lags':scatter_lags}); plt.show()

# Just going to keep tau_max at 10 for now, but this coudl certainly be changed.
# Set pc_alpha to None so the PCMCI will choose the significance level in the condition-selection step itself.
# and set alpha_level to 0.05 to be the threshold for significance in plotting (less conservative than theirs).
pcmci.verbosity = 1
results = pcmci.run_pcmci(tau_max=10, pc_alpha=None, alpha_level=0.05)

##
## Step 1: PC1 algorithm with lagged conditions
##

#Parameters:
#independence test = par_corr
#tau_min = 1
#tau_max = 10
#pc_alpha = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
#max_conds_dim = None
#max_combinations = 1

## Resulting lagged parent (super)sets:

#     Variable GPP_raw has 13 link(s):
#     [pc_alpha = 0.5]
#         (GPP_raw -1): max_pval = 0.00000, min_val =  0.153
#         (GPP_raw -2): max_pval = 0.00001, min_val =  0.119
#         (GPP_raw -5): max_pval = 0.00018, min_val =  0.102
#         (GPP_raw -3): max_pval = 0.00053, min_val =  0.094
#         (discharge -1): max_pval = 0.00111, min_val = -0.088
#         (GPP_raw -4): max_pval = 0.01252, min_val =  0.068
#         (GPP_raw -10): max_pval = 0.22051, min_val =  0.033
#         (PAR_sum -3): max_pval = 0.35614, min_val =  0.025
#         (GPP_raw -6): max_pval = 0.37801, min_val =  0.024
#         (GPP_raw -7): max_pval = 0.39233, min_val =  0.023
#         (discharge -5): max_pval = 0.40763, min_val =  0.023
#         (PAR_sum -7): max_pval = 0.45905, min_val =  0.020
#         (ER_raw -6): max_pval = 0.48657, min_val =  0.019

#     Variable ER_raw has 11 link(s):
#     [pc_alpha = 0.3]
#         (ER_raw -1): max_pval = 0.00000, min_val =  0.375
#         (ER_raw -2): max_pval = 0.00000, min_val =  0.178
#         (discharge -1): max_pval = 0.00000, min_val =  0.175
#         (ER_raw -3): max_pval = 0.00000, min_val =  0.142
#         (GPP_raw -1): max_pval = 0.00008, min_val = -0.107
#         (ER_raw -6): max_pval = 0.00138, min_val =  0.087
#         (ER_raw -8): max_pval = 0.10026, min_val =  0.045
#         (PAR_sum -6): max_pval = 0.10111, min_val = -0.045
#         (discharge -4): max_pval = 0.10142, min_val =  0.044
#         (GPP_raw -2): max_pval = 0.21606, min_val = -0.034
#         (GPP_raw -5): max_pval = 0.27665, min_val =  0.030

#     Variable temp.water has 9 link(s):
#     [pc_alpha = 0.2]
#         (temp.water -1): max_pval = 0.00000, min_val =  0.753
#         (temp.water -2): max_pval = 0.00000, min_val = -0.213
#         (PAR_sum -1): max_pval = 0.00000, min_val =  0.182
#         (temp.water -7): max_pval = 0.00027, min_val =  0.099
#         (temp.water -4): max_pval = 0.00049, min_val =  0.094
#         (temp.water -10): max_pval = 0.00110, min_val =  0.088
#         (GPP_raw -1): max_pval = 0.00129, min_val =  0.087
#         (PAR_sum -2): max_pval = 0.02964, min_val =  0.059
#         (temp.water -5): max_pval = 0.17882, min_val = -0.037

#     Variable discharge has 9 link(s):
#     [pc_alpha = 0.2]
#         (discharge -1): max_pval = 0.00000, min_val =  0.646
#         (GPP_raw -1): max_pval = 0.00256, min_val = -0.082
#         (temp.water -1): max_pval = 0.00602, min_val =  0.075
#         (discharge -2): max_pval = 0.01075, min_val = -0.069
#         (ER_raw -7): max_pval = 0.09857, min_val =  0.045
#         (temp.water -2): max_pval = 0.10608, min_val =  0.044
#         (temp.water -9): max_pval = 0.11333, min_val = -0.043
#         (PAR_sum -7): max_pval = 0.13370, min_val =  0.041
#         (discharge -4): max_pval = 0.16111, min_val =  0.038

#     Variable PAR_sum has 10 link(s):
#     [pc_alpha = 0.2]
#         (PAR_sum -1): max_pval = 0.00000, min_val =  0.279
#         (PAR_sum -4): max_pval = 0.00000, min_val =  0.130
#         (PAR_sum -7): max_pval = 0.00003, min_val =  0.113
#         (PAR_sum -3): max_pval = 0.00033, min_val =  0.097
#         (PAR_sum -5): max_pval = 0.00526, min_val =  0.076
#         (PAR_sum -10): max_pval = 0.00687, min_val =  0.073
#         (PAR_sum -2): max_pval = 0.01244, min_val =  0.068
#         (PAR_sum -8): max_pval = 0.08371, min_val =  0.047
#         (PAR_sum -9): max_pval = 0.13341, min_val =  0.041
#         (discharge -1): max_pval = 0.16307, min_val =  0.038

# ##
# ## Step 2: MCI algorithm
# ##

# Parameters:

# independence test = par_corr
# tau_min = 0
# tau_max = 10
# max_conds_py = None
# max_conds_px = None

# ## Significant links at alpha = 0.05:

#     Variable GPP_raw has 10 link(s):
#         (ER_raw  0): pval = 0.00000 | val =  0.249 | unoriented link
#         (discharge  0): pval = 0.00000 | val = -0.199 | unoriented link
#         (PAR_sum  0): pval = 0.00000 | val =  0.157 | unoriented link
#         (GPP_raw -1): pval = 0.00000 | val =  0.152
#         (GPP_raw -2): pval = 0.00001 | val =  0.119
#         (GPP_raw -5): pval = 0.00017 | val =  0.102
#         (GPP_raw -3): pval = 0.00061 | val =  0.093
#         (discharge -1): pval = 0.00327 | val = -0.080
#         (GPP_raw -4): pval = 0.01040 | val =  0.070
#         (temp.water -1): pval = 0.03177 | val = -0.059

#     Variable ER_raw has 15 link(s):
#         (ER_raw -1): pval = 0.00000 | val =  0.373
#         (GPP_raw  0): pval = 0.00000 | val =  0.249 | unoriented link
#         (ER_raw -2): pval = 0.00000 | val =  0.191
#         (discharge -1): pval = 0.00000 | val =  0.188
#         (ER_raw -3): pval = 0.00000 | val =  0.174
#         (GPP_raw -1): pval = 0.00000 | val = -0.128
#         (discharge -4): pval = 0.00084 | val = -0.091
#         (ER_raw -6): pval = 0.00103 | val =  0.089
#         (discharge  0): pval = 0.00177 | val = -0.085 | unoriented link
#         (GPP_raw -2): pval = 0.00655 | val = -0.074
#         (PAR_sum -6): pval = 0.01408 | val = -0.067
#         (ER_raw -8): pval = 0.01771 | val =  0.065
#         (discharge -2): pval = 0.02210 | val = -0.062
#         (GPP_raw -6): pval = 0.03424 | val =  0.058
#         (discharge -3): pval = 0.04186 | val =  0.055

#     Variable temp.water has 10 link(s):
#         (temp.water -1): pval = 0.00000 | val =  0.748
#         (temp.water -2): pval = 0.00000 | val = -0.320
#         (PAR_sum -1): pval = 0.00000 | val =  0.143
#         (temp.water -3): pval = 0.00000 | val =  0.136
#         (temp.water -4): pval = 0.00024 | val =  0.100
#         (temp.water -7): pval = 0.00045 | val =  0.095
#         (GPP_raw -1): pval = 0.00162 | val =  0.086
#         (PAR_sum  0): pval = 0.00243 | val =  0.083 | unoriented link
#         (PAR_sum -2): pval = 0.00319 | val =  0.080
#         (temp.water -5): pval = 0.01280 | val = -0.068

#     Variable discharge has 6 link(s):
#         (discharge -1): pval = 0.00000 | val =  0.646
#         (GPP_raw  0): pval = 0.00000 | val = -0.199 | unoriented link
#         (PAR_sum  0): pval = 0.00000 | val = -0.184 | unoriented link
#         (discharge -2): pval = 0.00003 | val = -0.113
#         (ER_raw  0): pval = 0.00177 | val = -0.085 | unoriented link
#         (GPP_raw -1): pval = 0.00324 | val = -0.080

#     Variable PAR_sum has 14 link(s):
#         (PAR_sum -1): pval = 0.00000 | val =  0.265
#         (discharge  0): pval = 0.00000 | val = -0.184 | unoriented link
#         (GPP_raw  0): pval = 0.00000 | val =  0.157 | unoriented link
#         (PAR_sum -4): pval = 0.00004 | val =  0.111
#         (PAR_sum -3): pval = 0.00186 | val =  0.085
#         (PAR_sum -7): pval = 0.00230 | val =  0.083
#         (temp.water  0): pval = 0.00243 | val =  0.083 | unoriented link
#         (discharge -10): pval = 0.00600 | val =  0.075
#         (temp.water -1): pval = 0.00874 | val = -0.071
#         (ER_raw -10): pval = 0.00951 | val = -0.071
#         (ER_raw -6): pval = 0.01881 | val =  0.064
#         (ER_raw -1): pval = 0.02557 | val =  0.061
#         (PAR_sum -5): pval = 0.02978 | val =  0.059
#         (PAR_sum -2): pval = 0.03289 | val =  0.058

# Not entirely sure what's happening here....
# but running in case it's needed below for the final plots.
q_matrix = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], tau_max=8, fdr_method='fdr_bh')
pcmci.print_significant_links(
        p_matrix = q_matrix,
        val_matrix = results['val_matrix'],
        alpha_level = 0.05)
graph = pcmci.get_graph_from_pmatrix(p_matrix=q_matrix, alpha_level=0.05, 
            tau_min=0, tau_max=10, selected_links=None)
results['graph'] = graph

# Process graph
# From the tutorial: "In the process graph, the node color denotes the auto-MCI value and the link colors the cross-MCI value. If links occur at multiple lags between two variables, the link color denotes the strongest one and the label lists all significant lags in order of their strength."

tp.plot_graph(
    val_matrix=results['val_matrix'],
    graph=results['graph'],
    var_names=var_names,
    link_colorbar_label='cross-MCI',
    node_colorbar_label='auto-MCI',
    ); plt.show()

# lol, ok.

# Time series graph
# Plot time series graph    
tp.plot_time_series_graph(
    figsize=(6, 4),
    val_matrix=results['val_matrix'],
    graph=results['graph'],
    var_names=var_names,
    link_colorbar_label='MCI',
    ); plt.show()

# ok, well that kind of worked :)

# End of script.
