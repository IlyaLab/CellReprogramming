# CellReprogramming
## Synthetic Data Analysis
### Full Network and Sampled Network Comparison
Please find the relevant files in SyntheticDataAnalysis/FullNetwork_SampledNetwork_Comparison folder.
1. Run PrepareAttractorSet.ipynb file (with R kernel). GA library is required.
Configurable parameters  are
n_genes= number of genes,
num_sub_bn_sets=number of samples per attractor set for  building sampled networks,
num_sel_bits=number of selected nodes for  building sampled networks,
num_bits=number of nodes.

2. Run full_network_constructBNs.sh from command line.
Configurable parameters are
Attractors_folder=The folder where the BN construction algorithm  will find attractors,
output_folder=The folder where results from BN construction algorithm be saved, 
num_bns= Number of times the BN construction algorithm will run, 
n_sets=Number of attractor sets.

3. Run sampled_network_constructBNs.sh from command line. 
Configurable parameters are
Attractors_folder=The folder where the BN construction algorithm  will find attractors,
output_folder=The folder where results from BN construction algorithm be saved, 
num_bns= Number of times the BN construction algorithm will run,
n_sets=Number of attractor sets

4. Run PrepareDataForPlots.jl from command line. DelimitedFiles,DataFrames, FileIO, JLD2,Statistics are the required packages.  
Configurable parameters are
n_genes=number of genes,
shift_result_folder, steady_state_correlation_folder, cell_state_probability_folder are the folders where the probablity shift from departure to destination states, correlation between steady state probability distribution of full and sampled networks, attractor probabilities will be stored. 

5. Run Plots.ipynb notebook (R kernel). ggpubr,  ggplot2, gridExtra,ggpubr are the required libraries. figure_folder is the parameter for the folder where the figures will be saved. 

### Sampled Network Consistency Analysis
Please find the relevant files in SyntheticDataAnalysis/SampledNetworkConsistency folder.
1. Run Consistency_Analysis_PrepareAtts.ipynb notebook (R kernel). GA library is required. parent_sampled_folder= the folder where the sampled attractors will ve saved, num_sub_bn_sets= number of times the BN construction algorith will run.  num_sel_bits= number of bits that will be subsampled for BN construction. 

2. Run cons_sampled_network_constructBNs.sh from command line. The configurable parameters are Attractors_folder =the folder where the BN construction will read attractors, output_folder= where the results from BN construction algorithm will be saved, n_sets= the number of runs for the same  attractor set for consistency analysis, num_bns=Number of times the BN construction algorithm will run.

3. Run PrepareDataForPlot.jl from command line. The required packages are FileIO, JLD2, DelimitedFiles, DataFrames. attractor_folder= the folder where the attractors will be read from, intervention_folder= the folder where the intervention results will be read from, shift_result_folder=the folder the probability shift results will be written. 

4. Run Plots.ipynb notebook (R kernel). ggplot2, ggpubr are the required libraries. Configurable parameters are figure_folder=the folder where the plots will be saved, prob_shift_folder=the folder where probability shift values will be read from. 


