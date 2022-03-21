# CellReprogramming
## Synthetic Data Analysis
### Full Network and Sampled Network Comparison
Please find the relevant files in SyntheticDataAnalysis/FullNetwork_SampledNetwork_Comparison folder
1. Run PrepareAttractorSet.ipynb file (with R kernel). GA library is required.
Configurable parameters 
n_genes= number of genes
num_sub_bn_sets=number of samples per attractor set for  building sampled networks
num_sel_bits=number of selected nodes for  building sampled networks
num_bits=number of nodes

2. Run full_network_constructBNs.sh from command line, 
Configurable parameters
Attractors_folder=The folder where the BN construction algorithm  will find attractors
output_folder=The folder where results from BN construction algorithm be saved. 
num_bns= Number of times the BN construction algorithm will run.  
n_sets=Number of attractor sets

3. Run sampled_network_constructBNs.sh from command line, 
Configurable parameters
Attractors_folder=The folder where the BN construction algorithm  will find attractors
output_folder=The folder where results from BN construction algorithm be saved. 
num_bns= Number of times the BN construction algorithm will run. 
n_sets=Number of attractor sets

4.Run PrepareDataForPlots.jl from command line.  DelimitedFiles,DataFrames, FileIO, JLD2,Statistics are the required packages 
Configurable parameters
n_genes=number of genes
shift_result_folder, steady_state_correlation_folder, cell_state_probability_folder are the folders where the probablity shift from departure to destination states, correlation between steady state probability distribution of full and sampled networks, attractor probabilities will be stored. 
