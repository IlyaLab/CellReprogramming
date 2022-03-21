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
num_bns= Number of Boolean networks that the BN construction algorithm will build. 
