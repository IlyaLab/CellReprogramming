attractors_folder="./Sampled_Network_Attractors"
BNfolder="./RESULTS"
parameters_folder="./Parameters"
julia saveJuliaParameters.jl  ${BNfolder} ${attractors_folder} ${parameters_folder} --threads 10

