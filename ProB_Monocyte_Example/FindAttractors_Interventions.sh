n_sets=500
attractor_file="Attractors.csv"
result_folder="RESULTS"
output_folder="Output"
attractor_folder="Sampled_Network_Attractors"
for j in $(seq 1 1 $n_sets); do
	julia  FindAttractors_and_Intervention.jl $j ${attractor_file}  ${attractor_folder}/$j/ind_matrix.txt ${result_folder}/$j  ${output_folder}/$j
done
