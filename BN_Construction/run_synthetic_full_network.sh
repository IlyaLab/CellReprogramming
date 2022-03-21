Attractors_folder="../SyntheticDataAnalysis/FullNetwork_SampledNetwork_Comparison/Full_AttractorSet"
output_folder="../SyntheticDataAnalysis/FullNetwork_SampledNetwork_Comparison/Full_Network/RESULTS"
num_bns=50
n_sets=100
if [ -d "$output_folder" ];
then
	rm -r $output_folder;
fi
for j in $(seq 1 1 $n_sets); do
	for i in $(seq 1 1 $num_bns); do
		mkdir -p ${output_folder}/$j/$i
	./program --att_info_file ${Attractors_folder}/attractor_info.txt --attractors_file ${Attractors_folder}/attractors$j.txt --output_directory ${output_folder}/$j/$i --num_genes 10
	done
done
