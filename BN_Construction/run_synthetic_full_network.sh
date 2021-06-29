Attractors_folder="../Synthetic_Data_Analysis/Full_Network/Full_AttractorSet"
output_folder="../Synthetic_Data_Analysis/Full_Network/RESULTS"
num_bns=49
n_sets=99
if [ -d "$output_folder" ];
then
	rm -r $output_folder;
fi
for j in $(seq 0 1 $n_sets); do
	for i in $(seq 0 1 $num_bns); do
		mkdir -p ${output_folder}/$j/$i
	./program --att_info_file ${Attractors_folder}/$j/attractor_info$i.txt --attractors_file ${Attractors_folder}/$j/attractors$i.txt --output_directory ${output_folder}/$j/$i --num_genes 10
	done
done
