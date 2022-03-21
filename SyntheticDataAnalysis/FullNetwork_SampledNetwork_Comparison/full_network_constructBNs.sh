Attractors_folder="Full_Network_Attractors"
output_folder="Full_Network_BNs"
num_bns=50
n_sets=100
if [ -d "$output_folder" ];
then
	rm -r $output_folder;
fi
for j in $(seq 1 1 $n_sets); do
	for i in $(seq 1 1 $num_bns); do
		mkdir -p ${output_folder}/$j/$i
	   ./../../BN_Construction/program --att_info_file ${Attractors_folder}/attractor_info.txt --attractors_file ${Attractors_folder}/attractors$j.txt --output_directory ${output_folder}/$j/$i --num_genes 10
	done
done
