Attractors_folder="Sampled_Network_Attractors"
output_folder="RESULTS"
mkdir $output_folder
num_bns=500
n_sets=500

for j in $(seq 1 1 $n_sets); do
	for i in $(seq 1 1 $num_bns); do
        mkdir -p ${output_folder}/$j/$i
	   ./../BN_Construction/program --att_info_file ${Attractors_folder}/$j/attractor_info.txt --attractors_file ${Attractors_folder}/$j/attractors$i.txt --output_directory ${output_folder}/$j/$i --num_genes 10
	done
done
