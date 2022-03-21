using FileIO, JLD2
using DelimitedFiles
using DataFrames
n_genes=10
# For attractor probabilities figure
struct steady_state 
    pair_no::Int8
    intv_bit::Int8
    val::Int8
    state_probs::DataFrame
end 

bin_departure=   "1010101110"
bin_destination= "1001010011"

departure=parse(Int, bin_departure, base = 2)
destination=parse(Int, bin_destination, base = 2)

attractor_folder="PBN_Construction/Attractors"
attractors_sampled=FileIO.load(attractor_folder *"/sampled_attractors_"*".jld2",  "atts")
attractors_sampled_prob=FileIO.load(attractor_folder *"/sampled_attractors_"*".jld2",  "probs")

departure_sampled_prob=zeros(1, size(attractors_sampled_prob,1))
destination_sampled_prob=zeros(1, size(attractors_sampled_prob,1))

for i in 1:size(attractors_sampled,1)
    if attractors_sampled[i,1][1]==departure
        departure_sampled_prob[i]=attractors_sampled_prob[i,1][1]
        destination_sampled_prob[i]=attractors_sampled_prob[i,1][2]
    else
       departure_sampled_prob[i]=attractors_sampled_prob[i,1][2]
       destination_sampled_prob[i]=attractors_sampled_prob[i,1][1]
    end
end
intervention_folder="PBN_Construction/States_After_Intervention"
sampled_intervention=FileIO.load(intervention_folder *"/sampled_steady_state"*".jld2",  "state_probs")

neg_pos_0=[];neg_pos_1=[];pos_neg_0=[]; pos_neg_1=[]; neg_neg_0=[]; neg_neg_1=[]; pos_pos_0=[]; pos_pos_1=[]

neg_ind='0'; pos_ind='1'


rel_bits_neg_pos=intersect(findall(x->x==neg_ind, bin_departure), findall(y->y==pos_ind, bin_destination))
rel_bits_pos_neg=intersect(findall(x->x==pos_ind, bin_departure), findall(y->y==neg_ind, bin_destination))
rel_bits_neg_neg=intersect(findall(x->x==neg_ind, bin_departure), findall(y->y==neg_ind, bin_destination))
rel_bits_pos_pos=intersect(findall(x->x==pos_ind, bin_departure), findall(y->y==pos_ind, bin_destination))
 
bin_departure=string(departure, base=2, pad=n_genes)
bin_destination=string(destination, base=2, pad=n_genes)

for i in 1:size(sampled_intervention,1)
    
    intv_val=sampled_intervention[i].val
    intv_bit=sampled_intervention[i].intv_bit
    
    sampled_states=zeros(Float64, 2^(n_genes-1))
    sampled_states[sampled_intervention[i].state_probs[:,"States"].+1]=sampled_intervention[i].state_probs[:,"Probabilities"]
    
    initial_departure_sampled_prob=departure_sampled_prob[sampled_intervention[i].pair_no]
    initial_destination_sampled_prob=destination_sampled_prob[sampled_intervention[i].pair_no]
   
    departure_new=string(bin_departure[1:(intv_bit-1)], bin_departure[(intv_bit+1):end])
    departure_int= parse(Int, departure_new, base = 2)
    destination_new= string(bin_destination[1:(intv_bit-1)], bin_destination[(intv_bit+1):end])
    destination_int= parse(Int, destination_new, base = 2)
   
    shift=((sampled_states[destination_int+1]- initial_destination_sampled_prob)+(initial_departure_sampled_prob- sampled_states[departure_int+1]))/2
    
    if (intv_val==0)
        if (intv_bit in rel_bits_neg_pos)
            push!(neg_pos_0, shift)
        elseif (intv_bit in rel_bits_pos_neg)
            push!(pos_neg_0, shift)
        elseif (intv_bit in rel_bits_neg_neg)
            push!(neg_neg_0, shift)
        elseif (intv_bit in rel_bits_pos_pos)
            push!(pos_pos_0, shift)
        end
        
    elseif (intv_val==1)
        if (intv_bit in rel_bits_neg_pos)
            push!(neg_pos_1, shift)
        elseif (intv_bit in rel_bits_pos_neg)
            push!(pos_neg_1, shift)
        elseif (intv_bit in rel_bits_neg_neg)
            push!(neg_neg_1, shift)
        elseif (intv_bit in rel_bits_pos_pos)
            push!(pos_pos_1, shift)
        end
    end
end

shift_result_folder="Probability_Shift_CSVs"
mkpath(shift_result_folder)

writedlm(shift_result_folder*"/neg_pos_1.csv",  neg_pos_1, ",", header=false)
writedlm(shift_result_folder*"/neg_pos_0.csv",  neg_pos_0, ",", header=false)
writedlm(shift_result_folder*"/pos_neg_1.csv",  pos_neg_1, ",", header=false)
writedlm(shift_result_folder*"/pos_neg_0.csv",  pos_neg_0, ",", header=false)
writedlm(shift_result_folder*"/neg_neg_0.csv",  neg_neg_0, ",", header=false)
writedlm(shift_result_folder*"/neg_neg_1.csv",  neg_neg_1, ",", header=false)
writedlm(shift_result_folder*"/pos_pos_0.csv",  pos_pos_0, ",", header=false)
writedlm(shift_result_folder*"/pos_pos_1.csv",  pos_pos_1, ",", header=false)
