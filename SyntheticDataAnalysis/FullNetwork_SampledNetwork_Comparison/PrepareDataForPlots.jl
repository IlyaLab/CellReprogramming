#using Pkg
#Pkg.add("CSV")
using DelimitedFiles
using DataFrames
using FileIO, JLD2
using Statistics

n_genes=10
# For attractor probabilities figure
struct steady_state 
    pair_no::Int8
    intv_bit::Int8
    val::Int8
    state_probs::DataFrame
end 

#original_attractors = DataFrame(CSV.File("100_cell_state_pairs.txt",header=false))>table
original_attractors=readdlm("100_cell_state_pairs.txt", Int64)
orig_att_mat=original_attractors
departure=original_attractors[:,1]
destination=original_attractors[:,2]

attractor_folder="PBN_Construction/Attractors"

attractors_full=FileIO.load(attractor_folder *"/full_attractors_"*".jld2",  "atts")
attractors_sampled=FileIO.load(attractor_folder *"/sampled_attractors_"*".jld2",  "atts")
attractors_full_prob=FileIO.load(attractor_folder *"/full_attractors_"*".jld2",  "probs")
attractors_sampled_prob=FileIO.load(attractor_folder *"/sampled_attractors_"*".jld2",  "probs")

departure_full_prob=zeros(size(attractors_full_prob,1))
destination_full_prob=zeros(size(attractors_full_prob,1))

departure_sampled_prob=zeros(1, size(attractors_sampled_prob,1))
destination_sampled_prob=zeros(1, size(attractors_sampled_prob,1))

#for i in 1:size(attractors_full,1) # check whether all attractors are found correctly
 #  print(sort(orig_att_mat[i,:])==sort(vec(attractors_full[i,:][1])))
 #  print("\n")
 #  print(sort(orig_att_mat[i,:])==sort(vec(attractors_sampled[i,:][1])))
 #  print("\n")
#end
# For cell state probability distribution plot

for i in 1:size(attractors_full,1)
    if attractors_full[i,1][1]==departure[i]
        departure_full_prob[i]=attractors_full_prob[i,1][1]
        destination_full_prob[i]=attractors_full_prob[i,1][2]
    else
       departure_full_prob[i]=attractors_full_prob[i,1][2]
       destination_full_prob[i]=attractors_full_prob[i,1][1]
    end
end

for i in 1:size(attractors_sampled,1)
    if attractors_sampled[i,1][1]==departure[i]
        departure_sampled_prob[i]=attractors_sampled_prob[i,1][1]
        destination_sampled_prob[i]=attractors_sampled_prob[i,1][2]
    else
       departure_sampled_prob[i]=attractors_sampled_prob[i,1][2]
       destination_sampled_prob[i]=attractors_sampled_prob[i,1][1]
    end
end

a=0
# For steady state correlation  plot
intervention_folder="PBN_Construction/States_After_Intervention"
full_intervention=FileIO.load(intervention_folder *"/full_steady_state"*".jld2",  "state_probs")
sampled_intervention=FileIO.load(intervention_folder *"/sampled_steady_state"*".jld2",  "state_probs")

neg_pos_1_full=[]
neg_pos_1_sampled=[]
pos_neg_0_full=[]
pos_neg_0_sampled=[]

all_cor=[]
att_cor=[]

x=n_genes*2

#b=2 .^ collect(StepRange(n_genes-2 , Int8(-1), 0))
for i in 1:size(full_intervention,1)
    if mod(i,x)==1      
        global pair_spec_full_states=[]
        global pair_spec_sampled_states=[]
     end   
   
    #bin_departure=digits(departure[full_intervention[i].pair_no], base = 2, pad = n_genes)
    #bin_destination= digits(destination[full_intervention[i].pair_no], base = 2, pad = n_genes)
    
    #diff=bin_destination - bin_departure
    #source_bits=findall(x->x==(-1), diff)
    #target_bits=findall(x->x==0, diff)
    #print(departure[full_intervention[i].pair_no])
    
    full_states=zeros(Float64, 2^(n_genes-1))
    sampled_states=zeros(Float64, 2^(n_genes-1))
    full_states[full_intervention[i].state_probs[:,"States"].+1]=full_intervention[i].state_probs[:,"Probabilities"]
    sampled_states[sampled_intervention[i].state_probs[:,"States"].+1]=sampled_intervention[i].state_probs[:,"Probabilities"]
    
    initial_departure_full_prob=departure_full_prob[full_intervention[i].pair_no]
    initial_destination_full_prob=destination_full_prob[full_intervention[i].pair_no]

    initial_departure_sampled_prob=departure_sampled_prob[sampled_intervention[i].pair_no]
    initial_destination_sampled_prob=destination_sampled_prob[sampled_intervention[i].pair_no]
    
    bin_departure=string(departure[full_intervention[i].pair_no], base=2, pad=n_genes)
    bin_destination=string(destination[full_intervention[i].pair_no], base=2, pad=n_genes)
 
    neg_ind='0'
    pos_ind='1'
    rel_bits_pos=intersect(findall(x->x==neg_ind, bin_departure), findall(y->y==pos_ind, bin_destination))
    rel_bits_neg=intersect(findall(x->x==pos_ind, bin_departure), findall(y->y==neg_ind, bin_destination))

    intv_val=full_intervention[i].val
    intv_bit=full_intervention[i].intv_bit
    
   if ((intv_val==1) & (intv_bit in rel_bits_pos))
        departure_new=string(bin_departure[1:(intv_bit-1)], bin_departure[(intv_bit+1):end])
        departure_int= parse(Int, departure_new, base = 2)
        destination_new= string(bin_destination[1:(intv_bit-1)], bin_destination[(intv_bit+1):end])
        destination_int= parse(Int, destination_new, base = 2)
        
        push!(all_cor, cor(full_states, sampled_states)) 
        
        push!(neg_pos_1_full, ((full_states[destination_int+1]- initial_destination_full_prob)+(initial_departure_full_prob-    full_states[departure_int+1]))/2)
        push!(neg_pos_1_sampled,((sampled_states[destination_int+1]- initial_destination_sampled_prob)+(initial_departure_sampled_prob- sampled_states[departure_int+1]))/2)
        
        append!(pair_spec_full_states, full_states)
        append!(pair_spec_sampled_states, sampled_states)
  
    elseif ((intv_val==0) & (intv_bit in rel_bits_neg))
        departure_new=string(bin_departure[1:(intv_bit-1)], bin_departure[(intv_bit+1):end])
        departure_int= parse(Int, departure_new, base = 2)
        destination_new= string(bin_destination[1:(intv_bit-1)], bin_destination[(intv_bit+1):end])
        destination_int= parse(Int, destination_new, base = 2)
        
        push!(all_cor, cor(full_states, sampled_states)) 

        push!(pos_neg_0_full, ((full_states[destination_int+1]- initial_destination_full_prob)+(initial_departure_full_prob -    full_states[departure_int+1]))/2)
        push!(pos_neg_0_sampled, ((sampled_states[destination_int+1]- initial_destination_sampled_prob)+(initial_departure_sampled_prob -    sampled_states[departure_int+1]))/2)

        append!(pair_spec_full_states, full_states)
        append!(pair_spec_sampled_states, sampled_states)
    end
    if mod(i,x)==0
         if size(pair_spec_full_states,1)==0
            global a+=1
        end
        if size(pair_spec_full_states,1)>0
           push!(att_cor, cor(pair_spec_full_states, pair_spec_sampled_states))
        end
    end
 end

shift_result_folder="Probability_Shift_CSVs"
mkpath(shift_result_folder)

steady_state_correlation_folder="Correlation_CSVs"
mkpath(steady_state_correlation_folder)

cell_state_probability_folder="CellState_Probability_CSVs"
mkpath(cell_state_probability_folder)
## write results into files.


writedlm(shift_result_folder*"/neg_pos_1_full.csv", neg_pos_1_full, ",", header=false)
writedlm(shift_result_folder*"/neg_pos_1_sampled.csv", neg_pos_1_sampled, ",", header=false)
writedlm(shift_result_folder*"/pos_neg_0_full.csv", pos_neg_0_full, ",", header=false)
writedlm(shift_result_folder*"/pos_neg_0_sampled.csv", pos_neg_0_sampled, ",", header=false)

writedlm(steady_state_correlation_folder*"/all_correlation.csv", all_cor, ",", header=false)
writedlm(steady_state_correlation_folder*"/attractor_grouped_correlation.csv", att_cor, ",", header=false)

writedlm(cell_state_probability_folder*"/full_network_att_probs.csv", attractors_full_prob, ",", header=false)
writedlm(cell_state_probability_folder*"/sampled_network_att_probs.csv", attractors_sampled_prob, ",", header=false)


