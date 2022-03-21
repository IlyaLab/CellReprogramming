using Pkg
Pkg.add("Random")
Pkg.add("FileIO")
Pkg.add("JLD2")
Pkg.add("StatsBase")
Pkg.add("HypothesisTests")
Pkg.add("Base")
#Pkg.add("JuliaInterpreter")
using FileIO, JLD2
using Random
using StatsBase
using HypothesisTests
using Base
#using JuliaInterpreter


mutable struct node 
    predictor_function::Array{Int8,2}
    prob::Array{Float64,1}    
end 


function NextState(prev_state, pred_fun_vector, maxp, num_genes)
    p=0.002
    y=falses(num_genes)
    x=prev_state
    gam=rand(num_genes).<p
    b=2 .^ collect(StepRange(maxp-1, Int8(-1), 0))
    ngenes=length(pred_fun_vector)
    if any(gam)
        y=xor.(x,gam) 
    else
      for i=1:ngenes  
        z=pred_fun_vector[i].prob  
        cdf = cumsum(z);
        u = rand(Float64,1);
        s= sum(cdf .< u) + 1;
        #s=sample(Random.GLOBAL_RNG, collect(StepRange(1, Int8(1), length(z))),  Weights(z))
        pr=pred_fun_vector[i].predictor_function
        m= pr[s,:]
        num_pred=m[2]
        preds=m[3:(num_pred+2)]
        ind=x[preds]'*b[end-num_pred+1:end]+1
        y[i]= m[maxp+2+ind]   
        end
    end
    return(y)
end


function FindAttractors(n_iter, n_burnin, parameters, n_genes, maxp)
        num_attractors=2
        inc_num=10000

        z=falses(n_iter, n_genes)
        state=bitrand(n_genes)
        for i=1:n_iter
            z[i,:]=NextState(state, parameters, maxp, n_genes)
            state=z[i,:] 
        end  
        b=2 .^ collect(StepRange(n_genes-1, Int8(-1), 0))
        new_z=z[n_burnin+1:n_iter,:]     
        kolm_test=true
        while kolm_test==true
           mid_point= trunc(Int, (n_iter- n_burnin)/2)
           first_part=new_z[1:mid_point,:]
           second_part=new_z[(mid_point+1):size(new_z,1),:]
           first=mapslices(xx ->xx'*b, first_part, dims=2)
           second=mapslices(xx ->xx'*b, second_part, dims=2)
           inds=collect(StepRange(1, Int8(100), size(first,1)))
           first_used=first[inds]
           second_used=second[inds]
           fs=size(first_used,1)
           first_e=zeros(Int64, 2^n_genes)
           els1=sort(unique(first_used))
           cts1=counts(first_used)
           first_e[els1.+1]=cts1[cts1.!=0]
           second_e=zeros(Int64, 2^n_genes)
           els2=sort(unique(second_used))
           cts2=counts(second_used)
           second_e[els2.+1]=cts2[cts2.!=0]
           kolm_test=maximum(abs.(first_e.-second_e)/fs)>1.358*sqrt(2/fs)
           a=falses(inc_num, n_genes)
           if kolm_test==true
           for i in 1:inc_num
                a[i,:]=NextState(state, parameters, maxp, n_genes)
                state=a[i,:]
            end
            new_z=[new_z;a]
            mid_point=mid_point+inc_num
            n_iter=n_iter+inc_num
         end
        end
       all_states= mapslices(xx ->xx'*b, new_z, dims=2)
       els=sort(unique(all_states))
       z=counts(all_states)
       cts=z[z.!= 0]
       cts=cts/sum(cts)
       inds=sortperm(cts, rev=true)
       att=reshape(els[inds[1:2]],1,:)
       att_prob=reshape(cts[inds[1:2]], 1,:)
       return(att,att_prob)
end

n_iter=100000
n_burnin=50000

n_genes=10
maxp=3
nBNsets=100
parameter_folder="Parameters"
i=0
full_atts=[]
full_att_probs=[]
sampled_atts=[]
sampled_att_probs=[]
for i=1:nBNsets
    print(i)
    parameters_full=FileIO.load(parameter_folder *"/full_parameters_"*string(i)*".jld2",  "parameters")
    att_full, att_prob_full=FindAttractors(n_iter, n_burnin, parameters_full, n_genes, maxp)
    push!(full_atts, att_full)
    push!(full_att_probs,att_prob_full)
    print("-")
    parameters_sampled=FileIO.load(parameter_folder *"/sampled_parameters_"*string(i)*".jld2",  "parameters")
    att_sampled, att_prob_sampled=FindAttractors(n_iter, n_burnin, parameters_sampled, n_genes, maxp)
    push!(sampled_atts, att_sampled)
    push!(sampled_att_probs,att_prob_sampled)
end

attractor_folder="Attractors"
mkpath(attractor_folder)
FileIO.save(attractor_folder *"/sampled_attractors_"*".jld2","atts",sampled_atts, "probs", sampled_att_probs)
FileIO.save(attractor_folder *"/full_attractors_"*".jld2","atts",full_atts, "probs", full_att_probs)
