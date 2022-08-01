#using Pkg
#Pkg.add("Random")
#Pkg.add("FileIO")
#Pkg.add("JLD2")
#Pkg.add("StatsBase")
#Pkg.add("HypothesisTests")
#Pkg.add("FreqTables")
#Pkg.add("DataStructures")
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("Tables")

import Base.Threads.@threads
using FileIO, JLD2
using Random
using StatsBase
using DataStructures
using CSV
using DataFrames
using Base
#using HypothesisTests
mutable struct node 
    predictor_function::Array{Int8,2}
    prob::Array{Float64,1}    
end 


function NextState(prev_state, pred_fun_vector, maxp, num_genes)
    p=0.0001
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
        z=falses(n_iter, n_genes)
        state=bitrand(n_genes)
        for i=1:n_iter
            z[i,:]=NextState(state, parameters, maxp, n_genes)
            state=z[i,:] 
        end  
        b=2 .^ collect(StepRange(n_genes-1, Int8(-1), 0))
        new_z=z[(n_burnin+1):n_iter,:]     
        all_states=mapslices(xx ->join(xx.+0), new_z, dims=2) 
       	c = counter(all_states)
        return(c)
end

n_iter=300000
n_burnin=100000
nsets=500
parameter_main_folder=ARGS[1]
found_attractor_folder=ARGS[2]

n_genes=30
maxp=3

parameter_folder=parameter_main_folder
attractor_folder=found_attractor_folder
mkpath(attractor_folder)

@threads for i in 1: nsets
	parameters_sampled=FileIO.load(parameter_folder *"/parameters_"*string(i)*".jld2",  "parameters")
	att_sampled=FindAttractors(n_iter, n_burnin, parameters_sampled, n_genes, maxp)
	CSV.write(string(attractor_folder, "/", i , "_states.csv"),  att_sampled)
end