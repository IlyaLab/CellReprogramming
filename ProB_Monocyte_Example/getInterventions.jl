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
#Pkg.add("BenchmarkTools")
#using BenchmarkTools

#Pkg.add("Tables")
using FileIO, JLD2
using Random
using StatsBase
using DataStructures
using CSV
using DataFrames
import Base.Threads.@threads

#using HypothesisTests
using Base
mutable struct node 
    predictor_function::Array{Int8,2}
    prob::Array{Float64,1}    
end 


 function  NextStateIntervention(prev_state, pred_fun_vector, maxp, num_genes, ind, val)
       p=0.0001
       y=falses(num_genes)
        x=prev_state
        gam=rand(num_genes).<p
        b=2 .^ collect(StepRange(maxp-1, Int8(-1), 0))
        if any(gam)
            y=xor.(x,gam) 
            y[ind]=val

      else
       for i=1:num_genes
          if i==ind
               y[i]=val
               continue
            else
           z=pred_fun_vector[i].prob  
          cdf = cumsum(z);
               u = rand(Float64,1);
               s= sum(cdf .< u) + 1;
                pr=pred_fun_vector[i].predictor_function
                m= pr[s,:]
                num_pred=m[2]
                preds=m[3:(num_pred+2)]
                ix=x[preds]'*b[end-num_pred+1:end]+1
                y[i]= m[maxp+2+ix]  
              end
            end
         end
        return(y)
     end
   
function FindInterventionStates(n_iter, n_burnin, int_bit, val, parameters, n_genes, maxp)
        z=falses(n_iter, n_genes)
        state=bitrand(n_genes)
        for i=1:n_iter
            z[i,:]=NextStateIntervention(state, parameters, maxp, n_genes, int_bit, val)
            state=z[i,:] 
        end  
        b=2 .^ collect(StepRange(n_genes-1, Int8(-1), 0))
        new_z=z[(n_burnin+1):n_iter,:]     
        new_z=hcat(new_z[:,1:int_bit-1], new_z[:,int_bit+1:end])
        all_states=mapslices(xx ->join(xx.+0), new_z, dims=2) 
	   c = counter(all_states)
       return(c)
end

n_iter=300000
n_burnin=100000


parameter_folder =ARGS[1]
intv_folder =ARGS[2]

n_genes=30
maxp=3
mkpath(intv_folder)
nBNsets=500
 for i=1:nBNsets
    parameters_sampled=FileIO.load(parameter_folder *"/parameters_"*string(i)*".jld2",  "parameters")
     for int_bit=1: n_genes
	 	for val=0:1
    		intv_sampled=FindInterventionStates(n_iter, n_burnin, int_bit, val, parameters_sampled, n_genes, maxp)
    		CSV.write(string(intv_folder, "/",i ,"_", int_bit, "_", val,  "_states.csv"),  intv_sampled)
	   end
	end
end

