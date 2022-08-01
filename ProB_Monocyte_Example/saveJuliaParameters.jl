# Save Parameters
#using Pkg
#Pkg.add("FileIO")
#Pkg.add("JLD2")
#Pkg.add("DelimitedFiles")
#Pkg.add("BenchmarkTools")

using BenchmarkTools
using DelimitedFiles
using FileIO, JLD2
using Random
import Base.Threads.@threads

function Function_Probability(A)
    U_A=unique(A, dims=1)
    func_count=zeros(Int64, size(U_A,1))
    for u=1:size(U_A,1)
        for d=1:size(A,1)
           if  A[d,:]==U_A[u,:]
                func_count[u]=func_count[u]+1
           end
        end
    end  
    prop=func_count/sum(func_count)
    return (U_A, prop)
end 

mutable struct node 
    predictor_function::Array{Int8,2}
    prob::Array{Float64,1}    
end 

function ComputeParameters(f, BNs_folder, Attractors_folder,  n_inner, n_genes, n)
    maxp=3
    a = node[]
    k=100 # this is from c file
    result_folder=BNs_folder
    data_gen=readdlm(string(Attractors_folder,"/", f, "/ind_matrix.txt"))
    real_ind= data_gen # .+1
    real_ind=convert(Array{Int8}, real_ind)
    mps=2^maxp
    LLEN= mps+n
    alls=zeros(Int8, n_inner*k*n, maxp+2+mps)
    for in_ind=1:n_inner
        folder=string(result_folder, "/",(f),"/" ,(in_ind)) #-1
        string1 = string(folder, "/", "bns_", 0, ".txt")
        bn=readdlm(string1,'\t',  '\n')
        for L = 1:k
          no_predictors = zeros(Int8,n);
            for j=1:n  
               while bn[LLEN*(L-1)+j+mps,no_predictors[j]+1]!=-10
                    no_predictors[j]= no_predictors[j]+1;
                    if no_predictors[j]== maxp
                    	break;
                    end
              end   
           end
           predictors=bn[LLEN*(L-1)+1+mps:LLEN*(L-1)+ mps+n,1:maxp] 
           func=bn[LLEN*(L-1)+1:LLEN*(L-1)+mps, 1:n]
           newf=fill(-1, (2^maxp, n))
           genes=real_ind[in_ind,:]
           for gn=1:n
           	num_pred=no_predictors[gn]
                for prd=1:num_pred
                	predictors[gn,prd]= real_ind[in_ind, predictors[gn, prd]]
                end
                ind=sortperm(predictors[gn,1:num_pred])
                predictors[gn, 1:num_pred]= predictors[gn,ind]
                x=bitstring.(collect(StepRange(0, Int8(1), 2^num_pred-1)))
                all_l=length(x[1])
                rel=SubString.(x, all_l-num_pred+1, all_l)
                x=collect.(rel)
                zz=transpose(Int8.(reduce(hcat,x)).-48)
                Base.permutecols!!(zz,ind)
                x=[2 .^ collect(StepRange(num_pred-1, Int8(-1), 0))]
                newf[zz*reduce(hcat, x).+1, gn]=func[1:2^num_pred, gn] 
             end
                all=hcat(reshape(genes, :, 1), reshape(no_predictors, :, 1), predictors,transpose(newf))
                alls[(in_ind-1)*k*n+(L-1)*n+1:(in_ind-1)*k*n+(L-1)*n+n, :]=all
              end
            end 
            for i=1: n_genes
                ss=alls[alls[:,1].==i,:]
                xx, yy=Function_Probability(ss)
                push!(a, node(xx, yy))
            end
            return(a)
        end


BNs_folder=ARGS[1]
Attractors_folder=ARGS[2]
parameters_folder=ARGS[3]

n_genes=30
n=10
n_inner=500
nsets=500
mkpath(parameters_folder)

@btime @threads for i=1:nsets
	params=ComputeParameters(i, BNs_folder, Attractors_folder,  n_inner, n_genes, n)
	FileIO.save(parameters_folder *"/parameters_"*string(i)*".jld2","parameters", params)
end

