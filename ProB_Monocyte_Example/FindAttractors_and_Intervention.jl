import Pkg;
Pkg.add("DelimitedFiles")
Pkg.add("Random");
Pkg.add("StatsBase")

using Random
using StatsBase
using DelimitedFiles

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

function NextStateIntervention(prev_state, pred_fun_vector, maxp, num_genes, ind, val)
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
            # s=sample(Random.GLOBAL_RNG, collect(StepRange(1, Int8(1), length(z))),  Weights(z))
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

function FindAttractors(n_iter, n_burnin, parameters, n_genes, maxp)
        z=falses(n_iter, n_genes)
        state=bitrand(n_genes)
        for i=1:n_iter
            z[i,:]=NextState(state, parameters, maxp, n_genes)
            state=z[i,:] 
        end  
        new_z=z[n_burnin+1:n_iter,:]
        b=2 .^ collect(StepRange(n_genes-1, Int8(-1), 0))
        all_states= mapslices(xx ->xx'*b, new_z, dims=2)
        els=sort(unique(all_states))
        z=counts(all_states)
        cts=z[z.!= 0]/sum(z)
        inds=sortperm(cts, rev=true)
        att=reshape(els[inds[1:2]],1,:)
        att_prob=reshape(cts[inds[1:2]], 1,:)
        return(att,att_prob)
end

function FindProbabilityShift(int_bit, val, params, n_genes, source, target, source_initial_prob, target_initial_prob)    
    n_iter=300000
    n_burnin=100000
    maxp=3
   
    # rel_bits=collect(StepRange(1, Int8(1), n_genes)
	#if(length(rel_bits)>0)
	 #shift=zeros(length(rel_bits))
      #for tr=1:length(rel_bits)  
           ind=int_bit
           z=falses(n_iter, n_genes)
           state=bitrand(n_genes)
           for i=1:n_iter
               z[i,:]=NextStateIntervention(state, params, maxp, n_genes, ind, val)
               state=z[i,:] 
            end 
            new_z=z[n_burnin+1:n_iter,:]
            new_z=hcat(new_z[:,1:ind-1], new_z[:,ind+1:end]) 
            b=2 .^ collect(StepRange(n_genes-2, Int8(-1), 0)) 
            all_states= mapslices(xx ->xx'*b, new_z, dims=2)
            els=sort(unique(all_states))
            z=counts(all_states)
            cts=z[z.!= 0]/sum(z)
          
            source_new=string(source[1:ind-1], source[ind+1:end])
            target_new=string(target[1:ind-1], target[ind+1:end])

            source_int=parse(Int, source_new, base=2)
            target_int=parse(Int, target_new, base=2)

           
	       source_ind2= findfirst(x->x==source_int, els)
            target_ind2= findfirst(y->y==target_int, els)
    
            if source_ind2==nothing
                source_final=0
            else
                source_final=cts[source_ind2]
            end
    
            if target_ind2==nothing
                target_final=0
            else
                target_final=cts[target_ind2]
           end 
            #xx=sortperm(cts, rev=true)[1]
            shift=((source_initial_prob -source_final) + (target_final - target_initial_prob))/2
  # end
    return(shift)
end 

mutable struct node 
    predictor_function::Array{Int8,2}
    prob::Array{Float64,1}    
end 

function ComputeParameters(full_or_sampled, att_file,  result_in_folder)
    maxp=3
    alg=full_or_sampled
    a = node[]
    k=100
    if  alg=="sampled"
       #result_folder="RESULTS"
        n=10
        data_gen= readdlm(att_file)
        real_ind= data_gen
        real_ind=convert(Array{Int8}, real_ind)
        n_inner=500

    else
        #result_folder="full_network/RESULTS"
        n=10
        n_inner=50
    end
    mps=2^maxp
    LLEN= mps+n
    n_genes=30
        alls=zeros(Int8, n_inner*k*n, maxp+2+mps)
        for in_ind=1:n_inner
            folder=string(result_in_folder,  "/" , (in_ind))
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
             if alg=="sampled"
                   genes=real_ind[in_ind,:]
               else
                   genes=collect(StepRange(1, Int8(1), n_genes))
             end
             for gn=1:n
                 num_pred=no_predictors[gn]
                 if alg=="sampled"
                    for prd=1:num_pred
                        predictors[gn,prd]= real_ind[in_ind, predictors[gn, prd]]
                    end   
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




f=parse(Int64, ARGS[1])
att_file=ARGS[2]
att_ind_file=ARGS[3]
input_folder=ARGS[4]
output_folder=ARGS[5]


#attractors=convert(Array{Int64,2}, readdlm("Attractors.csv"))
attractors= readdlm(att_file, ',', Int, '\n')

n_genes=30

f=parse(Int64, ARGS[1])
s_set=ComputeParameters("sampled", att_ind_file, input_folder) # f is from command line

sampled_attractors, sampled_att_probs = FindAttractors(300000, 100000, s_set, n_genes, 3)

source_ind=findfirst(x->x==attractors[(f+1),1], sampled_attractors)
target_ind=findfirst(x->x==attractors[(f+1),2], sampled_attractors)

source_init_prob=sampled_att_probs[source_ind]
target_init_prob=sampled_att_probs[target_ind]

atts_to_write=[sampled_attractors[source_ind] sampled_attractors[target_ind]]
probs_to_write=[source_init_prob target_init_prob]

writedlm(output_folder *  "/found_attractors.csv", atts_to_write, ',')
writedlm(output_folder *   "/attractor_probability.csv", probs_to_write, ',')

source=string(sampled_attractors[source_ind], base=2, pad=n_genes)
target=string(sampled_attractors[target_ind], base=2, pad=n_genes)

int_result= zeros(Float64, n_genes, 2)
for interv_bit=1:n_genes
	for  val=0:1
		int_result[interv_bit,(val+1)]= FindProbabilityShift(interv_bit,val, s_set,  n_genes, source, target, source_init_prob, target_init_prob)
	end
end

writedlm(output_folder *  "/intervention_result.csv", int_result, ',')

