%%
clear

nbns=100;
Res={};
folder='../Sampled_Intervention_SS/';
mkdir(folder);
addpath("../.")
for i=1:nbns
    indices_file=strcat('../Sampled_Attractor_Sets/', string(i-1), '/indices.txt');
    [varF, nf, nv, cij ,F]=GetPBNParameters(strcat('RESULTS/',string(i-1), '/'),indices_file );
    Res{i}{1}=varF;
    Res{i}{2}=nf;
    Res{i}{3}=nv;
    Res{i}{4}=cij;
    Res{i}{5}=F;
end
save(strcat(folder,'parameters.mat'),  'Res')

 

%%
 
 % now use BN/PBN package for PBN construction 
  
% y = pbnNextState(x,F,varF,nf,nv,cij,p)
% y = pbnNextState(x,F,varF,nf,nv,cij,p) - one step of a PBN
% This function performs one step of the PBN. The current state is x
% (binary vector of 1 x n, where n is # of genes) For other parameters, see
% e.g. pbnrnd.m

%function [F,varF,cij] = pbnRnd(n,nf,nv)

% [F,varF,cij] = pbnRnd(n,nf,nv) - random BN/PBN
% This function creates a random PBN with n nodes. The number of functions
% per node is defined in nf. Each function may have different number of
% variables, defined in nv. Function returns all the necessary parameters
% (excluding nf and nv), to define an independent PBN.
% INPUT:
% n     - The number of nodes in a random pbn.
% nf    - The number of functions per node. nf is a 1-by-n vector where
%         each element of nf is an integer larger than 0. If nf is a vector
%         of ones, then this function generates a random Boolean network.
% nv    - The number of variables in each Boolean function. The length of
%         nv must be equal to the sum of the elements in vector nf. Let cnf
%         be the cumulative sum of the number of functions in the network
%         (starting from the first node), i.e., cnf = cumsum(nf);. Then,
%         the number of variables of the functions for the first node are
%         nv(1:cnf(1)), for the second node nv(cnf(1)+1:cnf(2)), ..., and
%         for the last node nv(cnf(n-1)+1:cnf(n)).
% OUTPUT:
% F     - The functions for the random network. F has size
%         2^max(nv)-by-sum(nf). Truth tables of the functions for the first
%         node are defined in F(:,1:cnf(1)), for the second node
%         F(:,cnf(1)+1:cnf(2)), ..., and for the last node
%         F(:,cnf(n-1)+1:cnf(n)). Since the length of the truth tables may
%         vary between different functions, only the first 2^nv(i) bits are
%         relevant in the i:th column of F. Remaining "unnecessary" bits in
%         each column are set to -1. Let the f = F(;,i) be the i:th
%         function in F, and assume that it is a function of three
%         variables xi, xj, and xk (variables are defined in varF(1:3,i).
%         Then, f(1) defines the output for the input vector (000).
%         Correspondinly, f(2) defines the output for the input vector
%         (001), where xi = xj = 0 and xk = 1, i.e., where the third input
%         variable is equal to one. As another example, f(6) defines the
%         output for the input vector (101), where xi = 1, xj = 0, and
%         xk = 1.
% varF  - The variables of each function. varF is a max(nv)-by-sum(nf)
%         matrix. Variables of the functions for the first node are
%         varF(:,1:cnf(1)), for the second node varF(:,cnf(1)+1:cnf(1)),
%         ..., and for the last node varF(:,cnf(n-1)+1:cnf(n)). Since the
%         number of variables may vary between different functions, only
%         the first nv(i) elements are relevant in the i:th column of varF.
%         Remaining "unnecessary" elements in each column are set to -1.
% cij   - The selection probabilities for the functions. cij is a
%         max(nf)-by-n matrix. The selection probabilities of the functions
%         for the i:th node are cij(:,1:i). Since the number of functions
%         may vary between nodes, only the first nf(i) elements are
%         relevant in the i:th column of cij. Remaining "unnecessary"
%         elements in each column are set to -1.
%   
%y = pbnNextState(x,F,varF,nf,nv,cij,p)
%%
%y = pbnNextState(x,F,varF,nf,nv,cij,p)

      


before_intervention=zeros(1, 2^m)
for i=1:size(states,1)
    state_ind = bi2de(states(i,:),'left-msb')+1;
    before_intervention(state_ind)=state_probs(i);
end

%%%%% ---
att1=[1 1 1 1 1 0 0 0 1 0];
att2=[0 1 1 0 0 1 0 1 0 1];

int_zero=zeros(1,n_genes);
int_one=zeros(1, n_genes);

for ind=1:n_genes
    states_pert_zero=Old_PermanentIntervention(p,x,F,varF,nf,nv,cij,ind, 0, att1, att2);
    int_zero(ind)=states_pert_zero;
    states_pert_one=Old_PermanentIntervention(p,x,F,varF,nf,nv,cij,ind, 1,  att1, att2);
    int_one(ind)=states_pert_one;
end

    
    %%
int_res_neig=zeros(30,2)
int_res_neig2=zeros(30,2)

for ind=1:30
    for value=0:1
    states_pert_neig=PermanentIntervention_Neig(x,F,varF,nf,nv,cij,ind, value, mono_state,3);
    %int_res_neig(ind, value+1)=states_pert;
    int_res_neig2(ind, value+1)=states_pert_neig;
    end
end

% no zero values in any state probability
%[A,v] = pbnA(F,varF,nf,nv,cij,p)
%%
%size(unique(y, 'rows', 'stable'))

%%
%cnf = [0,cumsum(nf)];
% for gene=1:n_genes
 %   pmfrnd(cij(1:nf(gene),gene))
 %   k = cnf(i)+j
 %end
 %%
% attractors=['100001011100000',  '010101001101011', '101001101101000', '110001001101010']
  
%  function [funcs_, function_probs]  = RuleFunction(func)
%     [Mu,ia,ic] = unique(func, 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
%     h = accumarray(ic, 1);                              % Count Occurrences
%     maph = h(ic);                                       % Map Occurrences To ?ic? Values
%     Result = [func, maph]; 
%     B = unique(Result,'rows', 'stable');
%     %[A,B,C]=unique(preds, 'row', 'stable')
%     %unique( C)
%     funcs_=B(:, 1:(size(B,2)-1));
%     function_probs=B(:,size(B,2))/sum(B(:,size(B,2)));
%  end  

