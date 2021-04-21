function [DI,MI,WI] = dca(msa)


    pseudocount_weight = 0.5; % relative weight of pseudo count   (.5)
    theta = 0.2;              % threshold for sequence id in reweighting (.2-.3)

    
    N = size(msa,2);
    M = size(msa,1);
    msa = double(msa);
    q = max(msa(:));
    
    [Pij_true,Pi_true, Meff]=Compute_True_Frequencies(msa,M,N,q,theta);
    
    [Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q);
    
    C = Compute_C(Pij,Pi,N,q);
    
    invC = inv(C);
    clear C;
    
    W_frob_store_wgc = zeros(N,N,21,21);
    for i=1:N
        for j=1:N % only fill in the upper triangle of W, because only this part is used (residue pairs where i<j) in compute_energies_large.m                    
                % get the block of the correlation matrix that corresponds to sites i and j
                W_frob = ReturnW_frob_gap(invC, i, j, q);		
                % change the gauge in this block to the zero-sum gauge
                W_frob_store_wgc(i,j,:,:) = frob_link(W_frob,q);        
        end
    end
    WI = FrobNorm(W_frob_store_wgc);
    [MI,DI] = Compute_Results(Pij, Pi, Pij_true, Pi_true, invC, N, q);

    
end
