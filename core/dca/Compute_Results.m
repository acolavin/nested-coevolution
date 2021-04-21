function [MI_true_mat, DI_mf_pc_mat] = Compute_Results(Pij,Pi,Pij_true,Pi_true,invC, N,q)
% computes and prints the mutual and direct informations
    MI_true_mat = zeros(N);
    DI_mf_pc_mat = zeros(N);
    for i=1:(N-1)
        for j=(i+1):N
            % mutual information
            [MI_true,~,~] = calculate_mi(i,j,Pij_true,Pi_true,q);
            
            % direct information from mean-field
            W_mf = ReturnW(invC,i,j,q); 
            DI_mf_pc = bp_link(i,j,W_mf,Pi,q);
            MI_true_mat(i,j) = MI_true;
            DI_mf_pc_mat(i,j) = DI_mf_pc;
            MI_true_mat(j,i) = MI_true_mat(i,j);
            DI_mf_pc_mat(j,i) = DI_mf_pc_mat(i,j);
            
        end
    end
end