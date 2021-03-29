function boot_msa = calculate_bootstrap_msa_struct(msa_struct, nc_cutoff)
    mat = msa_struct.aa2int;
    pairwised = msa_struct.pwdist;
    curmat = calculate_bootstrap_msa(mat, nc_cutoff, pairwised);
    curmat_dca = map_to_dca(curmat);
    boot_msa.aa2int = curmat;
    boot_msa.dcaint = curmat_dca;
    boot_msa.nc_cutoff = nc_cutoff;
    
    boot_fields = fieldnames(boot_msa);
    
    for field = fieldnames(msa_struct)'
        
        if ~ismember(field, boot_fields)
            boot_msa.(field{1}) = NaN;
        end
        
    end
end