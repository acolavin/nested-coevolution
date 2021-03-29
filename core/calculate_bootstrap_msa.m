function boot_dca = calculate_bootstrap_msa(msasort_dca, cutoff, pairwised)
    boot_dca = msasort_dca;
    done = [];
    for s = 1:size(boot_dca,1)
        if length(find(done==s))==0
            curperm = find(pairwised(s,:)<=cutoff);

            for j = 1:size(boot_dca,2)

                boot_dca(curperm,j) = boot_dca(datasample(curperm,length(curperm)),j);
            end
            done = [done curperm];
        end
    end
end