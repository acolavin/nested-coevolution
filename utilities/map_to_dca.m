function msasort_dca = map_to_dca(msa_input)

if ~isinteger(msa_input)
    msasort_dca = aa2int(msa_input);
else
    msasort_dca = msa_input;
end
msasort_dca = msasort_dca + 1;
msasort_dca(msasort_dca>21) = 22;
msasort_dca(msasort_dca==22) = 1;

end