function [sortmsa,phylodist,idists] = phylosort(msa,refseq)



if ~ischar(msa)
    msa = int2aa(msa);
end

refseq = upper(refseq);

phylodist = zeros(1,size(msa,1));

for l=1:(size(msa,1)-1)
    phylodist(l) = seqpdist({refseq,msa(l,:)});
    
end

[phylodist,idists] = sort(phylodist,'ascend');

sortmsa = msa(idists,:);
end