function pairwised = seqpdist_RS(msa)

[rows, columns] = size(msa);
pairwised = zeros(rows,rows, 'single');
for i = 1:rows
    row_repeat = repmat(msa(i,:),rows,1);
    matches = row_repeat ~= msa;
    p = sum(matches,2)/columns;
    p = min(p,.949);
    pairwised(i,:) = -19/20*log(1-p*20/19);
end