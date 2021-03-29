function [stpr, stds, varargout] = calculate_rand_stpr(protein, pdb_exclude, pdb_extent, pdb_contact, nb)

s = length(protein.msa(1).seq);
fake_signal = reshape(find(ones(s)),[s,s]);

residue_keep = protein.pdb.residue_keep;

pred_signal = get_predictions(fake_signal, residue_keep);
pred_signal_flat = squareform(pred_signal);

buffer_inds = find_buffer_inds(pred_signal, pdb_exclude);

% make sure buffer_inds takes into account gaps in pdb.
pred_signal_flat(buffer_inds) = [];


stprs = zeros(nb, pdb_extent);
pdbdist_compare_flat = squareform(protein.pdb.pdbdist_compare);
pdbdist_compare_flat(buffer_inds) = [];

for i = 1:nb
    [~, topinds_random] = sort(pred_signal_flat(randperm(length(pred_signal_flat))),'descend');
    aa = arrayfun(@(x) sum(pdbdist_compare_flat(topinds_random(1:x))<= pdb_contact), 1:pdb_extent);
    stprs(i, :) = arrayfun(@(x) sum(pdbdist_compare_flat(topinds_random(1:x))<= pdb_contact), 1:pdb_extent);
end

stpr = mean(stprs, 1);
stds = std(stprs, 1);

tp = stpr;
fp = (1:length(stpr))-stpr;

f1s = tp./(tp+fp);

if nargout == 3
   varargout{1} = f1s; 
end


end

function signal = get_predictions(signal, residue_keep)
    signal = signal(residue_keep, :);
    signal = signal(:, residue_keep);
    signal(find(diag(ones(1,size(signal,1))))) = 0;
end

function buffer_inds = find_buffer_inds(signal, buffer)
    remove = triu(ones(length(signal)),buffer);
    buffer_inds = find(~squareform(remove + remove'));


end