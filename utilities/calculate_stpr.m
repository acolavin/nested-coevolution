function [stpr, pred_signal, varargout] = calculate_stpr(pdb_struct, signal, pdb_exclude, pdb_extent, pdb_contact)


residue_keep = pdb_struct.residue_keep;

pred_signal = get_predictions(signal, residue_keep);
pred_signal_flat = squareform(pred_signal);

buffer_inds = find_buffer_inds(pred_signal, pdb_exclude);

% make sure buffer_inds takes into account gaps in pdb.
pred_signal_flat(buffer_inds) = [];

[~,topinds_signal] = sort(pred_signal_flat,'descend');

pdbdist_compare_flat = squareform(pdb_struct.pdbdist_compare);
pdbdist_compare_flat(buffer_inds) = [];
stpr = arrayfun(@(x) sum(pdbdist_compare_flat(topinds_signal(1:x))<= pdb_contact), 1:pdb_extent);

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