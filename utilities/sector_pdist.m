function [pdists, pdists2, V] = sector_pdist(protein, metric, correction, ...
    windows, zero, center, transpose, random_sectors)

unique_cutoffs = unique([protein.msa.nc_cutoff]);
m1 = cell2mat(arrayfun(@(u) extract_signal(protein, metric, correction, u), ...
    unique_cutoffs(windows), 'UniformOutput',0));

% zero:
if zero
    m1(m1<0)=0;
end

if center
    m1 = m1-repmat(mean(m1,1),size(m1,1),1);
end

% covariance:
if transpose
    m_covariance = m1*m1';
else
    m_covariance = m1'*m1;
end
% disp(size(m_covariance))

% eigen decomposition
[V,D] = eig(m_covariance);
V2 = zeros(size(m1));
s = size(m1,1);

if ~transpose
    % average across windows if large transpose.
    for w = 1:length(windows)
        V2 = V2 + V((s*(w-1)+1):(s*w),:);
    end
    V = V2/length(windows);
end

%% For each of the top ns eigenvectors, find the average pairwise distance
ns = 35;
extent = 100;
pdists = zeros(ns,2,extent-1);
for i = 1:ns
    if random_sectors
        v = V(:,randi(size(V,1),1));
    else
        v = V(:,end-i+1);
    end
    [~,sorted_inds] = sort(v);
%     sorted_inds = sorted_inds(ismember(sorted_inds, protein.pdb.residue_keep));
    
    %%% OLD
    % Added by Esha
    % Doesnt disregard indicies that are in the seq but not in the
    % structure. Outputs NaNs instead. Leads to long plateaus in the
    % pairwise distance vs top residues graphs
%     pdb_inds = [];
%     for s = 1:length(sorted_inds)
%         resid_temp = find(protein.pdb.residue_keep==sorted_inds(s));
%         if ~isempty(resid_temp) && resid_temp ~= 0
%             pdb_inds(s) = resid_temp;
%         else
%             pdb_inds(s) = NaN;
%         end
%     end
    
    %%% NEW 20200529
    pdb_inds = [];
    for s = 1:length(sorted_inds)
        resid_temp = protein.pdb.pdbind(find(protein.pdb.seqind==sorted_inds(s)));
        if ~isempty(resid_temp)
            pdb_inds(s) = resid_temp;
        else
            pdb_inds(s) = NaN;
        end
    end
    
%     pdb_inds = arrayfun(@(x) find(protein.pdb.residue_keep==x), sorted_inds);
    
    pdists(i, 1, 1) = 0;
    pdists(i, 2, 1) = 0;
    for j = 2:extent
       
        s1 = pdb_inds(1:j); % this tries negative values
        s2 = pdb_inds(end-j+1:end); % this tries positive values
        
        s1(isnan(s1)) = [];
        s2(isnan(s2)) = [];
        pdists(i, 1, j) = mean(pdist(protein.pdb.pdbaca_coords(s1,:)));
        pdists(i, 2, j) = mean(pdist(protein.pdb.pdbaca_coords(s2,:)));
        
%         pdists(i, 1, j-1) = mean(pdist(protein.pdb.pdbaca_coords(protein.pdb.resid(s1),:)));
%         pdists(i, 2, j-1) = mean(pdist(protein.pdb.pdbaca_coords(protein.pdb.resid(s2),:)));
%         %pdists(i, 1, j-1) = mean(pdist(protein.pdb.pdbaca_coords(s1,:)));
%         %pdists(i, 2, j-1) = mean(pdist(protein.pdb.pdbaca_coords(s2,:)));
        
        
    end
        
end
pdists2 = NaN; %reshape(pdists,[ns*2, extent-1]);


end