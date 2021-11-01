%% Demo of Nested Coevolution
% for structural contiguity analysis
clear all; close all; clc;

%% Loading
% load the protein "struct" generated elsewhere, such as generated in the
% contact_prediction_demo.m script.
load('generated_data/NP_417717-10000-aligned.fasta-20181028.mat')

%% Parameters

pdb_exclude = 5; % residues (buffer)
pdb_extent = 500; % residues (extent)
pdb_contact = 5; % angstroms (angstroms)

pdb_code = '4CZE';%'3H8A'; %'4OBY';%'3H8A'; %'1WVN'; % '1WVN'
chain ='A';
unique_cutoffs = unique([protein.msa.nc_cutoff]);

%% Match protein to structure
protein.pdb = load_pdb_struct(pdb_code, chain, protein.msa(1).seq);

%% Calculate background model for pairwise distance
nb = 2000;
extent = 100;
pdists_rand = zeros(nb,extent-1);
pdb_inds = 1:length(protein.pdb.pdbaca_coords);
for b = 1:nb
    for j = 2:extent

        s1 = pdb_inds(randperm(length(pdb_inds),j));
        pdists_rand(b,j-1) = mean(pdist(protein.pdb.pdbaca_coords(s1,:)));

    end
end

%% Get 5th and 95th percentile of background model
q025 = arrayfun(@(x) quantile(pdists_rand(:,x), .025), 1:(extent-1));
q975 = arrayfun(@(x) quantile(pdists_rand(:,x), .975), 1:(extent-1));

%% Calculate nested and full coevolution

metric='nje';
unique_cutoffs = unique([protein.msa.nc_cutoff]);

windows = 1:(length(unique_cutoffs));
window_nc = 1:3;
[pdists_zero2, pdists_zero, V_nc] = sector_pdist(protein, metric, 'apc', ...
    windows(window_nc), true, false, true, false);
[pdists_full2, pdists_full, V_fc] = sector_pdist(protein, metric, 'apc', ...
    windows(end), true, false, true, false);

%% Compare sectors
figure('Color', [1,1,1])
imagesc(abs(V_nc(:,end:-1:end-20)'*V_fc(:,end:-1:end-20)))
xlabel('FC Vectors');
ylabel('NC Vectors');
% notice index nc2xfc16 light up.
%% Compare sectors
% and use this to decide cutoffs and signs 
nc_sector_index = 2;
fc_sector_index = 16;
figure('Color', [1,1,1]);
subplot(2,1,1); hold on;
plot(V_nc(:,end-nc_sector_index+1), 'DisplayName', 'NC');
plot(V_fc(:,end-fc_sector_index+1), 'DisplayName', 'FC');
xlabel('Position along protein');
ylabel('Signal');
title('MreB');
%% Without picking cutoffs, what is the average distance of top residues?
subplot(2,1,2); hold on;
fill([1:extent-1, (extent-1):-1:1],[q025, q975(end:-1:1)],[.7,.7,.7], ...
    'DisplayName','95th Percentile Chance')
% The second index in these plot commands pick out the negatively- or positively-sorted residue indices (1 or 2 respectively).
% These are chosen based on the previous graph.
plot(squeeze(pdists_zero2(nc_sector_index,1,:)), 'Color',[0,0,1], 'DisplayName','NC');
plot(squeeze(pdists_full2(fc_sector_index,2,:)), 'Color',[1,0,0], 'DisplayName','FC');
ylabel('Mean pairwise distance');
xlabel('Number of top positions');
legend()

% the second index of pdists_full2 determines whether to test the positive (2) or
% negative (1) side of the eigenvector. 

%% Get sector residues for further analysis (e.g. mapping onto a structure)
nc_sector_resids = protein.pdb.resid(cell2mat(arrayfun(@(x) find(protein.pdb.residue_keep==x), ...
    find(V_nc(:,end-nc_sector_index+1)>.02), 'UniformOutput',0)'))';
fc_sector_resids = protein.pdb.resid(cell2mat(arrayfun(@(x) find(protein.pdb.residue_keep==x), ...
    find(V_fc(:,end-fc_sector_index+1)<-.05), 'UniformOutput',0)'))';
