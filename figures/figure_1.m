%% Figure 1
% Alexandre Colavin, April 2018
% This figure demonstrates the changing landscape of apparent coevolution
% when the nested coevolution cutoff is varied.
clear all; close all;
%% Parameters

mat_file = '../generated_data/KH-1_PF00013_full_1WVN.fasta.mat';
pdb_code = '1WVN';
chain = 'A';

%% Initialize/load structure and load MSA
protein = load(mat_file);
protein = protein.protein;
%% define metrics
metric = 'di';
correction = 'off';


%% define colormap
n1 = 10;
n2 = 30;
mymap = [linspace(1,1,n1)', linspace(0,1,n1)', linspace(0,1,n1)'];
mymap2 = [linspace(1,0,n2)', linspace(1,0,n2)', linspace(1,0,n2)'];
mymap = [mymap; mymap2];
clim = [0, n2]/300;

%% Calculate raw signal
raw_signal = protein.coev(1).(metric);
if strcmp(correction, 'asc')
    raw_signal = raw_signal-asc(raw_signal);
elseif strcmp(correction, 'apc')
    raw_signal = raw_signal-apc(raw_signal);
end
unique_cutoffs = unique([protein.msa.nc_cutoff]);

%% Plot figure B/C/D
figure('Color', [1,1,1]);
cutoff = unique_cutoffs(5);
subplot(1,3,1);
imagesc(raw_signal);
subplot(1,3,2);
imagesc(raw_signal-extract_signal(protein, metric, correction, cutoff));
subplot(1,3,3);
imagesc(extract_signal(protein, metric, correction, cutoff));

for i = 1:3
    subplot(1,3,i);
    axis equal; axis tight;
    caxis(clim);
    colormap(mymap2);
end
%% Now lets plot all subwindows for figure E (and raw)

figure('Color',[1,1,1]);
for i = [1,2,3,4]
    subplot(1,5,i);
    imagesc(extract_signal(protein, metric, correction, unique_cutoffs(i*2+2)));
end


subplot(1,5,5);
imagesc(raw_signal);

for i = 1:5
    subplot(1,5,i);
    axis equal; axis tight;
    caxis(clim);
    colormap(mymap2);
end
%% Lets load PDB structure
protein.pdb = load_pdb_struct(pdb_code, chain, protein.msa(1).seq);


%% Now lets do F, showing how each window has different sTPR

% PDB parameters
pdb_exclude = 5; % residues
pdb_extent = 50; % residues
pdb_contact = 5; % angstroms
correction = 'apc';
metric='di';
raw_signal = protein.coev(1).(metric);
if strcmp(correction, 'asc')
    raw_signal = raw_signal-asc(raw_signal);
elseif strcmp(correction, 'apc')
    raw_signal = raw_signal-apc(raw_signal);
end
stpr_raw = calculate_stpr(protein.pdb, raw_signal, pdb_exclude, pdb_extent, pdb_contact);
[stpr_rand_mean, stpr_rand_std] = calculate_rand_stpr(protein.pdb, raw_signal, pdb_exclude, pdb_extent, pdb_contact, 500);
figure('Color', [1,1,1]);
hold on;

fill([1:pdb_extent, pdb_extent:-1:1], [stpr_rand_mean+stpr_rand_std, max(stpr_rand_mean(end:-1:1)-stpr_rand_std, zeros(1,pdb_extent))], [.5,.5, .5], 'DisplayName', 'Rand'); 
plot(1:pdb_extent, stpr_raw, 'DisplayName', 'Raw Signal', 'LineWidth', 2, 'Color', [0, 0, 0]); 

colors = jet(40);
for i = 1:14
    aa = extract_signal(protein, metric, correction, unique_cutoffs(i));
    stpr = calculate_stpr(protein.pdb, aa, pdb_exclude, pdb_extent, pdb_contact);
    plot(1:pdb_extent, stpr, 'DisplayName', sprintf('Corrected Signal NC=%1.2f', unique_cutoffs(i)), 'Color', colors(i,:));
    
end

plot([0, pdb_extent], [0, pdb_extent], '--')
