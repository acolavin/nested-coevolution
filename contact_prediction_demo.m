%% Demo of Nested Coevolution
% for DCA-based structural contact analysis
%% FASTA files and PDB Structures
fasta_file = 'msas/KH-1_PF00013_full_1WVN.fasta';
pdb_code = '1WVN';
chain = 'A';

% nc windows:
nc_cutoffs = [0.1:.1:0.8, 10];
num_cutoffs = length(nc_cutoffs);

% msa construction
jc_max = 1.0;

% Base coevolution metric to correct:
coev_method = 'di'; % {'di', 'fn', 'mi', 'nje'};

% Other NC correction: off, apc or asc:
correction = 'apc';

% number bootstrap for NC
coev_bootstraps = 4; % previously nb

% pdb
pdb_exclude = 5; % previously buffer
pdb_extent = 50; % previously extent
pdb_contact = 5; % previously angstroms


%% Generate bootstrapped MSAs according to NC cutoffs 
disp('Calculating coevolution')
% Creating protein structure, with appropriate bootstrapping.
protein.msa_file = fasta_file;
protein.msa = fasta2msa_struct(protein.msa_file, jc_max);
for nc_cutoff_ind = 1:length(nc_cutoffs)
    for boot_num = 1:coev_bootstraps
        nc_cutoff = nc_cutoffs(nc_cutoff_ind);
        protein.msa(end+1) = calculate_bootstrap_msa_struct(protein.msa(1), nc_cutoff);
    end
end

%% Measure coevolution for all MSAs
for boot_ind = 1:length(protein.msa)
    %coev_results{end+1} = calculate_coev_struct(protein.msa(msa_ind), protein.msa(1).pwdist, coev_methods);
    if boot_ind == 1
        protein.coev = calculate_coev_struct(protein.msa(boot_ind), protein.msa(1).pwdist, {coev_method});
    else
        protein.coev(end+1) = calculate_coev_struct(protein.msa(boot_ind), protein.msa(1).pwdist, {coev_method});
    end
    fprintf('Completed: (%i/%i) %2.2f%%\n',boot_ind, length(protein.msa), (boot_ind/length(protein.msa)*100))
end            

%% Plot corrected coevolution:
correction_to_plot = 'off';
coev_method_to_plot = 'di';

% define colormap
n1 = 10;
n2 = 10;
mymap = [linspace(1,1,n1)', linspace(0,1,n1)', linspace(0,1,n1)'];
mymap2 = [linspace(1,0,n2)', linspace(1,0,n2)', linspace(1,0,n2)'];
mymap = [mymap; mymap2];
clim = [0, n2]/300;

figure; 
n = length(nc_cutoffs);
for i = 1:n
    subplot(1, n, i);
    imagesc(extract_signal(protein, coev_method_to_plot, correction_to_plot, nc_cutoffs(i)));
    title(sprintf('Corrected Signal NC=%1.2f', nc_cutoffs(i)))
    axis equal; axis tight;
    caxis(clim);
    colormap(mymap2);
end
%% Extract sector signal
aveMetric = cell2mat(arrayfun(@(x) extract_signal(protein, coev_method, correction, x), nc_cutoffs, 'UniformOutput',0));
[vecs_corrected,score,latent,~,explained,~] = pca(aveMetric', 'Economy', false);
[vecs_raw,~,~,~,~,~] = pca(extract_signal(protein, coev_method, correction, Inf), 'Economy', false);

raw_ind = 3;
corrected_ind = 2;

figure;
subplot(1,2,1); hold on;
imagesc(vecs_corrected(:,1:30)' * vecs_raw(:, 1:30))
axis tight; axis equal
ylabel('Corrected Index')
xlabel('Raw Index')
title('Dot Product between vectors')

subplot(1,2,2); hold on;
plot(vecs_raw(:, raw_ind))
plot(vecs_corrected(:, corrected_ind))

%% load PDB
pdb_struct = load_pdb_struct(pdb_code, chain, protein.msa(1).seq);

%% Calculate structural contact prediction
base_metric = extract_signal(protein, coev_method, correction, Inf);

stprBASE = calculate_stpr(pdb_struct, base_metric, pdb_exclude, pdb_extent, pdb_contact);
[stprRAND, stpr_rand_std] = calculate_rand_stpr(pdb_struct, base_metric, pdb_exclude, pdb_extent, pdb_contact, pdb_extent);
for cutoff_ind = 1:num_cutoffs
    ave_metric = extract_signal(protein, coev_method, correction, nc_cutoffs(cutoff_ind));
    stprNC{cutoff_ind} = calculate_stpr(pdb_struct, ave_metric, pdb_exclude, pdb_extent, pdb_contact);
end
stpr.(correction).(coev_method).BASE = stprBASE;
stpr.(correction).(coev_method).RAND = stprRAND;
stpr.(correction).(coev_method).RANDSTD = stpr_rand_std;
stpr.(correction).(coev_method).NC = stprNC;

%% Plot
correction_to_plot = 'apc';

figure; hold on;
plot(1:pdb_extent, stpr.(correction_to_plot).(coev_method_to_plot).BASE, 'DisplayName', 'Raw Signal', 'LineWidth', 2, 'Color', [0, 0, 0]);
stpr_rand_mean = stpr.(correction_to_plot).(coev_method_to_plot).RAND;
stpr_rand_std = stpr.(correction_to_plot).(coev_method_to_plot).RANDSTD;
fill([1:pdb_extent, pdb_extent:-1:1], [stpr_rand_mean+2*stpr_rand_std, max(stpr_rand_mean(end:-1:1)-2*stpr_rand_std(end:-1:1), zeros(1,pdb_extent))], [.8,.8, .8], 'DisplayName', 'Rand+2*Std'); 
fill([1:pdb_extent, pdb_extent:-1:1], [stpr_rand_mean+stpr_rand_std, max(stpr_rand_mean(end:-1:1)-stpr_rand_std(end:-1:1), zeros(1,pdb_extent))], [.5,.5, .5], 'DisplayName', 'Rand+Std'); 

colors = jet(length(nc_cutoffs));
for i = 1:length(nc_cutoffs)
    plot(1:pdb_extent, stpr.(correction_to_plot).(coev_method_to_plot).NC{i}, 'DisplayName', sprintf('Corrected Signal NC=%1.2f', nc_cutoffs(i)), 'Color', colors(i,:));
end
plot([0, pdb_extent], [0, pdb_extent], 'k--', 'DisplayName', 'Perfect Recall')
legend('Location', 'northwest')