function nc_signal = extract_signal(protein, metric, correction, nc_cutoff)


raw_signal = protein.coev(1).(metric);

% Calculate mean nested coevolution
nc_null = zeros(length(protein.msa(1).seq));

if isinf(nc_cutoff)
    inds = find(arrayfun(@(x) isinf(x), [protein.msa.nc_cutoff]));
else
    inds = find([protein.msa.nc_cutoff]==nc_cutoff);
end

for ind = inds
   nc_null = nc_null + protein.coev(ind).(metric);
end

nc_null = nc_null/length(inds);

% Subtract correction from raw signal
if ~isinf(nc_cutoff)
    nc_signal = (raw_signal - nc_null);
else
    nc_signal = raw_signal;
end

if strcmp(correction, 'apc')
    nc_signal = nc_signal - apc(nc_signal);
elseif strcmp(correction, 'apc2')
    nc_signal = nc_signal - apc2(nc_signal);
elseif strcmp(correction, 'asc')
    nc_signal = nc_signal - asc(nc_signal);
elseif ~strcmp(correction, 'off')
    error('Method not implemented');
end


end