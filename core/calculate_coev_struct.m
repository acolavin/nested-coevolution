function coev_struct = calculate_coev_struct(msa_struct, pairwised, varargin)

if nargin == 4
    metrics = varargin{1};
else
    metrics = {'fn'};
end

if sum(ismember(metrics, {'fn', 'mi', 'di'})) > 0
    [DI, MI, FN] = dca(msa_struct.dcaint);
end

if ismember('nje', metrics)
%     if ~isinf(NC_cutoff)
%         NJE = SCALEDtime2(msa_struct.aa2int, Inf, pairwised) - SCALEDtime2(msa_struct.aa2int, NC_cutoff, pairwised); 
%     else
    NJE = SCALEDtime2(msa_struct.aa2int, Inf, pairwised);
%     end
end

for metrici = 1:length(metrics)
    metric = metrics{metrici};
    if strcmp(metric,'di')
        coev_struct.(metric) = DI;
    elseif strcmp(metric, 'mi')
        coev_struct.(metric) = MI;
    elseif strcmp(metric, 'fn')
        coev_struct.(metric) = FN;
    elseif strcmp(metric,  'nje')
        coev_struct.(metric) = NJE;
    else
        disp(metric)
        error('Unrecognized metric.')
    end

end