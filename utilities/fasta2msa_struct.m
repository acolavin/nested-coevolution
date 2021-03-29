function msa_struct = fasta2msa_struct(fileloc, jc_prune, varargin)

[raw_msa, headers] = fasta2msa(fileloc);

%% Define sequence as mode of MSA
if contains(headers(end), 'Reference sequence')
    seq = raw_msa(end,:);
else
    ref_loc = find(cellfun(@(x) contains(x,'ReferenceSequence'), headers));
    if length(ref_loc)==1
        seq = raw_msa(ref_loc,:);
    else
        seq = int2aa(mode(aa2int(raw_msa)));
    end
end

%% Only keep where mode is not gap
msap_raw = raw_msa(:, seq~='-');
msap = aa2int(msap_raw);
msap = unique(msap,'rows');
seq = seq(seq~='-');

%% Sort MSA and look at it.
[msasort,phylodist,~] = phylosort(msap_raw,seq);
num2keep = find(phylodist>=jc_prune, 1);

if isempty(num2keep)
    num2keep = size(msasort,1);
end

if nargin == 2
    msa_keep = msasort(1:num2keep, :);
else
    % Esha Atolia Sampling Mode!
    rng(1)
    k = varargin{1};
    if k == 1
        % keep everything up to cutoff
        msa_keep = msasort(1:num2keep,:);
    elseif k == 2
        % keep 
        num2keep = randi(num2keep,1,ceil(num2keep./10));
        msa_keep = msasort(num2keep, :);
    elseif k == 3
        rng(1);
        num2keep = randi(num2keep,1,ceil(num2keep./100));
        msa_keep = msasort(num2keep, :);
    elseif k == 4
        num2keep = ceil(num2keep./10);
        msa_keep = msasort(1:num2keep, :);
    elseif k == 5
        num2keep = ceil(num2keep./100);
        msa_keep = msasort(1:num2keep, :);
    end
end

[~,phylodist,~] = phylosort(msa_keep,seq);
num2keep = size(msa_keep,1);

%% Fill in the gaps
msa_keep2 = aa2int(msa_keep);
msa_keep2(find(msa_keep2>=21))=0;
pairwised = seqpdist_RS(msa_keep2);
%%
for s = 1:num2keep
    [~,curi] = sort(pairwised(:,s),'ascend');
    for p = 1:length(seq)
        if msa_keep2(s,p)==0
            msa_keep2(s,p) = msa_keep2(curi(find(msa_keep2(curi,p)~=0,1)),p);
        end
    end
end

%% Generate DCA-ready MSA by remapping residues
msasort_dca = map_to_dca(msa_keep2);

msa_struct.raw = raw_msa;
msa_struct.headers = headers;
msa_struct.seq = seq;
msa_struct.aa2int = msa_keep2;
msa_struct.dcaint = msasort_dca;
msa_struct.pwdist = pairwised;
msa_struct.nc_cutoff = Inf;
msa_struct.phylodist = phylodist;
end
