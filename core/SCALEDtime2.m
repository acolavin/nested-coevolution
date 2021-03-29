function Cp_spec = SCALEDtime2(msa,phylo_frac,varargin)
% SCALED will return a pairwise coevolution matrix
% given an MSA, a fractional bootstrap frequency, and a phylogenetic window
% for the bootstrapping.
%   msa:    integer-representation of the multiple sequence alignment
%   nb:     number of times to run the bootstrapping
%   phylo_window_frac:  the fraction of the median phylogenetic distance
%                       between sequences to scramble.
%   varargin:   a fourth argument - a square matrix of pairwise distances 
%               between sequences - can be supplied to prevent the pairwise
%               distance from being recalculated every time. 
% Cp_raw is the unfiltered bootstrapped pairwise matrix
% Cp_spec is with speciation above phylo_frac removed.
%%
if ischar(msa)
    msa = aa2int(msa);
    msa(find(msa>=21))=0; % change gaps and wierd aa to zeros
end


if nargin==3
    pairwised = varargin{1};
    pairwised2 = num2cell(pairwised,2);
else
    fprintf('Getting pairwise distance information.\n');
    pairwised = squareform(seqpdist(int2aa(msa),'UseParallel',true));
    pairwised2 = num2cell(pairwised,2);
end


toswitch = cellfun(@(x) find(x<phylo_frac),pairwised2,'UniformOutput',0);

curtree = 0;
subtree = ones(1,length(toswitch))*-1;

for n = 1:length(toswitch)
    if subtree(n) == -1
        subtree(toswitch{n}) = curtree;
        curtree = curtree+1;
    end
end

%%

seqlength = size(msa,2);

Cp_spec = reshape(misc5(uint32(msa),uint32(subtree),uint32(curtree)),seqlength*[1,1]);


