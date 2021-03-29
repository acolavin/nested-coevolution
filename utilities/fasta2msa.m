function [msa, varargout] = fasta2msa(filenameout,varargin)
rawdata = fastaread(filenameout);
%% Convert MSA to numbers

% prune nonmatching MSAs is desired.
if nargin == 2
   if varargin{1}
       modeseq = mode(cellfun(@length,{rawdata.Sequence}));
       rawdata = rawdata(cellfun(@length,{rawdata.Sequence})==modeseq);
   end
end

numseq = length(rawdata);
seqlength = length(rawdata(1).Sequence);
msa = char(zeros(numseq,seqlength));

for curseq = 1:numseq
    % this line converts characters to their number representation
    msa(curseq,:) = rawdata(curseq).Sequence;
end

headers = cell(1,size(msa,1));
if nargout == 2
    for curseq = 1:numseq
        headers{curseq} = rawdata(curseq).Header; 
    end
    varargout{1} = headers;
end

end