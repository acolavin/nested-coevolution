function pdb_struct = load_pdb_struct(pdbcode, chain, seq)

% Load structure
websave([pdbcode '.pdb'], ['https://files.rcsb.org/download/' pdbcode '.pdb']);
pdb = pdbread([pdbcode '.pdb']);
% pdb = getpdb(pdbcode);
if length(pdb.Model) > 1
    pdb.Model = pdb.Model(1);
end
if ~strcmp(chain, 'off')
    pdb.Sequence = pdb.Sequence(find([pdb.Sequence.ChainID]==chain));
    pdb.Model.Atom = pdb.Model.Atom(find([pdb.Model.Atom.chainID]==chain));
end

% Get residue information
isca = cellfun(@(x) strcmp(x, 'CA'), {pdb.Model.Atom.AtomName});
resnum = [pdb.Model.Atom.AtomSerNo];
resnum = resnum(isca);
% Get sequence and pwdist between residues
[pdbseq, pdbresid, pdbdist2] = extract_pdb_metrics(pdb);
pdbacacoord = pdbcacoords(pdb);
missing_pdb = find(sum(pdbdist2)==0);

% Align against reference sequence
aligninfo = localalign(seq, cell2mat(pdbseq));

% Perform mapping between sequences
pdbind = cumsum(aligninfo.Alignment{1}(3,:)~='-')+aligninfo.Start(2)-1;
seqind = cumsum(aligninfo.Alignment{1}(1,:)~='-')+aligninfo.Start(1)-1;
alignment_inddrop = unique([find(diff(pdbind)==0)+1, find(diff(seqind)==0)+1]);
refkeepind = find(~ismember(1:length(pdbind),alignment_inddrop)); % relative to alignment
keepindpdb = pdbind(refkeepind);
keepindseq = seqind(refkeepind);
to_remove_more = find(ismember(keepindpdb,missing_pdb));

keepindseq(to_remove_more) = [];
keepindpdb(to_remove_more) = [];

pdbdist_compare = pdbdist2(keepindpdb,:);
pdbdist_compare = pdbdist_compare(:, keepindpdb);

% Determine which residues are present in both structure and alignment
residue_keep = keepindseq;

%pdbdist_compare(find(diag(ones(1,length(pdbdist_compare)))))=0;

% Create structure of ... structure
pdb_struct.pdbind = keepindpdb;
pdb_struct.seqind = keepindseq;
pdb_struct.pdbdist = pdbdist2;
pdb_struct.pdbseq = pdbseq;
pdb_struct.resid = pdbresid(keepindpdb);
pdb_struct.residue_keep = residue_keep;
pdb_struct.pdbdist_compare = pdbdist_compare;
pdb_struct.pdbaca_coords = pdbacacoord;
end