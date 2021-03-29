function pdbcoords = pdbcacoords(pdb)

    resids = [pdb.Model.Atom.resSeq]; 
    resnames = {pdb.Model.Atom.AtomName}; 
    pdbseq = mat2cell(repmat('A',1,max(resids)),1,ones(1,max(resids)));
    pdbres = zeros(1,length(pdbseq));
    pdbcoords = zeros(length(pdbseq),3);
    
    for l = 1:max(resids)
        curi = find(resids==l);
        curj = find(ismember(resnames,{'CA'}));
        joint = curi(ismember(curi, curj));
        if length(joint)==1
            X1 = pdb.Model.Atom(joint).X;
            Y1 = pdb.Model.Atom(joint).Y;
            Z1 = pdb.Model.Atom(joint).Z;
            pdbcoords(l,:) = [X1;Y1;Z1];
        end
    end

end

