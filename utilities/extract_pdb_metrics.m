function [pdbseq,pdbres,pdbdist2] = extract_pdb_metrics(pdb)
    resids = [pdb.Model.Atom.resSeq]; 
    pdbseq = mat2cell(repmat('A',1,max(resids)),1,ones(1,max(resids)));
    pdbres = zeros(1,length(pdbseq));

    for l = 1:max(resids)
        curi = find(resids==l);

        if ~isempty(curi)
            X1 = [pdb.Model.Atom(curi).X];
            Y1 = [pdb.Model.Atom(curi).Y];
            Z1 = [pdb.Model.Atom(curi).Z];
            pdbseq{l} = aminolookup(pdb.Model.Atom(curi(1)).resName);
            pdbres(l) = l;
            for g = (l+1):length(resids)

                curj = find(resids==g);
                if ~isempty(curj)
                    X2 = [pdb.Model.Atom(curj).X];
                    Y2 = [pdb.Model.Atom(curj).Y];
                    Z2 = [pdb.Model.Atom(curj).Z];

                    curd = pdist2([X1;Y1;Z1]',[X2;Y2;Z2]');
                    pdbdist2(l,g) = min(curd(:));
                    pdbdist2(g,l) = pdbdist2(l,g);
                end
            end
        end
    %     disp(l)
    end