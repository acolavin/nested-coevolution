function FW = FrobNorm(W)
%Calculate the matrix of Frobenius norms of each block of W describing the direct couplings between 2 amino-acids
%W is a L*L*q*q matrix -> couplings between aas i and j are in W(i,j,:,:)

L=size(W,1);
q=size(W,3);
FW=zeros(L);

for i=1:L
    for j=1:L
        We=zeros(q);
        for k=1:q
            for l=1:q
                We(k,l)=W(i,j,k,l);
            end
        end            
        FW(i,j)=norm(We,'fro');
    end
end

end

