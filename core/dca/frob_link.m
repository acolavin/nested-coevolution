function W = frob_link(W,q)

%perform gauge change to the zero-sum gauge in the block of interest
for i = 1:q
    W(:,i) = W(:,i) - mean(W(:,i));
end
for i = 1:q
    W(i,:) = W(i,:) - mean(W(i,:));
end 

end
