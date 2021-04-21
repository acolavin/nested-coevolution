function output_asc = asc(input_mat)

output_asc = zeros(size(input_mat));
input_mat(find(diag(ones(1,size(input_mat,1))))) = 0;
mi_mean = mean(input_mat(find(triu(ones(size(input_mat)),1))));
s = size(input_mat,1);
for i = 1:size(input_mat,1)
    for j = i:size(input_mat,2)
        
        mix = (sum(input_mat(i,:))-input_mat(i,j))/(s-1);
        mjy = (sum(input_mat(:,j))-input_mat(i,j))/(s-1);
        
        output_asc(i,j) = mix + mjy - mi_mean;
        output_asc(j,i) = output_asc(i,j);
    end
end


end