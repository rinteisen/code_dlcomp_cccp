function x = Transfer_Binary(num_UE,num_TP,num_CH,serving_UE_index,Power)
% review down
% transfer UE index into variabel x
x = [];
for j=1:num_CH
    x_j = zeros(num_UE+1,num_TP);
    x_j(num_UE+1,:) = Power(:,j)';
    for b=1:num_TP
        if serving_UE_index(b,j) ~= 0
            x_j(serving_UE_index(b,j),b) = 1;
        end
    end
    x = cat(3,x,x_j);
end
end