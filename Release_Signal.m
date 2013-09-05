function [S I] = Release_Signal(num_CH,S,I,path_loss,Power,serving_UE_index)
% review done
for j=1:num_CH
    if serving_UE_index(j) ~= 0
        S(serving_UE_index(j),j) = S(serving_UE_index(j),j) - Power(j)*path_loss(serving_UE_index(j),1,j);
        path_loss(serving_UE_index(j),1,j) = 0;
        I(:,j) = I(:,j) - Power(j)*path_loss(:,1,j);
    else
        I(:,j) = I(:,j) - Power(j)*path_loss(:,1,j);
    end
end