function [serving_UE_index Power S I] = UE_Selection_RRH_v4(S,I,path_loss,candidate_UE_all,num_CH,BS_UE,power_channel,average_rate)
% review done
serving_UE_index = zeros(1,num_CH);
Power = zeros(1,num_CH);
for j = 1:num_CH
    candidate_UE = candidate_UE_all(1,1:find(candidate_UE_all(1,:,j)==0,1,'first')-1,j);
    if ~isempty(candidate_UE)
        H1 = path_loss(:,1,j);
        object_value = log2(1+power_channel(j)*H1(candidate_UE)./I(candidate_UE,j))./average_rate(candidate_UE);
        [max_object select_UE] = max(object_value);
        max_object = max_object + sum(log2(1+S(:,j)./(I(:,j)+power_channel(j)*H1))./average_rate);
        zero_object = sum(log2(1+S(:,j)./I(:,j))./average_rate);
        if zero_object > max_object
            % choose the UE with the largest first order deviation
            object_value_diff = H1(candidate_UE)./average_rate(candidate_UE)./I(candidate_UE,j);
            [~,select_UE] = max(object_value_diff);
            serving_UE_index(j) = candidate_UE(select_UE);
        else
            serving_UE_index(j) = candidate_UE(select_UE);
            Power(j) = power_channel(j);
            S(serving_UE_index(j),j) = S(serving_UE_index(j),j) + Power(j)*path_loss(serving_UE_index(j),1,j);
            path_loss(serving_UE_index(j),1,j) = 0;
            I(:,j) = I(:,j) + Power(j)*path_loss(:,1,j);
        end
    end
end
end