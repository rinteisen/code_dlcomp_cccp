function candidate_set = Candidate_Decide(num_UE,num_TP,num_CH,path_loss,average_rate,B_max,power_max,noise)
% review : some questions remain
threshold_consider = 0.01;
averge_rate = average_rate/max(average_rate);
candidate_set = zeros(num_TP,num_UE+1,num_CH);
[~,BS_node] = max(power_max);
power_max(BS_node) = power_max(BS_node)/4; % why /4 ?
power_max = power_max/num_CH;

for j=1:num_CH
    candidate_count = zeros(num_TP,1);
    for i=1:num_UE
        [temp TP] = sort(path_loss(i,:,j).*power_max,'descend');
        candidate_count(TP(1)) = candidate_count(TP(1)) + 1;
        candidate_set(TP(1),candidate_count(TP(1)),j) = i;
        count = 2;
        while(count<=B_max)
            if temp(count) > temp(1)*threshold_consider*averge_rate(i)
                candidate_count(TP(count)) = candidate_count(TP(count)) + 1;
                candidate_set(TP(count),candidate_count(TP(count)),j) = i;
                count = count + 1;
            else
                break;
            end
        end
        for b=1:num_TP
            if candidate_count(b) == 0
                [~,max_UE] = max(log2(1+path_loss(:,b,j)*power_max(b)./noise')./average_rate);
                candidate_set(b,1,j) = max_UE;
            end
        end
    end
end
end