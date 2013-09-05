function power_opt = power_update_v4(x0,num_UE,path_loss,noise,average_rate,power_bound)
round_max = 10;
trace_performance = zeros(1,round_max);
consider_TP = find(sum(x0(1:num_UE,:),1)>0);
power_opt = x0(num_UE+1,:);
num_TP = length(consider_TP);
H_i = zeros(num_TP,1);
serving_UE_index = zeros(1,num_TP);
v_opt=HetNetfun_power(x0,num_UE,1,noise,path_loss,average_rate);
v_ini = v_opt;

for b = 1:num_TP
    if sum(x0(num_UE+1,b,:)) > power_bound(1,b,:)+1e-6
        a = 'power_update_v4'
        x0(num_UE+1,b,:)
        power_bound(1,b,:)
    end
end

for b=consider_TP
    UE_index = find(x0(1:num_UE,b)==1);
    if ~isempty(UE_index)
        serving_UE_index(b) = UE_index;
    end
    H_i(b) = path_loss(UE_index,b);
end
if length(consider_TP) >= 2
    path_loss_i = path_loss(:,consider_TP);
    H = path_loss_i(:,:);
    resource = x0(1:num_UE,consider_TP);
    power_bound = power_bound(1,consider_TP);
    power = power_opt(1,consider_TP);
    average_rate = average_rate*ones(1,num_TP);
    serving_UE_index = serving_UE_index(consider_TP);
    H_i = H_i(consider_TP);
    power_first_round = max(power,power_bound);
    power_first_round = [power_first_round;power];
    power_first_round(2,1) = sqrt(max(power_bound(2:length(power_bound)))*power_bound(1));
    power_first_round = [power_first_round;max(power,power_bound)];
    power_first_round(3,1) = max(max(power_bound(2:length(power_bound)))*1.5,0);
    power_first_round = [power_first_round;max(power,power_bound)];
    power_first_round(4,1) = 0;
    power_first_round = [power_first_round;power];
    power_first_round = max(min(power_first_round,ones(size(power_first_round,1),1)*power_bound),0);
    power_second_round = zeros(size(power_first_round,1)+1,length(consider_TP));
    for iii = 1:size(power_first_round,1)
        power_i = power_first_round(iii,:);
        v_temp=HetNetfun_power(cat(1,resource,power_i),num_UE,1,noise,path_loss_i,average_rate(:,1));
        if v_temp>v_opt
            v_opt = v_temp;
            power = power_i;
        end
        for ii = 1:round_max
            S = (path_loss_i(:,:).*resource)*power_i'*ones(1,num_TP);
            I = (path_loss_i(:,:)*power_i' + noise' - S(:,1))*ones(1,num_TP);
            % sum_intf = xi between Eq.(47) and Eq.(48)
            sum_intf = sum(H.*S./average_rate./I./(S+I),1);
            sum_intf = sum_intf - (H_i.*S(serving_UE_index,1)./average_rate(serving_UE_index,1)./I(serving_UE_index,1)./(S(serving_UE_index,1)+I(serving_UE_index,1)))';
            if sum(sum_intf) == 0
                power_i = power_bound;
            else
                % power_i here: the first term in min of Eq.(48)
                power_i = ((H_i./sum_intf'./average_rate(serving_UE_index,1) - I(serving_UE_index,1) - S(serving_UE_index,1) - H_i.*power_i')./H_i)';
                power_i(power_i<0) = 0;
                power_i = min(power_i,power_bound);
            end
            if ~isempty(find(power_i<-1e-8,1))
                power_i
            end
            xtest = cat(1,resource,power_i);
            v_temp=HetNetfun_power(xtest,num_UE,1,noise,path_loss_i,average_rate(:,1));
            if v_temp>v_opt
                v_opt = v_temp;
                power = power_i;

            end
            trace_performance(ii) = v_temp;
            power_second_round(iii,:) = power_i;
            if ii~=1 && abs(v_temp_pre-v_temp)/v_temp < 0.0001
                break
            else
                v_temp_pre=v_temp;
            end
        end
    end
    power_opt(consider_TP) = power;
else
    power_opt(consider_TP) = power_bound(1,consider_TP);
end
end