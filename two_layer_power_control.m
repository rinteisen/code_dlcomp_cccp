function [ttt power] = two_layer_power_control(x0,power_max,path_loss,noise,average_rate,num_UE,num_TP,num_CH,iteration_max)

% Transfer binary matrix into combination of index
BS_node = 1;
power_bound = x0(num_UE+1,:,:);
for b=1:num_TP
    if (power_max(b)-sum(x0(num_UE+1,b,:)))<-1e-6
    end
    power_bound(1,b,:) = power_bound(1,b,:) + (power_max(b)-sum(x0(num_UE+1,b,:)))/num_CH;
if ~isempty(find(power_bound(1,b,:)<0,1))
    reshape(power_bound(1,b,:),1,num_CH)
end
end
serving_UE_index = zeros(num_TP,num_CH);
for j=1:num_CH
    for b=1:num_TP
        serving_UE = find(x0(1:num_UE,b,j)==1);
        if ~isempty(serving_UE)
            serving_UE_index(b,j) = serving_UE;
        end
    end
end
iteration = 0;
power = x0(num_UE+1,:,:);
v_opt = HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
ttt = zeros(1,iteration_max);
v_temp_pre = v_opt;
opt_time = 0;

% Start to run interaction between two layers
while iteration < iteration_max
    iteration = iteration + 1;
    ttt(iteration)=HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
    if ttt(iteration)>v_opt
        v_opt = ttt(iteration);
        power = x0(num_UE+1,:,:);
        opt_time = opt_time + 1;
    end
    % Run the water-filling at first time (about 1% gain in practice...)
    if iteration == 1
        p_type_input = 4;
    else
        p_type_input = 0;
    end
    % Master layer
    x0(num_UE+1,:,:) = cross_power_update(x0,num_UE,num_TP,num_CH,serving_UE_index,path_loss,noise,average_rate,BS_node,power_max,p_type_input);
    v_temp = HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
    for b=1:num_TP
        power_bound(1,b,:) = x0(num_UE+1,b,:) + (power_max(b)-sum(x0(num_UE+1,b,:)))/num_CH;
    end
    
    % Slave layer
    for j=1:num_CH
        x0(num_UE+1,:,j) = power_update_v4(x0(:,:,j),num_UE,path_loss(:,:,j),noise,average_rate,power_bound(1,:,j));
    end

    % Record the best solution ever meet
    if v_temp>ttt(iteration)
        ttt(iteration) = v_temp;
        power = x0(num_UE+1,:,:);
    end
    if ttt(iteration)>v_opt
        v_opt = ttt(iteration);
        power = x0(num_UE+1,:,:);
        opt_time = opt_time + 1;
    end

    if abs(v_temp_pre-ttt(iteration))/ttt(iteration) < 0.001
        break
    else
        v_temp_pre=ttt(iteration);
    end

end
end