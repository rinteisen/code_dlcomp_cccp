function [ttt power] = bound_update_iterative(x0,power_max,path_loss,noise,average_rate,num_UE,num_TP,num_CH,m_order,bound_mode,iteration_max,p_type)
% master layer
% p_type
% 0: TLPC method (proposed method)
% 1: PSO (not in this function)
% 2: Approximate convex (not in this function)
% 3: Complex RRH (RRH adjust power in the same manner as BS)
% 4: Water-filling (BS adjust power in the same manner as RRH)
% 5: Equal (channels share the unused power equally)
% 6: only water-filling (no numerical method)
% 7: Single channel (without cross-channel control)
% 9: Randomized algortihm

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
if p_type == 10
    p_type = 9;
end
while iteration < iteration_max
    iteration = iteration + 1;
    ttt(iteration)=HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
    if ttt(iteration)>v_opt
        v_opt = ttt(iteration);
        power = x0(num_UE+1,:,:);
        opt_time = opt_time + 1;
    end
    % master layer
    if p_type == 0 || p_type == 3 || p_type == 4 || p_type == 6 || p_type == 9
        if p_type == 6
            p_type_input = 4;
        elseif iteration == 1 && ~ismember(p_type,[9 10])
            p_type_input = 4;
        else
            p_type_input = p_type;
        end
        x0(num_UE+1,:,:) = cross_power_update(x0,num_UE,num_TP,num_CH,serving_UE_index,path_loss,noise,average_rate,BS_node,power_max,p_type_input);
    elseif p_type == 5
        for b=1:num_TP
            x0(num_UE+1,b,:) = x0(num_UE+1,b,:) + (power_max(b)-sum(x0(num_UE+1,b,:)))/num_CH;
        end
    end
    v_temp = HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
    for b=1:num_TP
        power_bound(1,b,:) = x0(num_UE+1,b,:) + (power_max(b)-sum(x0(num_UE+1,b,:)))/num_CH;
    end
    % end master layer

    if p_type == 9
        alpha_c = 0.8;
        for j=1:num_CH
            power_temp = power_update_v4(x0(:,:,j),num_UE,path_loss(:,:,j),noise,average_rate,power_bound(1,:,j));
            for b = 1:num_TP
                if rand(1) < alpha_c
                    x0(num_UE+1,b,j) = power_temp(1,b);
                end
            end
        end
    elseif p_type ~= 6
        for j=1:num_CH
            x0(num_UE+1,:,j) = power_update_v4(x0(:,:,j),num_UE,path_loss(:,:,j),noise,average_rate,power_bound(1,:,j));
        end
    end
    if v_temp>ttt(iteration)
        ttt(iteration) = v_temp;
        power = x0(num_UE+1,:,:);
    end
    if ttt(iteration)>v_opt
        v_opt = ttt(iteration);
        power = x0(num_UE+1,:,:);
        opt_time = opt_time + 1;
    end
    m_order = m_order*0.85;

    if abs(v_temp_pre-ttt(iteration))/ttt(iteration) < 0.0001
        break
    else
        v_temp_pre=ttt(iteration);
    end

end
if p_type == 9
    power = x0(num_UE+1,:,:);
end
end