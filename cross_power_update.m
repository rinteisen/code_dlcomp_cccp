function power = cross_power_update(x0,num_UE,num_TP,num_CH,serving_UE_index,path_loss,noise,average_rate,BS_node,power_max,p_type)
% review down
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

iteration_max = 100;
v_opt = HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
resource = x0(1:num_UE,:,:);
power = reshape(x0(num_UE+1,:,:),num_TP,num_CH);
power_i = power;
S = zeros(num_UE,num_CH);
I = noise'*ones(1,num_CH);
for j=1:num_CH
    I(:,j) = I(:,j) + path_loss(:,:,j)*power_i(:,j);
    for b=1:num_TP
        if serving_UE_index(b,j) ~= 0
            serving_UE = serving_UE_index(b,j);
            S(serving_UE,j) = S(serving_UE,j) + power_i(b,j)*path_loss(serving_UE,b,j);
        end
    end
    I(:,j) = I(:,j) - S(:,j);
end
v_temp_pre = 0;
TP_flag = zeros(1,num_TP);
RRH_node = 1:num_TP;
if p_type == 3
    BS_node = 1:num_TP;
    RRH_node = [];
elseif p_type == 4 || p_type == 9
    BS_node = [];
else
    RRH_node(BS_node) = [];
end
if p_type == 9
    iteration_max_2 = 1;
    iteration_max = 1;
else
    iteration_max_2 = num_CH*5;
end
m_order = 0.5;
end_flag = 0;
if power_max(b)<sum(power_i(b,:))-1e-6
    
end
for count = 1:iteration_max
    if p_type == 4
        beta = 0;
    else
        beta = 0.8;
    end
    [~,random_priority] = sort(rand(1,num_TP-length(BS_node)));
    RRH_node = RRH_node(random_priority);
    TP_priority = [BS_node RRH_node];
    for b = TP_priority
        if p_type == 7
            m_order = m_order*0.9;
        end
        if ismember(b,BS_node)
            
            if TP_flag(b) == 0
                [S I] = Release_Signal(num_CH,S,I,path_loss(:,b,:),power_i(b,:),serving_UE_index(b,:));
                if p_type == 3
                    [power_i(b,:) S I TP_flag(b)] = Power_Update_BS_v4(S,I,serving_UE_index(b,:),path_loss,average_rate,num_UE,...
                        m_order,power_max(b),b,power_i(b,:)*beta,power_i(b,:),iteration_max_2);
                else
                    [power_i(b,:) S I TP_flag(b)] = Power_Update_BS_v4(S,I,serving_UE_index(b,:),path_loss,average_rate,num_UE,...
                        m_order,power_max(b),b,power_i(b,:)*beta,power_i(b,:),iteration_max_2);
                end
            else
%                 [power_i(b,:) S I] = Power_Update_RRH_v4(num_UE,power_max(b),S,I,path_loss(:,b,:),serving_UE_index(b,:),average_rate,power_i(b,:)*beta);
                TP_flag(b) = TP_flag(b)-1;
            end
        else
            [S I] = Release_Signal(num_CH,S,I,path_loss(:,b,:),power_i(b,:),serving_UE_index(b,:));
            [power_i(b,:) S I] = Power_Update_RRH_v4(num_UE,power_max(b),S,I,path_loss(:,b,:),serving_UE_index(b,:),average_rate,power_i(b,:)*beta);
        end
if power_max(b)<sum(power_i(b,:))-1e-6
    
end
if ~isempty(find(power_i<-1e-8, 1)) || ~isreal(power_i)
    a='cross_power_update'
    power_i
    b
end
        v_temp = HetNetfun_power(cat(1,resource,reshape(power_i,1,num_TP,num_CH)),num_UE,num_CH,noise,path_loss,average_rate);
        if v_temp>v_opt
            v_opt = v_temp;
            power(b,:) = power_i(b,:);
        end
    end
    v_temp = HetNetfun_power(cat(1,resource,reshape(power_i,1,num_TP,num_CH)),num_UE,num_CH,noise,path_loss,average_rate);
    if v_temp_pre>v_temp
        end_flag = end_flag + 1;
    end
    if end_flag == 3
        break;
    end
    if count~=1 && abs(v_temp_pre-v_temp)/v_temp < 0.0002
        if p_type == 7
            if count > iteration_max/5
                break
            end
        else
            break
        end
    else
        v_temp_pre=v_temp;
    end
    
end
power = reshape(power,1,num_TP,num_CH);
end