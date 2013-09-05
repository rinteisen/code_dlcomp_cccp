function [Power_opt S_opt I_opt BS_flag] = Power_Update_BS_v4(S,I,serving_UE_index,path_loss,average_rate,num_UE,m_order_p,power_max,BS_node,Power_lower,Power,iteration_max)
% review : beta_2?
% gradient projection method
Power_opt = Power;
S_temp = S;
I_temp = I;
path_loss_temp = path_loss(:,BS_node,:);
consider_CH = find(serving_UE_index(1,:)~=0);
serving_UE_BS_c = serving_UE_index(consider_CH) + num_UE*(consider_CH-1);
S_temp(serving_UE_BS_c) = S(serving_UE_BS_c) + Power(consider_CH).*path_loss_temp(serving_UE_BS_c);
I_temp(:,consider_CH) = I(:,consider_CH) + (ones(num_UE,1)*Power(consider_CH)).*reshape(path_loss_temp(:,1,consider_CH),num_UE,length(consider_CH));
I_temp(serving_UE_BS_c) = I_temp(serving_UE_BS_c) - Power(consider_CH).*path_loss_temp(serving_UE_BS_c);
v_opt = sum(sum(log2(1+S_temp./I_temp),2)./average_rate);
v_temp_pre = v_opt;
v_ini = v_opt;
S_opt = S_temp;
I_opt = I_temp;
Power_previous = Power;
trace_performance = zeros(1,iteration_max);

Power_ini=Power;

for ii = 1:iteration_max
    m_order = min(m_order_p*power_max,power_max);
    first_order = zeros(1,length(consider_CH));
    power_plus = power_max-sum(Power_previous);
    boundary_length = power_plus/sqrt(length(consider_CH));
    m_order = max(m_order,boundary_length);
    % calculate first-order deviation
    for j=consider_CH
        selected_UE = serving_UE_index(j);
        if ~isempty(selected_UE)
            H = path_loss(selected_UE,1,j);
            R = 1/average_rate(selected_UE);
            first_order(j) = H/(I(selected_UE,j)+Power_previous(j)*H)/R - ...
                sum(S(:,j)./average_rate./(I(:,j)+Power_previous(j).*path_loss(:,BS_node,j))./(S(:,j)+I(:,j)+Power_previous(j).*path_loss(:,BS_node,j)).*path_loss(:,BS_node,j));
        end
    end
    % end calculation
    infeasible_flag = 1;
    if m_order<=boundary_length
        Power(consider_CH) = Power_previous(consider_CH) + power_plus/length(consider_CH);
        infeasible_flag = 0;
    end
    
    while infeasible_flag == 1
        num_CH = length(consider_CH);
        boundary_length = power_plus/sqrt(num_CH);
        if m_order<=boundary_length || num_CH == 1
            Power(consider_CH) = Power_previous(consider_CH) + power_plus/num_CH ;
            infeasible_flag = -3;
            break
        end
        if abs((max(first_order(consider_CH))-min(first_order(consider_CH)))/max(first_order(consider_CH))) < 1e-3
            Power(consider_CH) = Power_previous(consider_CH) + power_plus/num_CH;
            infeasible_flag = -2;
            break
        end
        N = ones(num_CH,1); % D = num_CH
        P = eye(num_CH) - N*(N'*N)^(-1)*N'; % Q = P
        delta_P = P*first_order(consider_CH)'; % brackets of Eq. (40)
        if first_order(consider_CH)*delta_P < 0
            delta_P = delta_P*(-1);
        end
        mu_2 = sqrt(num_CH/(num_CH*m_order^2-power_plus^2)*(delta_P'*delta_P)); % mu_2 = v, m_order = beta_1 
        if ~isreal(mu_2)
            mu_2
            Power(consider_CH) = Power_previous(consider_CH) + power_plus/num_CH;
            infeasible_flag = -1;
            break
        end
        delta_P = delta_P/mu_2 + power_plus/num_CH; % Eq. (40)
        Power(consider_CH) = Power_previous(consider_CH) + delta_P';
        [~,min_CH] = min(first_order(Power(consider_CH)<Power_lower(consider_CH)));
        temp_CH = find(Power(consider_CH)<Power_lower(consider_CH));
        min_CH = temp_CH(min_CH);
        if ~isempty(min_CH)
            power_plus = power_plus + Power_previous(consider_CH(min_CH)) - Power_lower(consider_CH(min_CH));
if power_plus < -1e-8

end
if ~isreal(sqrt(m_order^2-sum((Power_previous(consider_CH(min_CH)) - Power_lower(consider_CH(min_CH)))).^2))

end
            m_order = sqrt(m_order^2-sum((Power_previous(consider_CH(min_CH)) - Power_lower(consider_CH(min_CH)))).^2);
            % get rid of CH with too small power, and give them Power_lower
            Power(consider_CH(min_CH)) = Power_lower(consider_CH(min_CH));
            consider_CH(min_CH) = [];
if power_plus + sum(Power) - sum(Power(consider_CH)) > power_max + 1e-6
    
end
        else
            infeasible_flag = 0;
        end
    end
if sum(Power)> power_max + 1e-6
    
end
    consider_CH = find(serving_UE_index(1,:)~=0);
    S_temp(serving_UE_BS_c) = S(serving_UE_BS_c) + Power(consider_CH).*path_loss_temp(serving_UE_BS_c);
    I_temp(:,consider_CH) = I(:,consider_CH) + (ones(num_UE,1)*Power(consider_CH)).*reshape(path_loss_temp(:,1,consider_CH),num_UE,length(consider_CH));
%     I_temp(serving_UE_BS_c) = I(serving_UE_BS_c) - Power(consider_CH).*path_loss_temp(serving_UE_BS_c);
    I_temp(serving_UE_BS_c) = I_temp(serving_UE_BS_c) - Power(consider_CH).*path_loss_temp(serving_UE_BS_c);
    v_temp = sum(sum(log2(1+S_temp./I_temp),2)./average_rate);
    trace_performance(ii) = v_temp;
    if v_temp > v_opt
        v_opt = v_temp;
        S_opt = S_temp;
        I_opt = I_temp;
        Power_opt = Power;
        Power_previous = Power;
    else
        para = 0;
        Power_previous = (1-para)*Power+para*Power_opt;
    end
    if ii~=1 && abs(v_temp_pre-v_temp)/v_temp < 0.0001
        break
    else
        v_temp_pre=v_temp;
    end
    if infeasible_flag < 0
        break
    end
    m_order_p = m_order_p*0.8;
end
if v_opt == v_ini
    BS_flag = 3;
else
    BS_flag = 0;
end
end