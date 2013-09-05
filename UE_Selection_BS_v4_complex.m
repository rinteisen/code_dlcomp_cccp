function [serving_UE_index Power S I] = UE_Selection_BS_v4_complex(S,I,path_loss,candidate_UE_all,num_TP,num_CH,power_bound,average_rate,power_limit,Power,BS_node,temperature,ratio)
% review done
serving_UE_index = zeros(1,num_CH);
max_step = 5;
num_UE = size(S,1);
for j = 1:num_CH
    candidate_UE = candidate_UE_all(1,1:find(candidate_UE_all(1,:,j)==0,1,'first')-1,j);
    if ~isempty(candidate_UE)
        if power_bound(j) > 1e-6
            consider_type = 3;
            u0 = zeros(num_UE,consider_type);
            power_u0 = [Power(j) power_limit(j) power_bound(j)];
            for ii=1:consider_type
                I_est = I(:,j) + power_u0(ii)*path_loss(:,BS_node,j);
                u0(:,ii) = (-1)*log2(1+S(:,j)./I_est)./average_rate;
                S_est = S(:,j) + power_u0(ii)*path_loss(:,BS_node,j);
                u0(:,ii) = u0(:,ii) + log2(1+S_est./I(:,j))./average_rate;
            end
            u0 = u0(candidate_UE,:);
            max_UE = zeros(1,3);
            % choose UE with a propability distribution
            for ii=1:consider_type
                [u0_sorted u0_index] = sort(u0(:,ii),'descend');
                max_UE(ii) = u0_index(length(candidate_UE));
                for iii = 1:length(candidate_UE)-1
                    prob = rand(1);
                    pp = exp((u0_sorted(iii+1)-u0_sorted(iii))/u0_sorted(iii)/temperature*10);
                    if prob > pp/(1+pp)
                        max_UE(ii) = u0_index(iii);
                        break
                    end
                end
            end
%             [~,max_UE(1)] = max(u0(:,1));
%             [~,max_UE(2)] = max(u0(:,2));
%             [~,max_UE(3)] = max(u0(:,3));
            power_upper_UE = [sqrt(power_bound(j)*Power(j)) power_limit(j) power_bound(j)];
            power_lower_UE = [0 0 (power_bound(j)+max(power_upper_UE(1:2)))/2];
            standard_UE = 1;
            bias = max(abs(power_upper_UE(standard_UE)-power_upper_UE)/power_bound(j),0.1);
            same_UE = find(max_UE==max_UE(standard_UE));
            same_UE(same_UE==standard_UE) = [];
            max_UE(same_UE) = [];
    %         max_UE_u0(same_UE) = [];
            power_upper_UE(same_UE) = [];
            bias(same_UE) = [];

            candidate_UE = candidate_UE(max_UE);
            object_value = zeros(1,length(candidate_UE));
            object_power = zeros(1,length(candidate_UE));
            count = 0;
            Hi = path_loss(:,BS_node,j);
            pop_size = 5*num_TP;
            if ~isempty(power_limit)
                pop_size_limit = ceil(pop_size*0.8);
            else
                pop_size_limit = 0;
            end

            for i = candidate_UE
                count = count + 1;
                if i ~= candidate_UE(standard_UE)
                    power_upper = power_upper_UE(count);
                    power_lower = power_lower_UE(count);
                    pop = [0:pop_size-1]*(power_upper-power_lower)/(pop_size-1)+power_lower;
                    delta_step = power_upper/(pop_size-1);
                    if count == 1
                        pop = [pop [0:pop_size_limit-1]*(power_limit(j)-power_lower)/(pop_size_limit-1)];
                        delta_step_limit = power_limit(j)/(pop_size_limit-1);
                    else
                        delta_step_limit = 0;
                    end
                    S_m = S(:,j)*ones(1,length(pop));
                    I_m = I(:,j)*ones(1,length(pop));
                    Hi_m = Hi*ones(1,length(pop));
                    H1 = path_loss(i,BS_node,j);

                    for kk=1:max_step
        %                 delta_step_limit = delta_step_limit*0.5;
                        delta_step = delta_step*0.5;
                        der_pop = H1./(I(i,j)+pop*H1)/average_rate(i) - sum(S_m./(average_rate*ones(1,length(pop))).*Hi_m./(I_m+S_m+Hi*pop)./(I_m+Hi*pop),1);
                        if count == 1
                            delta_step_limit = delta_step_limit*0.5;
                            step = [delta_step*ones(1,pop_size) delta_step_limit*ones(1,pop_size_limit)];
                        else
                            step = delta_step*ones(1,pop_size);
                        end
                        temp = der_pop>0;
                        temp = 2*temp - 1;
                        pop = max(min(pop + temp.*step,power_upper),0);
                    end
                    pop = [pop max(min(Power(j),power_upper),power_lower)];
                    other_UE = sum(log2(1+[S_m S_m(:,1)]./([I_m I_m(:,1)]+Hi*pop)./(average_rate*ones(1,length(pop)))),1) - log2(1+S(i,j)./(I(i,j)+pop*H1))./average_rate(i);
                    [object_value(count) best_power] = max(log2(1+(S(i,j)+pop*H1)/I(i,j))/average_rate(i)+other_UE);
                    object_power(count) = pop(best_power);
                else
                    S_est = S(:,j);
                    S_est(i) = S_est(i) + Power(j)*path_loss(i,BS_node,j);
                    I_est = I(:,j) + Power(j)*path_loss(:,BS_node,j);
                    I_est(i) = I_est(i) - Power(j)*path_loss(i,BS_node,j);
                    object_value(count) = sum(log2(1+S_est./I_est)./average_rate);
                    object_power(count) = Power(j);
                end

            end
            max_u0 = max(object_value);
%             prob_selecting = exp(((object_value-max_u0)/max_u0+bias)/temperature);
            prob_selecting = exp((ratio*(object_value-max_u0)/max_u0+(1-ratio)*bias)*temperature);
            prob_selecting = prob_selecting/sum(prob_selecting);
            prob = rand(1);
            cdf = 0;

            for select_UE = 1:length(candidate_UE)
                cdf = cdf + prob_selecting(select_UE);
                if prob<=cdf
                    break
                end
            end
            if object_power(select_UE) == 0
                Power(j) = 0;
                [~,serving_UE_index(j)] = max(path_loss(candidate_UE,BS_node,j)./average_rate(candidate_UE)./(I(candidate_UE,j)+S(candidate_UE,j)));
                serving_UE_index(j) = candidate_UE(serving_UE_index(j));
            else
                serving_UE_index(j) = candidate_UE(select_UE);
                Power(j) = object_power(select_UE);
                S(serving_UE_index(j),j) = S(serving_UE_index(j),j) + Power(j)*path_loss(serving_UE_index(j),BS_node,j);
                path_loss(serving_UE_index(j),BS_node,j) = 0;
                I(:,j) = I(:,j) + Power(j)*path_loss(:,BS_node,j);
            end
        else
            Power(j) = 0;
            [~,serving_UE_index(j)] = max(path_loss(candidate_UE,BS_node,j)./average_rate(candidate_UE)./(I(candidate_UE,j)+S(candidate_UE,j)));
            serving_UE_index(j) = candidate_UE(serving_UE_index(j));
        end
if ~isreal(Power(j)) || Power(j)<-1e-8
    a='UE_Selection_BS_v4'
    Power
    j
end
    end
end
% Power = ones(size(Power))*4;
end