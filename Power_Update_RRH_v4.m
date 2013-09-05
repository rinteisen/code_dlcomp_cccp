function [Power S I] = Power_Update_RRH_v4(num_UE,power_max,S,I,path_loss,serving_UE_RRH,average_rate,Power_lower)
% review down
% water-filling
error_range = -1e-8;
Power = zeros(1,length(serving_UE_RRH));
consider_CH = find(serving_UE_RRH~=0);
infeaisble_flag = 1;
average_rate = average_rate';
Power_lower = max(0,Power_lower);
while infeaisble_flag
    serving_UE_RRH_c = serving_UE_RRH(consider_CH) + num_UE*(consider_CH-1);
    lambda_zeros = 1./(S(serving_UE_RRH_c)+I(serving_UE_RRH_c)+path_loss(serving_UE_RRH_c).*Power_lower(consider_CH)).*path_loss(serving_UE_RRH_c)./average_rate(serving_UE_RRH(consider_CH));
    inverse_r = sum(1./average_rate(serving_UE_RRH(consider_CH)));
    K = sum((S(serving_UE_RRH_c)+I(serving_UE_RRH_c))./path_loss(serving_UE_RRH_c));
    lambda = inverse_r/(power_max+K); % Eq.(44)
    Power(consider_CH) = 1./average_rate(serving_UE_RRH(consider_CH))/lambda...
            - (S(serving_UE_RRH_c)+I(serving_UE_RRH_c))./path_loss(serving_UE_RRH_c); % Eq.(43)
    [~,min_CH] = min(lambda_zeros(Power(consider_CH)<Power_lower(consider_CH)-error_range));
    if isempty(min_CH)
        infeaisble_flag = 0;
    else
        % get rid of CH with too small power, and give them Power_lower
        temp_CH = find(Power(consider_CH)<Power_lower(consider_CH)-error_range);
        min_CH = temp_CH(min_CH);
        Power(consider_CH(min_CH)) = Power_lower(consider_CH(min_CH));
        power_max = power_max - Power_lower(consider_CH(min_CH));
        consider_CH(min_CH) = [];
    end
end

consider_CH = find(serving_UE_RRH~=0);
serving_UE_RRH_c = serving_UE_RRH(consider_CH) + num_UE*(consider_CH-1);
S(serving_UE_RRH_c) = S(serving_UE_RRH_c) + Power(consider_CH).*path_loss(serving_UE_RRH_c);
I(:,consider_CH) = I(:,consider_CH) + (ones(num_UE,1)*Power(consider_CH)).*reshape(path_loss(:,1,consider_CH),num_UE,length(consider_CH));
I(serving_UE_RRH_c) = I(serving_UE_RRH_c) - Power(consider_CH).*path_loss(serving_UE_RRH_c);

end