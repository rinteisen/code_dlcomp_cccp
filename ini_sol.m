function [v_ini v_ded x] = ini_sol(num_UE,num_CH,num_TP,path_loss,noise,B_max,average_rate,power_max,mode,ber,K_ini)
x = zeros(num_UE+1,num_TP,num_CH);
power_max = power_max/num_CH;
strength = zeros(num_UE,num_TP);
v_ini = 0;
v_ded = 0;
for i=1:num_CH
    for j=1:num_TP
        strength(:,j) = power_max(j)*path_loss(:,j,i);
    end
    [s_temp tp_org] =sort(strength,2,'descend');
    s_temp = s_temp(:,1:B_max);
    tp_org = tp_org(:,1:B_max);
    rate = discrete_rate(sum(s_temp,2)./noise',mode,ber);
    [v_org UE_org] = max(rate./average_rate);
    tp_org = tp_org(UE_org,:);
    v_ded = v_ded + v_org;
    
    UE_reu = zeros(1,num_TP);
    v_reu = zeros(num_TP,1);
    tp_reu = [1:num_TP]';
    for j=1:num_TP
        s_temp = min(power_max)*path_loss(:,j,i);
        i_temp = min(power_max)*sum(path_loss(:,:,i),2);
        i_temp = i_temp - s_temp;
        [v_temp tp_temp] = sort(discrete_rate(s_temp./(noise'+i_temp),mode,ber)./average_rate,'descend');
        if B_max > 1
            v_reu(j) = v_reu(j) + max(v_temp);
            UE_reu(j) = tp_temp(1);
        else
            count = 1;
            UE_temp = tp_temp(count);
            while sum(ismember(UE_reu,UE_temp)==1) > 0
                count = count + 1;
                UE_temp = tp_temp(count);
            end
            v_reu(j) = v_reu(j) + v_temp(count);
            UE_reu(j) = tp_temp(count);
        end
%         [v_temp UE_reu(j)] = sort(discrete_rate(s_temp./(noise'+i_temp),mode,ber)./average_rate,'descend');
    end
    
    if sum(v_org) > sum(v_reu)*K_ini
        max_UE = UE_org;
        max_tp = tp_org;
        max_v = v_org;
        power_m = power_max;
    else
        max_UE = UE_reu;
        max_tp = tp_reu;
        max_v = v_reu;
        power_m = min(power_max)*ones(1,num_TP);
    end
    for j=1:length(max_UE)
        x(max_UE(j),max_tp(j,:),i) = 1;
        x(num_UE+1,max_tp(j,:),i) = power_m(max_tp(j,:));
    end
    power_r = path_loss(:,:,i)*x(num_UE+1,:,i)';
    signal_r = (path_loss(:,:,i).*x(1:num_UE,:,i))*x(num_UE+1,:,i)';
    sinr = signal_r./(noise'+power_r-signal_r);
    rate = discrete_rate(sinr,mode,ber);
    v_ini = v_ini + rate./average_rate;
end
end