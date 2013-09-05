function tt = main_v4_cccp_sqp(c_str,m_str,b_str,iter_max_str)
%%
% Input
% c_str: data index
% m_str: number of total UE
% b_str: number of total CH
%-----------------------------------
% Parameters
% r_type: type of resource allocation method. Assign 0 for proposed method (refer to main_v4_para for more detailed)
% p_type: type of power allocation method. Assign 0 for proposed method (refer to main_v4_para for more detailed)
% B_max: maximum number of CoMP TPs
% mode: type of rate (only used for initial solution function and having no effect now)
% ber: used for rate considering capacity gap (only used for initial solution function and having no effect now)
% K_ini: bias for initial solution function
% Following variables are unused for proposed method and should can be deleted
% bound_mode
% 
%%
cc = eval(c_str);
cc = cc - 1;
basepath = '../data/data_0602/';
iter_max = eval(iter_max_str);

round = load(fullfile(basepath,'iteration_time.mat'));

r_type = 0;
p_type = 0;
sample = round.sample;
mode = 2;
B_max = 5;
K_ini = 3;
thermal_noise = 10^(-151.4/10);
num_CH_all = round.num_CH;
num_UE = eval(m_str);
num_CH = eval(b_str);
bound_mode = 7;
ber = round.ber;
num_TP = round.num_TP;
T = 4*num_UE;

fairness = zeros(sample,T);
a_capacity = zeros(sample,T);
objective = zeros(sample,T);
ratio_performance = zeros(1,T);
ratio_serving_UE = zeros(num_TP,T);
proportional_fair_index = zeros(1,T);
average_rate_UE = zeros(num_UE,T,sample);
serving_time = zeros(num_UE,T,num_CH,sample);
history_solution = zeros(num_UE+1,num_TP,num_CH,T);
iteration_max = num_CH*6;
time = zeros(1,T);
round_time = zeros(1,T);
for i=1:sample
    tic
    int2str(i+cc*sample)
    data_op = load(fullfile(basepath,['channel_matlab_' int2str(i+cc*sample) '.mat']));
    power_max = data_op.power_max;
    noise = (data_op.noise(1:num_UE)-thermal_noise)/num_CH*num_CH_all+thermal_noise;
    path_loss_all = data_op.pathloss_UEtoCell(1:num_UE,:,1:num_CH,:);
    average_rate = 0.01*ones(num_UE,1);
    real_average = zeros(num_UE,1);
    
    capacity_last_round = 0;
    round_max = num_CH*2;
    for t=1:T
        path_loss = path_loss_all(:,:,:,t);
            %% inital value
        [v_ini v_ded x0] = ini_sol(num_UE,num_CH,num_TP,path_loss,noise,B_max,average_rate,power_max,mode,ber,K_ini);
        x = x0;
        x_j = sum(x0(1:num_UE,:,:),2);
        x_j(x_j>0) = 1;
        round = 0;
        power_bound = zeros(1,num_TP,num_CH);
        for b=1:num_TP
            power_bound(1,b,:) = power_max(b)/num_CH;
        end
        v_opt = sum(v_ini);
             %% interative solve the problem
        qi = power_max*max(path_loss,[],3)';
        mu_t = sum(log2(1+qi./noise)./average_rate');
        while round < round_max
            round = round + 1;
            % solve resource allocation sub-problem
            % x_pri = x0(1:num_UE,:,1)
            if t<=0
                x0 = TP_based_selection(x0,num_UE,num_CH,num_TP,path_loss,noise,average_rate,power_bound,B_max,1,power_max);
            else
                x0 = ResourceAllocation_CCCP_v1(x0,num_UE,num_CH,num_TP,path_loss,noise,average_rate,power_bound,mu_t,iter_max);
            end
            % x_post = x0(1:num_UE,:,1)
            v_temp = HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
            if v_temp > v_opt
                v_opt = v_temp;
                x = x0;
            end
            % solve power allocation sub-problem
            % [~,x0(num_UE+1,:,:)] = two_layer_power_control(x0,power_max,path_loss,noise,average_rate,num_UE,num_TP,num_CH,iteration_max);
            [~,x0(num_UE+1,:,:)] = PowerControl_CCCP_v0(x0,power_max,path_loss,noise,average_rate,num_UE,num_TP,num_CH,iter_max);
            % reallocation power budget
            for b=1:num_TP
                power_bound(1,b,:) = x0(num_UE+1,b,:) + (power_max(b)-sum(x0(num_UE+1,b,:)))/num_CH;
            end
            power_bound(power_bound<0) = 0;
            nor_power = x0(num_UE+1,:,:);
            nor_power(nor_power<0) = 0;
            x0(num_UE+1,:,:) = nor_power;
            v_temp = HetNetfun_power(x0,num_UE,num_CH,noise,path_loss,average_rate);
            
            x_i = sum(x0(1:num_UE,:,:),2);
            x_i(x_i>0) = 1;
            if v_temp > v_opt
                v_opt = v_temp;
                x = x0;
            end
            % exit criterion
            if isempty(find(x_i ~= x_j, 1))
                break
            else
                x_j = x_i;
            end
            
        end
        fprintf('round = %d\n',round);
        round_time(t) = round;
        rate = zeros(num_UE,1);
        for j=1:num_CH
            zero_power = x(num_UE+1,:,j)==0;
            x(:,zero_power,j) = 0;
            x_i = sum(x(1:num_UE,:,j),2);
            x_i(x_i>0) = 1;
            serving_time(:,t,j,i)= x_i;
            serving_UE = sum(x_i);
            if serving_UE~=0
                ratio_serving_UE(serving_UE,t) = ratio_serving_UE(serving_UE,t) + 1;
            end
            serving_time(:,t,j,i)= x_i;
            power_r = path_loss(:,:,j)*x(num_UE+1,:,j)';
            signal_r = (path_loss(:,:,j).*x(1:num_UE,:,j))*x(num_UE+1,:,j)';
            sinr = signal_r./(noise'+power_r-signal_r);
            rate = rate + log2(1+sinr);
        end
        history_solution(:,:,:,t) = x;
        ratio_performance(t) = ratio_performance(t) + sum(rate./average_rate)/sum(v_ini)/sample;
        if sum(rate./average_rate)/sum(v_ini) < 1-0.01
            rate = v_ini.*average_rate;
        end
        objective(i,t) = sum(rate./average_rate);
        min_avg_rate = min(average_rate)*0.5;
        average_rate = max((t-1)/t*average_rate + rate/t,min_avg_rate);
        
            %% calculate parameters we want
        real_average = (t-1)/t*real_average + rate/t;
        average_rate_UE(:,t,i) = real_average;
        fairness(i,t) = sum(real_average)^2/num_UE/(real_average'*real_average);
        a_capacity(i,t) = capacity_last_round + sum(rate);
        capacity_last_round = a_capacity(i,t);
        proportional_fair_index(t) = proportional_fair_index(t) + sum(log(average_rate))/sample;
        time(t) = toc;
%         time_temp=time(t)
        [int2str(i) '-' int2str(t)]
    end
end
tt = 1;
ave_fairness = sum(fairness/sample,1);
ave_capacity = sum(a_capacity/sample,1);
resultpath = './fmincon_sqp';
save(fullfile(resultpath,['result_UE_' m_str '_CH_' b_str '_Iter_' iter_max_str '_round_' int2str(cc/sample+1)]),'ave_fairness','ave_capacity','ratio_performance',...
    'ratio_serving_UE','objective','average_rate','proportional_fair_index','history_solution','average_rate_UE','serving_time','time','round_time');
end


