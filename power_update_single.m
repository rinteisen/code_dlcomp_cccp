function power = power_update_single(x0,num_UE,num_TP,path_loss,noise,average_rate,power_max,B_max)
% x0 = maximum_serving_node(x0,B_max,path_loss);
v_temp1=HetNetfun_power(x0,num_UE,1,noise,path_loss,average_rate);
w_max = 0.9;
w_min = 0.4;
c_min = 0.4;
c_max = 0.8;
r_max = 0.1;
r_min = 0.05;
shift_max = 0.4;
shift_min = 0.1;
population = 100;
max_generation = 50;
thresold_end = 0.02;
p_end = 0.9*population;
optimal_pop = ceil(population*0.1);
low_pop = ceil(population*0.3);

power_low = (max(power_max)-sum(power_max))/(num_TP-sum(power_max==0));
if power_low == 0
    power_low = max(power_max);
end
power_low = min(power_low,power_max);
power_coefficient = sum(x0(1:num_UE,:));
power_coefficient(power_coefficient>0) = 1;
power = x0(num_UE+1,:,:);
% power(1) = 0;
optimal_global = HetNetfun_power(x0,num_UE,1,noise,path_loss,average_rate);
optimal_local = zeros(1,population);
x_pop = rand(population,num_TP);
x_pop = x_pop.*(ones(population,1)*power_coefficient);
x_pop(1:population-optimal_pop-low_pop,:) = x_pop(1:population-optimal_pop-low_pop,:).*(ones(population-optimal_pop-low_pop,1)*power_max);
x_pop(population-optimal_pop-low_pop+1:population-low_pop,:) = x_pop(population-optimal_pop-low_pop+1:population-low_pop,:).*(ones(optimal_pop,1)*power_max);
x_pop(population-low_pop+1:population,:) = x_pop(population-low_pop+1:population,:).*(ones(low_pop,1)*power_low);

v_pop = zeros(population,num_TP);
x_local = x_pop;
for i=1:population
    optimal_local(i) = HetNetfun_power(cat(1,x0(1:num_UE,:,:),x_pop(i,:,:)),num_UE,1,noise,path_loss,average_rate);
end

generation = 0;
while(1)
    w = w_max - (w_max-w_min)*generation/max_generation;
    c1 = c_max - (c_max-c_min)*generation/max_generation;
    c2 = (c_max-c_min)*generation/max_generation + c_min;
    r = r_max - (r_max-r_min)*generation/max_generation;
    shift = shift_max - (shift_max-shift_min)*generation/max_generation;
    random_term = r*[(rand(population-low_pop,num_TP)-0.5).*(ones(population-low_pop,1)*power_max);(rand(low_pop,num_TP)-0.5).*(ones(low_pop,1)*power_low)];
    v_pop = w*v_pop + c1*(rand(population,num_TP)-shift).*(x_local-x_pop) + c2*(rand(population,num_TP)-shift).*(ones(population,1)*power-x_pop) + random_term;
    x_temp = x_pop + v_pop;
    x_temp(x_temp<0) = 0;
    x_pop = min(x_temp,ones(population,1)*power_max);
    power_r = path_loss(:,:)*x_pop';
    signal_r = (path_loss(:,:).*x0(1:num_UE,:))*x_pop';
    x = signal_r./(noise'*ones(1,population)+power_r-signal_r);
    rate = log2(1+x);
    v_temp = sum(rate./(average_rate*ones(1,population)),1);
    consider_population = find(v_temp>=optimal_local);
    for i=consider_population
        optimal_local(i) = v_temp(i);
        x_local(i,:) = x_pop(i,:);
        if v_temp(i) > optimal_global
%             a='PSO'
%             v_temp(i)/optimal_global
%             [x_pop(i,:);power]
            optimal_global = v_temp(i);
            power = x_pop(i,:);
        end
    end
    generation = generation + 1;
    if optimal_global > 0
        ll_end = sum(optimal_local>=optimal_global*(1-thresold_end));
    else
        ll_end = sum(optimal_local>=optimal_global*(1+thresold_end));
    end
    if ll_end >= p_end*population
        break;
    end
    if generation >= max_generation
        break;
    end
end
% optimal_global
% optimal_local
% x_local(190:200,:,:)
% power
% generation