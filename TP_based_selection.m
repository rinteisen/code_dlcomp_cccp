function x0 = TP_based_selection(x0,num_UE,num_CH,num_TP,path_loss,noise,average_rate,power_bound,B_max,BS_Node,power_max)

% Tranfer binary matrix into combination of index
Power = zeros(num_TP,num_CH);
serving_UE_index = zeros(num_TP,num_CH);
S = zeros(num_UE,num_CH);
I = noise'*ones(1,num_CH);
object_optimal = zeros(1,num_CH);
for j=1:num_CH
    Power(:,j) = x0(num_UE+1,:,j)';
    I(:,j) = I(:,j) + path_loss(:,:,j)*Power(:,j);
    for b=1:num_TP
        serving_UE = find(x0(1:num_UE,b,j)==1);
        if ~isempty(serving_UE)
            serving_UE_index(b,j) = serving_UE;
            S(serving_UE,j) = S(serving_UE,j) + Power(b,j)*path_loss(serving_UE,b,j);
        end
    end
    I(:,j) = I(:,j) - S(:,j);
    object_optimal(j) = sum(log2(1+S(:,j)./I(:,j))./average_rate);
end

serving_UE_index_optimal = serving_UE_index;
Power_optimal = Power;
candidate_set = Candidate_Decide(num_UE,num_TP,num_CH,path_loss,average_rate,B_max,power_max,noise);
RRH_Node = 1:num_TP;
RRH_Node(BS_Node) = [];
temperature = 3;
serving_UE_index_pre = serving_UE_index;
loop_flag = zeros(1,num_CH);
consider_CH = 1:num_CH;

% Derive initial value of limited power
limit_power = ones(1,num_CH);
for b = RRH_Node
    for j = 1:num_CH
        candidate_UE = candidate_set(b,1:find(candidate_set(b,:,j)==0,1,'first')-1,j);
        if ~isempty(candidate_UE)
            u0 = log2(1+power_bound(1,b,j)*path_loss(candidate_UE,b,j)./(2*noise(candidate_UE)'))./average_rate(candidate_UE);
            [~,max_UE] = max(u0);
            limit_power(j) = limit_power(j)*max(min(noise(candidate_UE(max_UE))/path_loss(candidate_UE(max_UE),b,j),power_bound(1,b,j)*2),0);
        end
    end
end
limit_power = min(limit_power.^(1/length(RRH_Node)),reshape(power_bound(1,BS_Node,:),1,num_CH));

ratio_ini = 0.1;
ratio = ratio_ini;
count_temp_size = num_TP*3;

for count_temp = 1:count_temp_size
    ratio = ratio + (1-ratio_ini)/count_temp_size;
    TP_priority = BS_Node;
    [~,random_priority] = sort(rand(1,length(TP_priority)));
    TP_priority = TP_priority(random_priority);
    for b = TP_priority
        % Release influence of TP b
    [S(:,consider_CH) I(:,consider_CH)] = Release_Signal(length(consider_CH),S(:,consider_CH),I(:,consider_CH),path_loss(:,b,consider_CH),Power(b,consider_CH),serving_UE_index(b,consider_CH));
    power_upper = reshape(power_bound(1,b,consider_CH),1,length(consider_CH));
    % Update lmited power
    if count_temp == 1
        limit_power_pre = limit_power;
    else
        limit_power = ones(1,num_CH);
        for j = consider_CH
            consider_TP_size = 0;
            for bb = 1:num_TP
                if serving_UE_index(bb,j) ~= 0 && bb ~= b
                    limit_power(j) = limit_power(j)*max(min(noise(serving_UE_index(bb,j))/path_loss(serving_UE_index(bb,j),b,j),power_bound(1,bb,j)*2),0);
                    consider_TP_size = consider_TP_size + 1;
                end
            end
            limit_power(consider_CH) = min(limit_power(consider_CH).^(1/consider_TP_size),reshape(power_bound(1,b,consider_CH),1,length(consider_CH)));
            limit_power = sqrt(limit_power.*limit_power_pre);
            limit_power_pre = limit_power;
        end
    end
    % Select macro BS serving UE
    [serving_UE_index(b,consider_CH) Power(b,consider_CH) S(:,consider_CH) I(:,consider_CH)] = UE_Selection_BS_v4(S(:,consider_CH),I(:,consider_CH),path_loss(:,:,consider_CH),...
        candidate_set(b,:,consider_CH),num_TP,length(consider_CH),power_upper,average_rate,min(limit_power(consider_CH),power_upper),Power(b,consider_CH),b,temperature,ratio^2);
    end

    temperature = temperature*0.5;
    TP_priority = RRH_Node;
    [~,random_priority] = sort(rand(1,length(TP_priority)));
    TP_priority = TP_priority(random_priority);
    % Select RRH serving UE (random priority)
    for b = TP_priority
        [S(:,consider_CH) I(:,consider_CH)] = Release_Signal(length(consider_CH),S(:,consider_CH),I(:,consider_CH),path_loss(:,b,consider_CH),Power(b,consider_CH),serving_UE_index(b,consider_CH));
        power_upper = reshape(power_bound(1,b,consider_CH),1,length(consider_CH));
        [serving_UE_index(b,consider_CH) Power(b,consider_CH) S(:,consider_CH) I(:,consider_CH)] = UE_Selection_RRH_v4(S(:,consider_CH),I(:,consider_CH),path_loss(:,b,consider_CH),...
            candidate_set(b,:,consider_CH),length(consider_CH),serving_UE_index(BS_Node,consider_CH),power_upper,average_rate);

    end
    % Transfer to desired format (binary matrix)
    x0 = Transfer_Binary(num_UE,num_TP,num_CH,serving_UE_index,Power);
    
    % Record the optimal solution ever meet
    for j=consider_CH
        object = HetNetfun_power(x0(:,:,j),num_UE,1,noise,path_loss(:,:,j),average_rate);
        if object > object_optimal(j)
            object_optimal(j) = object;
            serving_UE_index_optimal(:,j) = serving_UE_index(:,j);
            Power_optimal(:,j) = Power(:,j);
        end
        if serving_UE_index_pre(:,j) == serving_UE_index(:,j)
            loop_flag(j) = loop_flag(j) + 1;
        else
            loop_flag(j) = 0;
        end
        if loop_flag(j) == 2
            consider_CH(consider_CH==j) = [];
        end
    end
    
    if isempty(consider_CH)
        break
    end
end

x0 = Transfer_Binary(num_UE,num_TP,num_CH,serving_UE_index_optimal,Power_optimal);
end