function [serving_UE_index S I] = UE_Selection_RRH_ctd(S,I,path_loss,candidate_UE_all,num_CH,num_TP,BS_node,BS_UE,power_channel,average_rate,noise,ccc,bbb)

RRH_set = 1:num_TP;
RRH_set(RRH_set==BS_node) = [];
RRH_set_size = length(RRH_set);
serving_UE_index = zeros(RRH_set_size,num_CH);
least_ratio = 0;
for j = 1:num_CH
    least_power = max(power_channel(:,j))*1e-10;
    rate_matrix = zeros(RRH_set_size);
    UE_matrix = zeros(RRH_set_size);
    A = zeros(RRH_set_size,RRH_set_size*RRH_set_size);
    for b = RRH_set
        candidate_UE = candidate_UE_all(b,1:find(candidate_UE_all(b,:,j)==0,1,'first')-1,j);
        H1 = path_loss(:,b,j);
        candidate_UE = candidate_UE(candidate_UE~=BS_UE(j));
        consider_RRH = find(RRH_set==b);
        A_temp = zeros(RRH_set_size);
        A_temp(:,consider_RRH) = 1;
        A_temp(consider_RRH,:) = 1;
        A(consider_RRH,:) = reshape(A_temp,1,RRH_set_size*RRH_set_size);
        if ~isempty(candidate_UE)
            interference = sum(path_loss(candidate_UE,:,j).*(ones(length(candidate_UE),1)*power_channel(:,j)'),2) - ...
                path_loss(candidate_UE,b,j)*power_channel(b,j) + noise(candidate_UE)';
            [rate_matrix(consider_RRH,consider_RRH) UE_matrix(consider_RRH,consider_RRH)] = ...
                max(log2(1+power_channel(b,j)*H1(candidate_UE)./interference)./average_rate(candidate_UE));
            UE_matrix(consider_RRH,consider_RRH) = candidate_UE(UE_matrix(consider_RRH,consider_RRH));
            other_TP = RRH_set(RRH_set~=b);
            for bb = other_TP
                consider_RRH_2 = find(RRH_set==bb);
                [rate_matrix(consider_RRH,consider_RRH_2) UE_matrix(consider_RRH,consider_RRH_2)] = ...
                    max(log2(1+(power_channel(bb,j).*path_loss(candidate_UE,bb,j)+power_channel(b,j)*H1(candidate_UE))...
                    ./(interference-power_channel(bb,j).*path_loss(candidate_UE,bb,j)))./average_rate(candidate_UE));
                UE_matrix(consider_RRH,consider_RRH_2) = candidate_UE(UE_matrix(consider_RRH,consider_RRH_2));
                if rate_matrix(consider_RRH_2,consider_RRH) == -1
                    rate_matrix(consider_RRH,consider_RRH_2) = -1;
                elseif rate_matrix(consider_RRH,consider_RRH_2) > rate_matrix(consider_RRH_2,consider_RRH)
                    rate_matrix(consider_RRH_2,consider_RRH) = rate_matrix(consider_RRH,consider_RRH_2);
                    UE_matrix(consider_RRH_2,consider_RRH) = UE_matrix(consider_RRH,consider_RRH_2);
                else
                    rate_matrix(consider_RRH,consider_RRH_2) = rate_matrix(consider_RRH_2,consider_RRH);
                    UE_matrix(consider_RRH,consider_RRH_2) = UE_matrix(consider_RRH_2,consider_RRH);
                end
            end
        else
            UE_matrix(consider_RRH,consider_RRH) = BS_UE(j);
        end
    end
    rate_matrix_bs = [];
    if BS_UE(j)~=0
        RSRP_BS_UE = power_channel(:,j).*path_loss(BS_UE(j),:,j)';
        BS_signal = RSRP_BS_UE(BS_node);
        [~,sorted_index] = sort(RSRP_BS_UE,'descend');
        sorted_index(sorted_index==BS_node) = [];
        CoMP_RRH = sorted_index(1:2);
        rate_matrix_bs = -1*ones(1,4);
        if ~isempty(CoMP_RRH)
            interference = sum(path_loss(BS_UE(j),RRH_set,j).*power_channel(RRH_set,j)') + noise(BS_UE(j));
            
            consdier_RRH = [];
            consdier_TP = [BS_node consdier_RRH];
            interference_temp = interference - sum(path_loss(BS_UE(j),consdier_RRH,j).*power_channel(consdier_RRH,j)');
            rate_matrix_bs(1) = log2(1+sum(power_channel(consdier_TP,j).*path_loss(BS_UE(j),consdier_TP,j))/...
                interference_temp/average_rate(BS_UE(j)));
            
            if RSRP_BS_UE(CoMP_RRH(1))>least_ratio*BS_signal
                consdier_RRH = CoMP_RRH(1);
                consdier_TP = [BS_node consdier_RRH];
                interference_temp = interference - sum(path_loss(BS_UE(j),consdier_RRH,j).*power_channel(consdier_RRH,j)');
                rate_matrix_bs(3) = log2(1+sum(power_channel(consdier_TP,j).*path_loss(BS_UE(j),consdier_TP,j)')/...
                    interference_temp/average_rate(BS_UE(j)));
                
                if RSRP_BS_UE(CoMP_RRH(2))>least_ratio*BS_signal
                    consdier_RRH = CoMP_RRH(1:2)';
                    consdier_TP = [BS_node consdier_RRH];
                    interference_temp = interference - sum(path_loss(BS_UE(j),consdier_RRH,j).*power_channel(consdier_RRH,j)');
                    rate_matrix_bs(4) = log2(1+sum(power_channel(consdier_TP,j).*path_loss(BS_UE(j),consdier_TP,j)')/...
                        interference_temp/average_rate(BS_UE(j)));

                    consdier_RRH = CoMP_RRH(2);
                    consdier_TP = [BS_node consdier_RRH];
                    interference_temp = interference - sum(path_loss(BS_UE(j),consdier_RRH,j).*power_channel(consdier_RRH,j)');
                    rate_matrix_bs(2) = log2(1+sum(power_channel(consdier_TP,j).*path_loss(BS_UE(j),consdier_TP,j)')/...
                        interference_temp/average_rate(BS_UE(j)));
                end
            end
        end
    end
    f = [];
    A_map = zeros(RRH_set_size);
    range = 0:RRH_set_size;
    range = (2*RRH_set_size+1-range).*range/2;
    alpha = zeros(RRH_set_size*(RRH_set_size+1)/2,1);
    for b = 1:RRH_set_size
        if power_channel(RRH_set(b),j)*path_loss(UE_matrix(b,b),j) < noise(RRH_set(b))*least_ratio*10
            rate_matrix(b,b+1:RRH_set_size) = -1;
            rate_matrix(b+1:RRH_set_size,b) = -1;
        else
            ignore_RRH = find(abs(rate_matrix(b,b+1:RRH_set_size)-rate_matrix(b,b))/rate_matrix(b,b)<least_ratio*10);
            rate_matrix(b,b+ignore_RRH) = -1;
            rate_matrix(b+ignore_RRH,b) = -1;
        end
        alpha(range(b)+1)=1;
        f = [f rate_matrix(b,b:RRH_set_size)];
        A_map((b-1)*RRH_set_size+[b:RRH_set_size]) = 1;
    end

    A(:,A_map==0) = [];
    if BS_UE(j)~=0
        alpha = [alpha;1;zeros(3,1)];
        A = [A;zeros(1,length(f))];
        A = [A [zeros(RRH_set_size,4);ones(1,4)]];
        A(RRH_set==CoMP_RRH(1),length(f)+[3 4]) = 1;
        A(RRH_set==CoMP_RRH(2),length(f)+[2 4]) = 1;
        f = [f rate_matrix_bs];
    end
    lb = zeros(length(f),1);
    A_b = ones(size(A,1),1);
    alpha = bintprog(-1*f,A,A_b,[],[],alpha,optimset('Display','off'));

    alpha_RRH = zeros(RRH_set_size);
    for b = 1:RRH_set_size
        alpha_RRH(b,b:RRH_set_size) = alpha(range(b)+1:range(b+1))';
        alpha_RRH(b:RRH_set_size,b) = alpha(range(b)+1:range(b+1));
        aa = alpha_RRH(b,:);
        aaa = aa(aa>0);
        if isempty(aaa)
            alpha_RRH(b,:) = 0;
            alpha_RRH(b,b) = 1;
        else
            aa(aa<prod(aaa)^(1/length(aaa))) = 0;
            alpha_RRH(b,:) = aa;
        end
    end
    for b = 1:RRH_set_size
        if ~isempty(find(alpha_RRH(b,:)==1,1)) && power_channel(RRH_set(b),j) > least_power
            serving_UE_index(b,j) = UE_matrix(b,find(alpha_RRH(b,:)==1,1));
        elseif power_channel(RRH_set(b),j) > least_power
            serving_UE_index(b,j) = UE_matrix(b,b);
        else
            RRH_index = RRH_set(b);
            candidate_UE = candidate_UE_all(RRH_index,1:find(candidate_UE_all(RRH_index,:,j)==0,1,'first')-1,j);
            object_value_diff = path_loss(candidate_UE,RRH_index,j)./average_rate(candidate_UE)./I(candidate_UE,j);
            [~,select_UE] = max(object_value_diff);
            serving_UE_index(b,j) = candidate_UE(select_UE);
        end
    end
    if BS_UE(j)~=0 
        type_s = find(alpha(length(alpha)-3:length(alpha))==1,1);
        if isempty(type_s)
            type_s = 0;
        end
        switch type_s
            case 2
                serving_UE_index(RRH_set==CoMP_RRH(2),j) = BS_UE(j);
            case 3
                serving_UE_index(RRH_set==CoMP_RRH(1),j) = BS_UE(j);
            case 4
                serving_UE_index(RRH_set==CoMP_RRH(1),j) = BS_UE(j);
                serving_UE_index(RRH_set==CoMP_RRH(2),j) = BS_UE(j);
            otherwise
        end
    end
    if size(serving_UE_index,1)==5
    end
    for b=1:RRH_set_size
        RRH_UE = serving_UE_index(b,j);
        RRH_index = RRH_set(b);
        S(RRH_UE,j) = S(RRH_UE,j) + power_channel(RRH_index,j)*path_loss(RRH_UE,RRH_index,j);
        I(:,j) = I(:,j) + power_channel(RRH_index,j)*path_loss(:,RRH_index,j);
        I(RRH_UE,j) = I(RRH_UE,j) - power_channel(RRH_index,j)*path_loss(RRH_UE,RRH_index,j);
    end
end
end