function utility = HetNetfun_power(xi,num_UE,num_CH,noise,path_loss,average_rate)
rate = zeros(num_UE,1);
for j=1:num_CH
    power_r = path_loss(:,:,j)*xi(num_UE+1,:,j)';
    signal_r = (path_loss(:,:,j).*xi(1:num_UE,:,j))*xi(num_UE+1,:,j)';
    x = signal_r./(noise'+power_r-signal_r);
    rate = rate + log2(1+x);
%     rate = rate + discrete_rate(x,mode,ber);
end
utility = sum(rate./average_rate);
end