function y = discrete_rate(x,mode,ber)
table_snr = [-5.5 -3.5 -2.2 -1.0 1.3 3.4 5.2 7.0 10.5 11.5 14 16 17 26.8];
level_snr = length(table_snr);
table_snr = 10.^(table_snr/10);
table_bit = [0.25 0.4 0.5 2/3 1 4/3 1.6 2 8/3 3.2 4 4.5 4.8 6];
vector_length = length(x);
capacity_gap = -1.5/log(ber);
if mode == 1
    y = zeros(size(x));
    for i=1:vector_length
        for j=1:level_snr
            if x(i) < table_snr(j)
                break;
            else
                y(i) = table_bit(j);
            end
        end
    end
elseif mode == 2
    y = log2(1+x);
elseif mode == 3
    y = log2(1+x*capacity_gap);
end