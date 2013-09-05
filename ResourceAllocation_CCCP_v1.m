%% ResourceAllocation_CCCP_v1: exp
% ignore B_max
function x_out = ResourceAllocation_CCCP_v1(x0,nUEs,nCHs,nTPs,pathloss,noise,rate_ave,pwr_bound,mu_max,iter_max)
    % qi = power_max*max(pathloss,[],3)';
    % mu_max = sum(log2(1+qi./noise)./average_rate');
    x_out = x0;
	nRounds = 5;
	M = 1e-3;
	for j = 1:nCHs
		xp_j = x0(:,:,j);
		v0 = HetNetfun_power(xp_j,nUEs,1,noise,pathloss(:,:,j),rate_ave);
		S = pathloss(:,:,j).*repmat(x0(end,:,j),[nUEs,1]); % G.*P
		% Phi_fixed = sum(log2(noise'+sum(S,2))./rate_ave);
		%% mu_1->inf, mu_3->1
	    x_j = xp_j(1:nUEs,:);
	    x_j(x_j==0) = M;
	    x_j(x_j==1) = 1-M;
		mu = struct('v1',mu_max/1000,'v3',max(1.0,0.5*nUEs));
		r = struct('v1',(mu_max/mu.v1)^(1/nRounds),'v3',(1/mu.v3)^(1/(nRounds-1	)));
	    count = 0;
	    x_j_pri = x_j;
		while mu.v1 < mu_max*(1-M)
			count = count + 1;
			% normalization due to smaller mu.v3
			exp_g = 1.8;
			ind_over = find(sum(x_j,1)>mu.v3);
			for k = ind_over
				x_j_exp = x_j(:,k).^exp_g;
				x_j(:,k) = x_j_exp/sum(x_j_exp);
			end
			% cccp
			[~,x_j] = subproblem(x_j,nUEs,nTPs,noise,rate_ave,S,mu,M,iter_max,count);
			x_j_cur = x_j;
			% if isempty(find(x_j_cur-x_j_pri))
			% 	fprintf('there is no scheduling change!\n');
			% else
			% 	x_j_pri
			% 	count
			% 	mu.v1
			% 	x_j_cur
			% end
			x_j_pri = x_j_cur;
			% x_j(x_j<=0) = M;
			% x_tmp
			% mu.v1
			mu.v1 = mu.v1*r.v1;
			mu.v3 = max(1.0,mu.v3*r.v3);
		end
		% xp_j
		xp_j(1:nUEs,:) = x_j > 0.5;
		% xp_j
		% [~,xp_j(end,:)] = PowerControl_CCCP_v0(xp_j,pwr_bound(1,:,j),pathloss(:,:,j),noise,rate_ave,nUEs,nTPs,1);
		% x(end,:) = power_update_single(x,nUEs,nTPs,pathloss(:,:,j),noise,rate_ave,pwr_bound(1,:,j),5);
		v = HetNetfun_power(xp_j,nUEs,1,noise,pathloss(:,:,j),rate_ave);
		if v > v0
			% fprintf('channel = %d: gain = %e\n',j,(v-v0)/v0);
			x_out(:,:,j) = xp_j;
		end
	end

end

%% subproblem: for each channel
function [iter,x_opt_mat] = subproblem(x0_mat,nUEs,nTPs,noise,rate_ave,S,mu,M,iter_max,count)
	% max_min = [log(M),log(1-M)];
	% mu.v1 = mu.v1/max(abs(log(exp(max_min).*(1-exp(max_min)))));

	x0_mat = x0_mat'; % mat 2 vec
	x0 = x0_mat(:); % mat 2 vec
	z = log(1-x0);
	A = repmat(eye(nTPs)/log(M),[1,nUEs]);
	b = mu.v3*ones(nTPs,1);
	lb = log(M)*ones(size(z));
	ub = log(1-M)*ones(size(z));

	c_opt = 1e10;
	c_pri = 1e10;
	z_opt = z;

	% cccp procession
	for i = 1:iter_max
		% calculate f_cav'
		% f_cav_1st = mu.v1*(2-1./(1-exp(z)));
		f_cav_1st = mu.v1*(1-2*z);

		% options = optimset('Display','off');
		options = optimset('Display','off','GradObj','on','Algorithm','sqp','MaxIter',iter_max);
		[z,c_tmp,exitflag,output,grad,hessian] = fmincon(@HetNetfun_CCCP_inv,z,A,b,[],[],lb,ub,[],options);
		% c_tmp
		iter = output.iterations;
		% vec2mat(1-exp(z),nTPs)
		if c_tmp < c_opt
			c_opt = c_tmp;
			z_opt = z;
		end
		if abs((c_tmp-c_pri)/c_pri) < M & i > 1
			break;
		else
			c_pri = c_tmp;
		end
	end
	iter = i;
	x_opt_mat = vec2mat(1-exp(z_opt),nTPs);

	%% HetNetfun_CCCP_inv: function description
	function [y dy] = HetNetfun_CCCP_inv(x)
		% min convex:phi_interference + concave:A_penalty
		S_x = S.*exp(vec2mat(x,nTPs));
		y = sum(log2(noise' + sum(S_x,2))./rate_ave) + f_cav_1st'*x;

		% gradient 
		if nargout>1
			den = (log(2)*rate_ave).*(noise'+sum(S_x,2));
			dvex_mat = S_x./repmat(den,[1,nTPs]);
			dvex_mat = dvex_mat';
			dy = f_cav_1st + dvex_mat(:);
		end

	end
	
end

