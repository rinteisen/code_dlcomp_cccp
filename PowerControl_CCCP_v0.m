%% PowerControl_CCCP_v0: function description
function [iter,p_opt_mat] = PowerControl_CCCP_v0(x0,pwr_max,pathloss,noise,rate_ave,nUEs,nTPs,nCHs,iter_max)
	p = reshape(x0(end,:,:),nTPs*nCHs,1);
	p_opt = p;
	c_opt = 1e10;
	c_pri = 1e10;
	Gy = (1-x0(1:nUEs,:,:)).*pathloss; % (1-x)*G
	A = repmat(eye(nTPs),[1,nCHs]);
	lb = zeros(size(p));
	ub = repmat(pwr_max',[nCHs,1]);
	for i = 1:iter_max
		% calculate f_cav'
		f_cav_1st = zeros(nTPs*nCHs,1);
		for c = 1:nCHs
			i_tp = (c-1)*nTPs + (1:nTPs);
			NPGy_R = log(2)*rate_ave.*(noise' + Gy(:,:,c)*p(i_tp));
			f_cav_1st(i_tp) = sum(Gy(:,:,c)./repmat(NPGy_R,[1,nTPs]),1)';
		end
		% find P_jk with constraint
		options = optimset('Display','off','Algorithm','sqp');
		[p,c_tmp,exitflag,output,grad,hessian] = fmincon(@HetNetfun_CCCP_inv,p,A,pwr_max',[],[],lb,ub,[],options);
		% c_tmp
		% check stop criterion
		p = abs(p);
		x0(end,:,:) = reshape(vec2mat(p,nTPs)',1,nTPs,nCHs);
		% f_tmp = HetNetfun_power(x0,nUEs,nCHs,noise,pathloss,rate_ave);
		% if f_tmp > f_opt
		% 	f_opt = f_tmp;
		% 	pwr_opt = pwr;
		% end
		if c_tmp < c_opt
			c_opt = c_tmp;
			p_opt = p;
		end
		% fprintf('performance cur : %e and step size: %e\n',c_tmp,output.stepsize);
		if abs((c_tmp-c_pri)/c_pri) < 1e-5 & i > 1
			break;
		else
			c_pri = c_tmp;
		end
	end
	p_opt_mat = reshape(vec2mat(p_opt,nTPs)',1,nTPs,nCHs);
	iter = i;


	function [y dy] = HetNetfun_CCCP_inv(x)
		y = f_cav_1st'*x;
		NPG = zeros(nUEs,nCHs);
		for j = 1:nCHs
			NPG(:,j) = noise' + pathloss(:,:,j)*x((j-1)*nTPs+(1:nTPs));
			y = y - sum(log2(NPG(:,j))./rate_ave);
		end

		% gradient 
		if nargout>1
			dy = f_cav_1st;
			for j = 1:nCHs
				i_tp = (j-1)*nTPs + (1:nTPs);
				NPG_R = log(2)*rate_ave.*NPG(:,j);
				dy(i_tp) = dy(i_tp) - sum( pathloss(:,:,j)./repmat(NPG_R,[1,nTPs]),1 )';		
			end
		end	
	end


end