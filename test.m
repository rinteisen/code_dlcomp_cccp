%% test: function description
function [x] = test()
	myfun = @(p,c,d)...
		sum(c./(c+d*p))-1;  % parameterized function
	% c = [1 2 3 4 5];                    % parameter
	% d = [1 2 3 4 5];
	fun = @(p) myfun(p,[1 2 3 4 5],[1 2 3 4 5]);    % function of x alone
	x = fzero(fun,-2);
end
