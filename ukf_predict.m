function [x_new, P_new] = ukf_predict(x_old, P_old)
 % ukf_predict generates 2n+1 sigma points, it propagates them using
 % state_update and outputs the mean and the covariance of the new sigma 
 % points distribution. 
 %
 % Input		
 % x_old		old state vector nx1
 % P_old		old state covariance matrix nxn
 %
 % Output
 % x_new		new state vector nx1
 % P_new		new state covariance matrix nxn

 %generation of sigma points and weights
[sp, wm, wc] = sigmapoint_gen(x_old, P_old);

%propagation using state update function
[x_new, P] = sigmapoint_propagate(sp, wm, wc);

% forced simmetry
P_new = (P + P')/2;

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [sp, wm, wc] = sigmapoint_gen(x_old, P_old)
% function that generates sigma points and weights

	n		= size(x_old,1);	% state dimension in which sigma points are generated
	nsp		= 2*n+1;			% number of sigma points
	P       = P_old;			% covariance matrix
	
	% parameters for simmetric distribution
	alpha	= 1;
	beta	= 2;
	k		= 0;
    
	lambda = alpha^2 * (n + k) -n;
	
	% generation of first order weights (weights mean)
	wm = 1/(2*(n + lambda)) * ones(nsp, 1);
    wm(1) = lambda/(n + lambda);
	
	% generation of second order weights (weights covariance)
	wc = 1/(2*(n + lambda)) * ones(nsp, 1);
	wc(1) = lambda/(n + lambda) + (1 - alpha^2 + beta);
	
	% generation sigma points sp
	dist	= zeros(n,n);
	sp		= zeros(n, nsp);
	
	% calculation of dist = sqrt((n+lambda)*P)
	[U, SIGMA, ~] = svd((n+lambda)* P);
	dist = U * sqrt(SIGMA);

	for i=1:nsp
		if i == 1
			sp(:,i) = x_old;
		elseif i <= 1+n
			sp(:,i) = x_old + dist(:,i-1);
		elseif i > 1+n && i<=nsp
			sp(:,i) = x_old - dist(:,i-1-n);
		end
	end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [x_new, P_new] = sigmapoint_propagate(sp, wm, wc)
% this function, using state_update, propagates the sp distribution
	
	n = size(sp(:,1),1);	% state dimension
	nsp = n*2 + 1;			% number of sigma points
	
	% init of the new sp distribution
	sp_new = zeros(n, nsp);	
	
	% propagation
	for j = 1:nsp
		 sp_new(:,j) = state_update(sp(:,j));
	end
	
	%%ANGLE in state
	angle_state = 0; % turn 1 for a state with only angles
	if angle_state == 1
		% mean 
		% note: to correctly compute mean and covariance of angular entity is 
		% better to obtain first the mean of cosin and sin and
		% then obtain the mean angle using atan2

		mean_c = zeros(n, 1);	% cosin mean
		mean_s = zeros(n, 1);	% sin mean
		for j = 1:nsp
			c = cos(sp_new(:,j));
			s = sin(sp_new(:,j));

			mean_c = mean_c + wm(j)*c;
			mean_s = mean_s + wm(j)*s;
		end

		x_new = atan2(mean_s, mean_c); % atan2 of mean

		% cov
		cc = zeros(n,n);		% cosin covariance
		cs = zeros(n,n);		% sin covariance
		for j = 1:nsp
			c = cos(sp_new(:,j) - x_new);
			s = sin(sp_new(:,j) - x_new);

			cc = cc + wc(j)*(c*c');
			cs = cs + wc(j)*(s*s');
		end
		P_new = atan2(cs,cc);	% atan2 of covariance
	
	else
		% no angle in state
		%mean
		mean = zeros(n, 1);
		for j = 1:nsp
			mean = mean + wm(j)*sp_new(:,j);
		end
		x_new = mean;
		% covar
		covar = zeros(n,n);
		for j = 1:nsp
			covar = covar + wc(j)*((sp_new(:,j) - x_new)*(sp_new(:,j) - x_new)');
		end
		P_new = covar;

	end

	
end

