function [y_virt, S_virt, C] = ukf_virtmeas(x, P)
 % ukf_virtmeas generates 2n+1 sigma points, it propagates them using
 % state2meas and outputs the mean and the covariance of the new sigma 
 % points distribution and the crosscorelation between state and virtual 
 % measurements. 
 %
 % Input		
 % x			state vector nx1
 % P			state covariance matrix nxn
 %
 % Output
 % y_virt		virtual measurements vector mx1
 % S_virt		virtual measurements covariance matrix mxm
 % C			state-measurements crosscorelation matrix nxm
 
% generation of sigmapoints and their weights
[sp, wm, wc] = sigmapoint_gen(x, P);

% propagation of state vector using state2meas
[~, y_virt, esse, C] = spstate2meas(sp, x, wm, wc);

% forced simmetry							
S_virt = (esse + esse')/2;

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [sp, wm, wc] = sigmapoint_gen(x, P)
% function that generates sigma points and weights

	n		= size(x,1);		% state dimension in which sigma points are generated
	nsp		= 2*n+1;			% number of sigma points
	P       = P;				% covariance matrix
	
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
	
	% computation of dist = sqrt((n+lambda)*P)
	[U, SIGMA, ~] = svd((n+lambda)* P);
	dist = U * sqrt(SIGMA);

	for i=1:nsp
		if i == 1
			sp(:,i) = x;
		elseif i <= 1+n
			sp(:,i) = x + dist(:,i-1);
		elseif i > 1+n && i<=nsp
			sp(:,i) = x - dist(:,i-1-n);
		end
	end


end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [sp_meas_new, y_virt, S_new, C] = spstate2meas(sp, x, wm, wc)
% this function propagates sigmapoint using state2meas and then computes
% mean, covariance matrix and crosscorelation matrix
	
	n  = size(x,1);		% state dimension
	nsp = 2*n + 1;		% number of sigma points
	
	% init sigma points propagated
	n_meas = size(state2meas(x),1);
	sp_meas_new = zeros(n_meas, nsp);
	
	% propagation through meas function
	for j = 1:nsp
		sp_meas_new(:,j) = state2meas(sp(:,j));
	end
	
	% virt meas mean
	y_virt = zeros(n_meas, 1);
	for j = 1:nsp
		y_virt = y_virt + wm(j)*sp_meas_new(:,j);
	end

	% virt meas matrix covariance
	S_new = zeros(n_meas, n_meas);  
	for j = 1:nsp
		 S_new = S_new + wc(j)*(sp_meas_new(:,j) - y_virt)*(sp_meas_new(:,j) - y_virt)';
	end
	
	% crosscov state - meas virt
	
	angle_state = 0; % turn 1 for a state with only angles
	
	C = zeros(n, n_meas);
	if angle_state == 1
		% angles in the state
		for j = 1:nsp
			% angular difference obtained from cosin and sin of difference
			% between angles
			ang = (sp(:,j)-x);				% angle diff
			a = atan2(sin(ang),cos(ang));	% atan2 of angle diff

			C = C + wc(j)*( a * (sp_meas_new(:,j)- y_virt(:,1))');
		end
	else
		% no angle in the state
		for j = 1:nsp
			diff = (sp(:,j)-x);
			C = C + wc(j)*( diff * (sp_meas_new(:,j)- y_virt(:,1))');
		end
	end
	
end