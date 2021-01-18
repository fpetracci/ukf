%% intro
clear, clc;

animation_mode = 0; % put = 1 if you want an animation of the system in action

%% system definition - mass spring damper
% m xdd + b xd + k x = F
% y = x
m = 1;		% mass [kg]
k = 10;		% spring coeff [N/m]
b = 5;		% damper coeff [Ns/m]
F = 100;	% external force [N]


% state space matrices
A = [ 0 1 ; -k/m -b/m];
B = [0; 1];
C = [1 0];
D = 0;

sys_tc = ss(A,B,C,D);

%% discretization & time management
delta_t = 0.05; % Sampling period [s]
sys_td = c2d(sys_tc, delta_t);
global A_td B_td C_td D_td
A_td = sys_td.A;
B_td = sys_td.B;
C_td = sys_td.C;
D_td = sys_td.D;

% get dimensions
[n_out, n_state] = size(C_td);  % number of outputs, number of states
n_in = size(B_td,2);			% number of inputs

% time
t_end =  1.10 * max(damp(sys_tc));		% duration of simulation [s]
time	 = 0:delta_t:t_end;		% time vector
n_sample = length(time);		% number of samples

%% Simulation parameters

%--------------------------------------------------------------------------
% noise parameters
m_sigma	= 0.1;		% meas standard dev
s_sigma = 0.1;		% state standard dev

%--------------------------------------------------------------------------
% initial state
x0 = [1; zeros(n_state-1,1)];		

%--------------------------------------------------------------------------
%  inputs, uncomment the one to choose
u0 = F;
u = zeros(n_in, n_sample);

% constant input
for i = 1:n_sample
	u(:,i)  = F * ones(n_in, 1);		
end

% % sinusoidal input
% n_repeats = 3;
% for i = 1:n_sample
% 	u(:,i)  = F*sin(time(i)/t_end * 2*n_repeats* pi) * ones(n_in, 1);		
% end

% % constant + sinusoidal input
% n_repeats = 20;
% for i = 1:n_sample
% 	u(:,i)  = (F + F/2*sin(time(i)/t_end * 2*n_repeats* pi)) * ones(n_in, 1);		
% end

%%  preallocations

n = n_sample;	% number of time steps

%signals outside the filter	
x_real = zeros(n_state, n);	% real state (with noise)
y_meas = zeros(n_out, n);	% measurements (with noise)

% filter inner signals
x_init	   = x0;							% initial state
P_init = eye(n_state).*(s_sigma^2 + 0.1);	% initial noise state covariance matrix
y_virt	   = zeros(n_out, n);				% virtual measurements		
e		   = zeros(n_out, n);				% innovation / estimation error = y_meas - y_virt
R = eye(n_out).*m_sigma^2;					% noise measure covariance matrix
xCorrected = zeros(n_state, n);				% correction state
PCorrected = zeros(n_state, n_state, n);	% covar matrix of corrected state
xPredicted = zeros(n_state, n);				% prediction state
PPredicted = zeros(n_state, n_state, n);	% covar matrix of predicted state

%% iterations

for t = 1:n_sample
	%----------------------------------------------------------------------
	% mechanical system
	if t == 1
		x_old = x0;
	end
	x_real(:, t) = state_update(x_old, u(:,t))		+ s_sigma * randn(n_state,1);			
	y_meas(:, t) = state2meas(x_real(:, t), u(:,t))	+ m_sigma * randn(n_out,1);	
	x_old = x_real(:, t);
	
	%----------------------------------------------------------------------
	% filter
	
	if t == 1
		% virtual measurements
		[y_virt(:, t), S_first, C] = ukf_virtmeas(x_init, P_init, u(:,t));
		S = R + S_first;
		% innovation / estimate's error
		e(:, t) = y_meas(:, t) - y_virt(:, t);
		% correction step
		[xCorrected(:, t) , PCorrected(:, :, t)] = ukf_correct(...
								x_init, P_init, ...
								e(:, t), S, C);
	else
		% correction step
		[y_virt(:, t), S_first, C] = ukf_virtmeas(xPredicted(:, t-1), ...
											PPredicted(:, :, t-1), ...
											u(:,t));
		S = R + S_first;	% additive noise on measures
		e(:, t) = y_meas(:, t) - y_virt(:, t);
		[xCorrected(:, t), PCorrected(:, :, t)] = ukf_correct(...
								xPredicted(:, t-1), PPredicted(:, :, t-1),...
								e(:, t), S, C);
	end
		% prediction step
		[xPredicted(:, t), PPredicted(:, :, t)] = ukf_predict(...
							xCorrected(:, t), PCorrected(:, :, t), u(:,t));
	
	%----------------------------------------------------------------------
	% calculation of up and lower bounds of uncertainty
		if t == 1
			upstd = [];
			dnstd = [];
			yupstd = [];
			ydnstd = [];
		end
		upstd = cat(2, upstd, xCorrected(:,t) + 3 * diag(sqrt(PCorrected(:,:,t))));
		dnstd = cat(2, dnstd, xCorrected(:,t) - 3 * diag(sqrt(PCorrected(:,:,t))));
		yupstd = cat(2, yupstd, y_virt(:,t) + 3 * diag(sqrt(S)));
		ydnstd = cat(2, ydnstd, y_virt(:,t) - 3 * diag(sqrt(S)));	
	
	%----------------------------------------------------------------------
	% animation
	if animation_mode == 1
		
		figure(1)
		% clear figure
		if t == 1
			clf
		else
			delete(p1);
			delete(p2);
			delete(p3);
			delete(p4);
		end	
		
		for i = 1:n_state
			subplot(n_state+n_out, 1, i)
			p1 = plot(time(1:t), x_real(i,1:t), 'ko','DisplayName', 'real state');
			hold on
			p2 = plot(time(1:t), xCorrected(i,1:t), 'b','linewidth', 1.2, 'DisplayName', 'estimated state');
			p3 = plot(time(1:t), upstd(i,1:t), 'r', 'DisplayName', 'up');
			p4 = plot(time(1:t), dnstd(i,1:t), 'r', 'DisplayName', 'low');
% 			legend
			grid on
			xlabel('Time [s]')
			ylabel(['State ' num2str(i)])
			axis tight
			title(['Real and Estimate of state ' num2str(i)])
		end
		for i = n_state+1:n_state+n_out
			ii = i - n_state;
			subplot(n_state+n_out, 1, i)
			p1 = plot(time(1:t), y_meas(ii,1:t), 'ko','DisplayName', 'real meas');
			hold on
			p5 = plot(time(1:t), y_virt(ii,1:t), 'b','DisplayName', 'virtual meas');
			p6 = plot(time(1:t), yupstd(ii,1:t), 'r', 'DisplayName', 'up');
			p7 = plot(time(1:t), ydnstd(ii,1:t), 'r', 'DisplayName', 'low');
% 			legend
			grid on
			xlabel('Time [s]')
			ylabel(['Output ' num2str(ii)] )
			axis tight
			title(['Real and virtual output ' num2str(ii)])
		end
		
		drawnow
	end
	
end

%% final plots

% plot: real and estimated states
figure(2)
clf
for i = 1: n_state
	subplot(n_state, 1, i)
	plot(time, x_real(i,:), 'r', 'DisplayName', 'real state')
	hold on
	plot(time, xCorrected(i,:), 'b','linewidth', 1.2, 'DisplayName', 'estimated state')
	plot(time, upstd(i,:), ':b', 'DisplayName', 'up bound');
	plot(time, dnstd(i,:), ':b', 'DisplayName', 'low bound');
	legend('Location', 'best')
	grid on
	xlabel('Time [s]')
	ylabel('State')
	axis tight
	title(['Real and Estimate of state ' num2str(i)])
end

% plot: measured and virtual outputs
figure(3)
clf
for i = 1:n_out
	subplot(n_out, 1, i)
	plot(time, y_meas(i,:), 'r','DisplayName', 'Measured outputs')
	hold on
	plot(time, y_virt(i,:), 'b','linewidth', 1.2, 'DisplayName', 'Virtual outputs')
	plot(time, yupstd(i,:), ':b', 'DisplayName', 'up bound');
	plot(time, ydnstd(i,:), ':b', 'DisplayName', 'low bound');
	legend('Location', 'best')
	grid on
	xlabel('Time [s]')
	ylabel('Output')
	axis tight
	title(['Measured and Virtual output ' num2str(i)])
end
