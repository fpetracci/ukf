%% intro
clear, clc;

% time interval definition
t_ini	= 0;
t_fin	= 100;
n_step	= 1000;
d_step	= (t_fin - t_ini)/n_step;
time	= t_ini+d_step:d_step:t_fin;

%% signal generation
% initial values
n_state = 4;	% number of state
n_out = 3;		% number of output, to be consistent with the one inside state2meas fcn

x_0 = 10*ones(n_state, 1);			% initial state values
u = 0 ;								% input

% state noise
s_gaus	= randn(n_state, n_step);	% gaussian random signal having dev standard = 1
s_sigma	= 0.01;						% standard deviation for noise signal
s_noise	= s_sigma * s_gaus;

% measure noise
m_gaus	= randn(n_out, n_step);		% gaussian random signal having dev standard = 1
m_sigma	= 0.05;						% standard deviation for noise signal
m_noise	= m_sigma * m_gaus;

% state and output vectors pre-allocation
x_real = zeros(n_state, n_step);
y_meas = zeros(n_out, n_step);

% real signal generation and measure generation (simulated)
for i = 1 : n_step
	if i == 1
		x_old = x_0;
	end
	x_real(:, i) = state_update(x_old, u) + s_noise(:, i);			% xd = f(x) + s_noise
	y_meas(:, i) = state2meas(x_real(:, i), u) + m_noise(:, i);	% y = h(x) + m_noise
	x_old = x_real(:, i);
end

%% filter initialization

y = y_meas;		% measures vector
n_out = size(y,1);	% measurements dimension
n = size(y,2);		% number of time steps
n_state = 4;		% state dimension

e		   = zeros(n_out, n);
xCorrected = zeros(n_state, n);
PCorrected = zeros(n_state, n_state, n);
xPredicted = zeros(n_state, n);
PPredicted = zeros(n_state, n_state, n);

x_init = x_0;
P_init = eye(n_state).*s_sigma^2;	% noise state covariance matrix
R = eye(n_out).*m_sigma^2;			% noise measure covariance matrix
y_virt = zeros(n_out, n);

%% iterazione

for t = 1:n
	if t == 1
		% virtual measurements
		[y_virt(:, t), S_first, C] = ukf_virtmeas(x_init, P_init, u);
		S = R + S_first;
		% innovation / estimate's error
		e(:, t) = y(:, t) - y_virt(:, t);
		% correction step
		[xCorrected(:, t) , PCorrected(:, :, t)] = ukf_correct(...
								x_init, P_init, ...
								e(:, t), S, C);
	else
		% correction step
		[y_virt(:, t), S_first, C] = ukf_virtmeas(xPredicted(:, t-1), ...
											PPredicted(:, :, t-1), u);
		S = R + S_first;	% additive noise on measures
		e(:, t) = y(:, t) - y_virt(:, t);
		[xCorrected(:, t), PCorrected(:, :, t)] = ukf_correct(...
								xPredicted(:, t-1), PPredicted(:, :, t-1),...
								e(:, t), S, C);
	end
		% prediction step
		[xPredicted(:, t), PPredicted(:, :, t)] = ukf_predict(...
							xCorrected(:, t), PCorrected(:, :, t), u);
								
end

%% plots

% plot: real and estimated states
figure(1)
clf
for i = 1: n_state
	subplot(n_state, 1, i)
	plot(time, x_real(i,:), 'DisplayName', 'real state')
	hold on
	plot(time, xCorrected(i,:), 'linewidth', 1.2, 'DisplayName', 'estimated state')
	legend('Location', 'best')
	grid on
	xlabel('Time [s]')
	ylabel('State')
	axis tight
	title(['Real and Estimate of state ' num2str(i)])
end

% plot: measured and virtual outputs
figure(2)
clf
for i = 1:n_out
	subplot(n_out, 1, i)
	plot(time, y_meas(i,:), 'DisplayName', 'Measured outputs')
	hold on
	plot(time, y_virt(i,:), 'linewidth', 1.2, 'DisplayName', 'Virtual outputs')
	legend('Location', 'best')
	grid on
	xlabel('Time [s]')
	ylabel('Output')
	axis tight
	title(['Measured and Virtual output ' num2str(i)])
end