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
	x_real(:, i) = state_update(x_old) + s_noise(:, i);			% xd = f(x) + s_noise
	y_meas(:, i) = state2meas(x_real(:, i)) + m_noise(:, i);	% y = h(x) + m_noise
	x_old = x_real(:, i);
end

%% init filtro

n = n_step;
e = zeros(n_out, n);

xCorrected = zeros(n_state, n);
PCorrected = zeros(n_state, n_state, n);
xPredicted = zeros(n_state, n);
PPredicted = zeros(n_state, n_state, n);

x_init = x_0;
P_init = eye(n_state).*s_sigma^2;	% noise state covariance matrix
R = eye(n_out).*m_sigma^2;			% noise measure covariance matrix
y = y_meas;							% measures vector
y_virt = zeros(n_out, n);
%% iterazione

for t = 1:n
	if t == 1
		% predizione misure virtuali
		[y_virt(:, t), S_first, C] = ukf_virtmeas(x_init, P_init);
		S = R + S_first;
		% errore di stima
		e(:, t) = y(:, t) - y_virt(:, t);
		% correzione stato
		[xCorrected(:, t) , PCorrected(:, :, t)] = ukf_correct(...
								x_init, P_init, ...
								e(:, t), S, C);
	else
		% correction step
		[y_virt(:, t), S_first, C] = ukf_virtmeas(xPredicted(:, t-1), ...
											PPredicted(:, :, t-1));
		S = R + S_first;	% additive noise on measures
		e(:, t) = y(:, t) - y_virt(:, t);
		[xCorrected(:, t), PCorrected(:, :, t)] = ukf_correct(...
								xPredicted(:, t-1), PPredicted(:, :, t-1),...
								e(:, t), S, C);
	end
		% prediction step
		[xPredicted(:, t), PPredicted(:, :, t)] = ukf_predict(...
							xCorrected(:, t), PCorrected(:, :, t));
								
end

%% plots

% real state range for plots
x_min = min(min(x_real));
x_max = max(max(x_real));

% estimated state range for plots
x_c_min = min(min(xCorrected));
x_c_max = max(max(xCorrected));

% state range for plot
x_m = min([x_min, x_c_min, 0]);
x_M = max([x_max, x_c_max, 0]);

% plot: real states and outputs
figure(1);
clf
subplot(1,2,1)
plot(time, x_real)
title('real state')
subplot(1,2,2)
plot(time, y_meas)
title('measured output')

% plot: real and estimated states
figure(2)
clf
plot(time, x_real, 'DisplayName', 'real state')
hold on
plot(time, xCorrected, 'linewidth', 1.2, 'DisplayName', 'estimated state')
legend
grid on
xlabel('Time [s]')
ylabel('State')
ylim([x_m, x_M])
title('Real and Estimate states')

% plot: measured and virtual outputs
figure(3)
clf
plot(time, y_meas, 'DisplayName', 'Measured outputs')
hold on
plot(time, y_virt, 'linewidth', 1.2, 'DisplayName', 'Virtual outputs')
legend
grid on
xlabel('Time [s]')
ylabel('Output')
axis tight
title('Measured and Virtual outputs')