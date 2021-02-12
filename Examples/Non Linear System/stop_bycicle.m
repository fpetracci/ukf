% Questo script plotta alcuni plot interessanti

time	= out.state_ts.Time;
x_M		= out.state_ts.Data(:,1);
y_M		= out.state_ts.Data(:,2);
phi		= out.state_ts.Data(:,3);
theta_p	= out.state_ts.Data(:,4);
v_P		= out.state_ts.Data(:,5);
v_P_dot	= out.state_ts.Data(:,6);

x_rif   = out.rif.Data(:,1);
y_rif   = out.rif.Data(:,2);

x_e		= out.err.Data(:,1);
y_e		= out.err.Data(:,2);

s		= out.curv_coord.Data(:,1);
ds		= out.curv_coord.Data(:,2);
dds		= out.curv_coord.Data(:,3);

%% Uscite
figure(1)
clf
plot(time, x_M, 'r', 'DisplayName', 'x_p')
hold on
plot(time, y_M, 'b', 'DisplayName', 'y_p')
if size(x_rif,1) == 1
	plot(time, ones(length(time),1)*x_rif, 'm', 'DisplayName', 'x_{rif}')
	plot(time, ones(length(time),1)*y_rif, 'c', 'DisplayName', 'y_{rif}')
else
	plot(time, x_rif, 'm', 'DisplayName', 'x_{rif}')
	plot(time, y_rif, 'c', 'DisplayName', 'y_{rif}')
end
title('Uscite')
axis tight
grid on
legend
xlabel('Tempo [s]')
ylabel('Uscite [m]')

%% Stato: velocità
figure(2)
clf
plot(time, v_P, 'b', 'DisplayName', 'velocità')
hold on
plot(time, ds, 'r', 'DisplayName', 'vel lungo la traiettoria')
title('Velocità')
xlim([0  inf])
ylim([-inf inf])

grid on
legend
xlabel('Tempo [s]')
ylabel('velocità [m/s]')

%% Stato: Theta
figure(3)
clf
plot(time, theta_p*180/pi, 'r', 'DisplayName', '\theta')
title('\theta')
axis tight
grid on
legend
xlabel('Tempo [s]')
ylabel('\theta [grad]')

%% Stato: Phi
figure(4)
clf
plot(time, phi*180/pi, 'r', 'DisplayName', '\phi angolo sterzo')
title('\phi')
axis tight
grid on
legend
xlabel('Tempo [s]')
ylabel('\phi [grad]')

%% errore di inseguimento

figure(5)
clf
plot(time, x_e, 'r', 'DisplayName', 'error_x')
hold on
plot(time, y_e, 'b', 'DisplayName', 'error_y')
title('Errore di tracking')
axis tight
grid on
legend
xlabel('Tempo [s]')
ylabel('Errore [m]')

%% Mappa
figure(6)
clf

clf
if size(x_rif,1) == 1
	plot(x_rif, y_rif, 'rd', 'Linewidth', 1.0, 'DisplayName', 'Riferimento' )
else
	plot(x_rif, y_rif, 'r--', 'Linewidth', 1.0, 'DisplayName', 'Riferimento' )
end

hold on

if tipo_perc == 3 %clotoidi
	plot(xp, yp, 'co', 'DisplayName', 'Punti Inseriti')
elseif tipo_perc == 1 %punto
	plot(x_M_finale, y_M_finale, 'co', 'DisplayName', 'Punto finale')
end

plot(x_M_iniziale, y_M_iniziale, 'b*', 'DisplayName', 'Punto Iniziale' )

plot(x_M, y_M, 'b', 'DisplayName', 'Percorso')

title('Mappa')
axis equal
gap = 5;
xlim([min(min(x_M), min(x_rif)) - gap , max(max(x_M), max(x_rif)) + gap ])
ylim([min(min(y_M), min(y_rif)) - gap , max(max(y_M), max(y_rif)) + gap ])
grid on
legend
xlabel('Asse X [m]')
ylabel('Asse Y [m]')

%% mappa2
if 0
	figure(6)
	clf

	for i=1:10:length(time)
		if i~=1
			delete(freccia);
			delete(p1);
		else
			plot(x_M_iniziale, y_M_iniziale, 'b*', 'DisplayName', 'Punto Iniziale' )
		end
		hold on
		d = abs(v_P(i));
		x_v = d*cos(theta_p(i));
		y_v = d*sin(theta_p(i));
		freccia=quiver(x_M(i), y_M(i), x_v, y_v, 'r');
		p1 = plot(x_M(1:i), y_M(1:i), 'b', 'DisplayName', 'Percorso');
		xlim([min(x_M)-max(v_P) max(x_M)+max(v_P)])
		ylim([min(y_M)-max(v_P) max(y_M)+max(v_P)])
		grid on
		title(['Animazione, tempo ' num2str(time(i)) ' di ' num2str(time(end))])
		legend
		xlabel('Asse X [m]')
		ylabel('Asse Y [m]')

		drawnow
	end
end

