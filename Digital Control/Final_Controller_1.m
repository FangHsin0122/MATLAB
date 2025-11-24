close all; clear all;

Ts_sec = 33.4/1000; %sec
t_x_sec = [0, 33.4, 66.7, 100.1, 133.5, 166.9, 200.2, 233.6, 267.0, 300.4, 333.7, 367.1, 400.5, 433.8, 467.2, 500.6, 534.0, 567.3, 600.7, 634.1, 667.5, 700.8, 734.2]/1000;  % secst

pos_y_pixels = [0, 5, 16, 8, 5, 12, 10, 7, 8, 11, 10, 8, 9, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9].';

input_step_rt = 110; %counts

figure(1)
plot(t_x_sec, pos_y_pixels, 'o-', 'DisplayName', 'Measured Step Response');
xlabel('Time (s)');
ylabel('Output Positions (pixels)');
title('Measured Step Resp (Origin)');
legend('Location','Best');
grid on;
%hold on;


u = input_step_rt * ones(length(pos_y_pixels), 1); % step input
data = iddata(pos_y_pixels, u, Ts_sec); % Identification dataset
t = (0 : length(pos_y_pixels)-1)* Ts_sec; % time vector in SECONDS for the non-delayed data
%t

figure(2)
% --- Compare with measured step response ---
plot(t, pos_y_pixels, 'ko', 'DisplayName', 'Measured Data(c2d)');
hold on;

sys_C22 = tfest(data, 2, 2); % continuous time controller
step_sys_D22 = c2d(sys_C22, Ts_sec, 'zoh') % discrete time controller
[y_model, t_model] = step(input_step_rt * step_sys_D22, t_x_sec(end));
%
plot(t_model, y_model, '-', 'DisplayName', 'Model - 2p2z');

xlabel('Time (s)');
ylabel('Output Positions (pixels)');
title('Measured Data vs My Model (Negelect Delay)');
legend('Location','Best');
grid on;


figure(3)
rlocus(step_sys_D22)

figure(4)
bode(step_sys_D22)

figure(5)
nyquist(step_sys_D22)

%%
%{
ζ = zeta: damping ratio
Tr: rise time = 1.8/ omega_n, for zeta = 0.5
Tp: peak time = pi/ (omega_n * sqrt(1 - zeta^2))
Ts: settling time = 4/ (zeta* omeaga_n)
Mp: overshoot = exp(- pi* zeta / ( sqrt(1 - zeta^2)))

r: magnitude = exp(- zeta* omega_n* T)
theta: phase = omega_n* T* ( sqrt(1 - zeta^2))
THEN, 
zeta = -ln(r)/(sqrt(theta^2 + (ln(r))^2)
omega_n = 1/T * (sqrt(theta^2 + (ln(r))^2)
tau: time constant = 1/(zeta* omega_n) = -T/ln(r)
%}

%% Get step response charateristic

StepINFO = stepinfo(y_model,t_model)

Mp = StepINFO.Overshoot;
Tsettling = StepINFO.SettlingTime;
Tr = StepINFO.RiseTime;

y_final = y_model(end)

%% P controller

Kp = 10;
D = Kp;

