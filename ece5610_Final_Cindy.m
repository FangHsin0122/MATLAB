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

sys_C22 = tfest(data, 2, 2);
step_sys_C22 = c2d(sys_C22, Ts_sec, 'zoh')
[y_model, t_model] = step(step_sys_C22 * input_step_rt, t_x_sec(end));
plot(t_model, y_model, '-', 'DisplayName', 'Model - 2p2z');
hold on;

sys_C21 = tfest(data, 2, 1);
step_sys_C21 = c2d(sys_C21, Ts_sec, 'zoh')
[y_model, t_model] = step(step_sys_C21 * input_step_rt, t_x_sec(end));
plot(t_model, y_model, '-', 'DisplayName', 'Model - 2p1z');
hold on;

xlabel('Time (s)');
ylabel('Output Positions (pixels)');
title('Measured Data vs My Model (Negelect Delay)');
legend('Location','Best');
grid on;


figure
rlocus(step_sys_C22)
figure
rlocus(step_sys_C21)


figure
sys_D22 = tfest(data, 2, 2, 'Ts', Ts_sec);
[y_model, t_model] = step(sys_D22 * input_step_rt, t_x_sec(end));
plot(t_model, y_model, '-', 'DisplayName', 'Model - 2p2z(Discrete)');
hold on;

sys_D21 = tfest(data, 2, 1, 'Ts', Ts_sec);
[y_model, t_model] = step(sys_D21 * input_step_rt, t(end));
plot(t_model, y_model, '-', 'DisplayName', 'Model - 2p1z(Dis)');
hold on;

xlabel('Time (s)');
ylabel('Output Positions (pixels)');
title('Measured Data vs My Model (Negelect Delay)');
legend('Location','Best');
grid on;




%% calculate which model fit our date the best
% --- Simulate each model with the SAME input and time ---
y22 = lsim(sys_D22, u, t);
y21 = lsim(sys_D21, u, t);
%{
figure;
lsim(sys_S22, u, t)
figure;
lsim(sys_S21, u, t)
%}
% --- Compute errors ---
e22 = nonze_y - y22;
e21 = nonze_y - y21;

SSE22 = sum(e22.^2); % Sum of Squared Errors
SSE21 = sum(e21.^2);

MSE22 = mean(e22.^2); % Mean Squared Error
MSE21 = mean(e21.^2);

RMSE22 = sqrt(MSE22); % Root Mean Squared Error
RMSE21 = sqrt(MSE21);

fprintf('RMSE (2p2z): %.4f\n', RMSE22);
fprintf('RMSE (2p1z): %.4f\n', RMSE21);
%% My Model -> 2p2z with 2.0083 RMSE is close enough to the data
figure()

npoles = 2;
nzeros = 2;

sys_D = tfest(data, npoles, nzeros, 'Ts', Ts_sec);
[y_model, t_model] = step(sys_D * input_step_rt, t(end));
sys_D
plot(t_model, y_model, '-', 'DisplayName', 'My Model');
xlabel('Time (s)');
ylabel('Output Positions (pixels)');
legend('Location','Best');
grid on;

%% Q2
%{
(a) zero steady-state error to a step
(b) SPECS:
i. input to the head ≤ 200 counts;
ii. ≤ 10% overshoot; 
iii. insensitive to high-frequency noise;
iv. minimum settling time subject to all constraints

Once I have a controller that works without delay:
(c) Build the Simulink diagram with:
Controller block C(z), Plant block G(z), Step input = 9 pixels, 
Saturation block on u if we want to model 200-count limit explicitly, 
Verify SPECS.

(d) Add a pure time delay:

Continuous: 
1. multiply by exp(−sT_d)
2. Discrete: use a delay of N samples → block z^(−N) before G(z)
3. re-run the same tests, see how overshoot, settling time, and stability change
4. if necessary, retune K and z_0 (usually make bandwidth smaller)
so the specs are met even with delay, and then comment on the effect

%}

pole(sys_D)
%% (a) controllor C     
c_rt = 9;   % 9-pixel step input

z  = tf('z', Ts_sec);

[num, den] = tfdata(sys_D, 'v');
%num
%den
Gz = ( num(2)*z+num(3) ) / (den(1)*z^2 + den(2)*z + den(3));
z0 = 0.9; % initial guess for the zero
K  = 2; % can be tuned
Cz  = K * (z - z0) / (z - 1);

figure
rlocus(Gz)
figure
L = Cz*Gz;
rlocus(L)

%% (b) closed loop
Tcl = feedback(Cz*Gz, 1); % closed-loop r -> y

dcgain(Tcl) % should be finite and step error → 0

figure
step(c_rt*Tcl)
info = stepinfo(c_rt*Tcl) %ii. check overshoot


c_ut = c_rt* Cz / ( 1 + Cz*Gz)
