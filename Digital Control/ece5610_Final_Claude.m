% System Identification for Bisight Head/Image Processing System
clear all; close all; clc;

% Data
Ts_msec = 33.4; % sampling time in msec
Ts = Ts_msec / 1000; % convert to seconds
t_x_msec = [0, 33.4, 66.7, 100.1, 133.5, 166.9, 200.2, 233.6, 267.0, 300.4, ...
    333.7, 367.1, 400.5, 433.8, 467.2, 500.6, 534.0, 567.3, 600.7, 634.1, ...
    667.5, 700.8, 734.2, 767.6, 801.0];
t = t_x_msec / 1000; % convert to seconds
pos_y_pixels = [0, 0, 0, 5, 16, 8, 5, 12, 10, 7, 8, 11, 10, 8, 9, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9];
input_step = 110; % counts

% Plot measured data
figure(1);
plot(t, pos_y_pixels, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Time (sec)');
ylabel('Output Positions (pixels)');
title('Measured Step Response');
grid on;
hold on;

%% Step 1: Identify time delay and steady-state value
% Steady-state output (final value)
steady_state = mean(pos_y_pixels(end-3:end)); % average last few points
fprintf('Steady-state output: %.2f pixels\n', steady_state);

% Find approximate time delay (when output first becomes non-zero)
delay_idx = find(pos_y_pixels > 0, 1);
time_delay = t(delay_idx);
fprintf('Estimated time delay: %.3f sec (%.1f msec)\n', time_delay, time_delay*1000);

% DC gain: output/input
K_dc = steady_state / input_step;
fprintf('DC Gain: %.6f pixels/count\n', K_dc);

%% Step 2: Remove time delay from data for system identification
% Create time vector starting from when response begins
t_adjusted = t(delay_idx:end) - time_delay;
y_adjusted = pos_y_pixels(delay_idx:end);

%% Step 3: Fit a second-order system
% Try fitting to: y(t) = K*(1 - exp(-t/tau1) - (t/tau2)*exp(-t/tau2))
% This is an overdamped second-order response

% Initial guess for parameters
tau1_init = 0.05; % time constant 1
tau2_init = 0.10; % time constant 2

% Fitting using least squares
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
x0 = [K_dc, tau1_init, tau2_init];

% Objective function: sum of squared errors
objective = @(params) sum((y_adjusted - second_order_response(t_adjusted, params)).^2);

% Fit the model
x_fit = fminunc(objective, x0, options);

K_fit = x_fit(1);
tau1_fit = x_fit(2);
tau2_fit = x_fit(3);

fprintf('\n--- Fitted Second-Order System ---\n');
fprintf('DC Gain K: %.6f\n', K_fit);
fprintf('Time constant 1 (tau1): %.4f sec\n', tau1_fit);
fprintf('Time constant 2 (tau2): %.4f sec\n', tau2_fit);

% Generate fitted response
y_fit = second_order_response(t_adjusted, x_fit);

% Plot comparison
plot(t_adjusted + time_delay, y_fit, 's-', 'LineWidth', 2, 'MarkerSize', 4);
legend('Measured Data', 'Fitted Second-Order Model', 'Location', 'best');

%% Step 4: Derive transfer function
% For the form above, we can convert to poles
% Two poles: p1 = -1/tau1, p2 = -1/tau2
p1 = -1/tau1_fit;
p2 = -1/tau2_fit;

fprintf('\n--- Poles ---\n');
fprintf('p1 = %.4f\n', p1);
fprintf('p2 = %.4f\n', p2);

% Transfer function: G(s) = K*wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
% Or in factored form: G(s) = K / ((s/p1 + 1)(s/p2 + 1))
% Expanded: G(s) = K*p1*p2 / (s^2 + (p1+p2)*s + p1*p2)

wn_sq = p1 * p2;
wn = sqrt(wn_sq);
zeta = -(p1 + p2) / (2 * wn);

fprintf('\n--- Standard Second-Order Form ---\n');
fprintf('Natural frequency wn: %.4f rad/sec\n', wn);
fprintf('Damping ratio zeta: %.4f\n', zeta);

% Create transfer function
num = [K_fit * p1 * p2];
den = [1, -(p1+p2), p1*p2];
sys_continuous = tf(num, den);

fprintf('\nContinuous-time transfer function:\n');
disp(sys_continuous);

%% Step 5: Discretize for digital control implementation
sys_discrete = c2d(sys_continuous, Ts);
fprintf('\nDiscrete-time transfer function (Ts = %.4f sec):\n', Ts);
disp(sys_discrete);

% Plot step response of continuous-time model
figure(2);
t_sim = linspace(0, t_adjusted(end) + time_delay, 200);
y_sim = step(sys_continuous * time_delay, t_sim); % account for delay
plot(t_sim, y_sim, 'LineWidth', 2);
xlabel('Time (sec)');
ylabel('Output (pixels)');
title('Continuous-Time Second-Order System Step Response');
grid on;

%% Function Definition
function y = second_order_response(t, params)
    % Overdamped second-order response: y(t) = K*(1 - exp(-t/tau1) - (t/tau2)*exp(-t/tau2))
    K = params(1);
    tau1 = params(2);
    tau2 = params(3);
    
    y = K * (1 - exp(-t/tau1) - (t/tau2) .* exp(-t/tau2));
    y = max(y, 0); % avoid negative values
end



%% GPT
Tcl = feedback(C*G, 1);      % r -> y
Scl = feedback(1, C*G);      % r -> e
Gur = minreal(C * Scl);      % r -> u

t = 0:Ts_sec:2;              % simulate 2 seconds (adjust if needed)
r = 9 * ones(size(t));       % 9-pixel step

[y, ~] = lsim(Tcl, r, t);
[u, ~] = lsim(Gur, r, t);

max_u = max(abs(u))
%% 
L = Cz*Gz;
bode(L)    % inspect high-frequency gain
