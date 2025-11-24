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

xlabel('Time (s)');
ylabel('Output Positions (pixels)');
title('Measured Data vs My Model (Negelect Delay)');
legend('Location','Best');
grid on;


figure(3)
rlocus(step_sys_C22)

