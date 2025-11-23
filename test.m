%%  nov 18th
z




%%
%numerator = [0.1873, 0.1752];
%denominator = [1, -0.80562, 0.831];
ts = 1;
%sys = tf(numerator,denominator,ts)
%pole(sys)

%%
z = tf('z',ts);
G = 3*( 1 - exp(-2*ts)) / ( z - exp(-2*ts) );
F = z / (z-1);
%k = 10;
D = 1 + ( (0.1*z) / ( z - 1 ) );
sys = ( F * G) / (1 + D*G )
pole(sys)
%sys_min  = minreal(sys);            % exact
sys_tol  = minreal(sys,1e-6)       % near-cancellation with tolerance
pole(sys_tol)

%%
z = tf('z',ts);
Gp = 3*( 1 - exp(-2*ts)) / ( z - exp(-2*ts) );
G = 3*( 1 - exp(-2*ts)) / ( z - exp(-2*ts) );
F = z / (z-1);
%k = 10;
D = 1 + ( (0.1*z) / ( z - 1 ) );
sys = ( F * G) / (1 + D*G )
pole(sys)
%sys_min  = minreal(sys);            % exact
sys_tol  = minreal(sys,1e-6)       % near-cancellation with tolerance
pole(sys_tol)

%%
z = tf('z',ts);
Gp = 3*( 1 - exp(-2*ts)) / ( z - exp(-2*ts) );
G = 3*( 1 - exp(-2*ts)) / ( z - exp(-2*ts) );
ZOH = (z-1) / z;
%k = 10;
D = 1 + ( (0.1*z) / ( z - 1 ) );
sys = ( F * G) / (1 + D*G )
pole(sys)
%sys_min  = minreal(sys);            % exact
sys_tol  = minreal(sys,1e-6)       % near-cancellation with tolerance
pole(sys_tol)

%%
s = tf('s'); % create s domain variable
ts = 1;

Gs = 1/(s + 1); % continous time system
Gz = c2d(Gs, ts);

z = tf('z',ts);
D = (z-0.9)/(z-1);
sysC = D*Gz/(1+D*Gz)
pole(sysC)

sysC_tol  = minreal(sysC)       % near-cancellation with tolerance
pole(sysC_tol)
zero(sysC_tol)
[numC, denoC] = tfdata(sysC_tol, 'v');
omega_n = sqrt(denoC(3))
r = pole(sysC_tol)
log(r(1))
Ts = -4/ log(r(1))
%{
%If you prefer zeros/poles/gain:
[zC, pC, kC] = zpkdata(sysC_tol, 'v');
zC
pC
kC
%}

%% nov 14th - nyquist criteria
close all;

T = 1;
z = tf('z', T);
Cd = T/2 * ((z+1)/(z-1)^2);

C = (z-0.5)/(z-0.75);
figure(1)
rlocus(C*Cd);

figure(2)
%C = 1/(z-0.5);
%rlocus(C*Cd);

%% nov 14th - phase lead compensator example
num = 1;
den = poly([0 -1 -5]);

pzmap(num, den);
pause
rlocus(num, den);
pause
Knc = rlocfind(num, den)

% it's too slow, try lead compensation -> pole at z = 0.75, alphz = 0.1

nc = [1 0.75]; dc = [1 7.5];
num1 = conv(num, nc);
den1 = conv(den,dc);

rlocus(num1, den1)
pause

K1 = rlocfind(num1, den1)