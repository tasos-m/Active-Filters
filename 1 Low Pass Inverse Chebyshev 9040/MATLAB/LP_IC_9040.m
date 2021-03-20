%% Anastasios Mouratidis 9040
%% Low Pass Inverse Chebyshev
clear;
close all;

a1 = 9;
a2 = 0;
a3 = 4;
a4 = 0;

m = 2; % a3 = 4
fp = 1.1*(3+m)*1000;
fs = 1.7*fp;
amin = 25 + (max(1,a3)-5)*(3/4);
amax = 0.5 + (max(1,a4)-5)/16;

ws = 2*pi*fs;
wp = 2*pi*fp;

Ws = 1;
Wp = wp/ws;

n = acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(1/Wp);
n = ceil(n);

e = 1/sqrt(10^(amin/10)-1);
a = (1/n)*(asinh(1/e));

% sixnotita imiseias isxios 
whp = 1/(cosh((1/n)*acosh(1/e)));
%% Gwnies Butterworth

y1 = 0;
y2 = 36;
y3 = -36;
y4 = 72;
y5 = -72;

%% Poloi
p1 = -sinh(a)*cosd(y1) + 1i*cosh(a)*sind(y1);
p2 = -sinh(a)*cosd(y2) + 1i*cosh(a)*sind(y2);
p3 = -sinh(a)*cosd(y3) + 1i*cosh(a)*sind(y3);
p4 = -sinh(a)*cosd(y4) + 1i*cosh(a)*sind(y4);
p5 = -sinh(a)*cosd(y5) + 1i*cosh(a)*sind(y5);

Wo_1  = sqrt((real(p1))^2+(imag(p1))^2);
Wo_23 = sqrt((real(p2))^2+(imag(p2))^2);
Wo_45 = sqrt((real(p4))^2+(imag(p4))^2);

Q1  = 1/(2*cos(atan(imag(p1)/real(p1))));
Q23 = 1/(2*cos(atan(imag(p2)/real(p2))));
Q45 = 1/(2*cos(atan(imag(p4)/real(p4))));

%% Antistrofi Polwn
w1  = 1/Wo_1;
w23 = 1/Wo_23;
w45 = 1/Wo_45;
%% Midenika
Z1 = sec(pi/10);
Z2 = sec(3*pi/10);
Z3 = sec(5*pi/10); % to Z3 einai sto apeiro ara den metasximatizetai

%% Piknotis
C = 10^(-7);
%% Proti monada
%  0.3776 < k11 < 1 kai epilegw k11 = 0.9
% W0 = 1;
Wz1 = Z2/w23;
k11 = 0.9;
R11 = 2/((k11*Wz1^2)-1);
R12 = 1/(1-k11);
R13 = (1/2)*((k11/(Q23^2))+(k11*Wz1^2)-1);
R14 = 1/k11;
R15 = 1;
R16 = 1;
C11 = k11/(2*Q23);
C12 = 2*Q23;
k1 = 1/(0.5*((k11/(Q23^2))+(k11*Wz1^2)+1));
%% Klimakopoihsh
kf1 = ws * w23;
km1 = C11 / (kf1*C);
C11new = C;
C12new = C12/(km1 * kf1);
R11new = R11 * km1;
R12new = R12 * km1;
R13new = R13 * km1;
R14new = R14 * km1;
R15new = R15 * km1;
R16new = R16 * km1;
%% Deuteri monada
% 0.6137 < k21 < 1 kai epilegw k21 = 0.9
Wz2 = Z1/w45;
k21 = 0.9;
R21 = 2/((k21*Wz2^2)-1);
R22 = 1/(1-k21);
R23 = (1/2)*((k21/(Q45^2))+(k21*Wz2^2)-1);
R24 = 1/k21;
R25 = 1;
R26 = 1;
C21 = k21/(2*Q45);
C22 = 2*Q45;
k2 = 1/(0.5*((k21/(Q45^2))+(k21*Wz2^2)+1));
%% Klimakopoihsh
kf2 = ws * w45;
km2 = C21 / (kf2*C);
C21new = C;
C22new = C22/(km2 * kf2);
R21new = R21 * km2;
R22new = R22 * km2;
R23new = R23 * km2;
R24new = R24 * km2;
R25new = R25 * km2;
R26new = R26 * km2;
%% Triti monada 
C31 = 1;
R31 = 1 / (C31 * w1);
%% Klimakopoihsh
kf3 = ws;
km3 = C31 / (kf3*C);
C31new = C;
R31new = R31 * km3;
%% Rythmisi kerdous sta 0dB
ktotal_low = k1*k2*(Wz1^2)*(Wz2^2);
a_k = 1/ktotal_low;
%% Synartiseis Metaforas

T1 = tf([k1 0 k1*(Z2*ws)^2], [1, (w23*ws)/Q23, (w23*ws)^2]);
T2 = tf([k2 0 k2*(Z1*ws)^2], [1, (w45*ws)/Q45, (w45*ws)^2]);
T3 = tf(w1*ws, [1, w1*ws]);

T_total = a_k * T1 * T2 *T3;

plot_transfer_function(T1, [fp fs]);
plot_transfer_function(T2, [fp fs]);
plot_transfer_function(T3, [fp fs]);
plot_transfer_function(T_total, [100 fp fs]);

Invsys=inv(T_total);
plot_transfer_function(Invsys,[100 fp fs])

ltiview({'bodemag'}, T1, T2, T3, T_total)

%% Fourier Analysis

fsig = 2*1000;

Tp = 10* 1/fsig;
f_s = 10^6; 
dt= 1/f_s;
t = 0:dt:(Tp - dt);

x = (square(2*pi*2000*(t),40)+1)/2; 

figure;
plot(t,x);
axis([0 inf -0.5 1.5]);
title('Input signal');


y = lsim(T_total,x,t);
figure;
plot(t,y);
axis([0 inf -0.5 1.5]);
title('Output signal');


Xft = fft(x);
Lx = length(x);

P_2 = abs(Xft/Lx);
P_1 = P_2(1:Lx/2+1);
P_1(2:end-1) = 2*P_1(2:end-1);

fu = f_s*(0:(Lx/2))/Lx;
figure;
plot(fu, P_1);
axis([0 80000 0 inf]);
title('Single-Sided Amplitude Spectrum of of Input signal');

Yft = fft(y);
Ly = length(y);

P2Y = abs(Yft/Ly);
P1Y = P2Y(1:Ly/2+1);
P1Y(2:end-1) = 2*P1Y(2:end-1);

fy = f_s*(0:(Ly/2))/Ly;
figure;
plot(fy, P1Y);
axis([0 10^4 0 inf]);
title('Single-Sided Amplitude Spectrum of Output signal');


