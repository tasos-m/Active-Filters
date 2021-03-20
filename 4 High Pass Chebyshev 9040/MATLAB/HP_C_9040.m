%% Anastasios Mouratidis 9040
%% High Pass Chebyshev
clear;
close all;

a1 = 9;
a2 = 0;
a3 = 4;
a4 = 0;

m = 0; % a3 = 0
fp = (3 + m)*1000;
fs = fp/1.8;

amin = 25 + a3*(4/9);
amax = 0.5 + a4*(0.25/9);

wp = 2*pi*fp;
ws = 2*pi*fs;

%% Metasximatismos prodiagrafwn

Wp = 1;
Ws = wp/ws;

%% sintelestes a, e, taksi filtrou
e = sqrt(10^(amax/10)-1);

n = (acosh(((10^(amin/10)-1)/(10^(amax/10)-1))^(1/2)))/acosh(Ws);
n = ceil(n);

a = (asinh(1/e))/n;
%% sixnotita hmiseias isxios
Whp = cosh((acosh((10^(amax/10)-1)^(-1/2)))/n);
%% Gwnies Butterworth

y1 = 0;
y2 = 36;
y3 = -36;
y4 = 72;
y5 = -72;
%% poloi protipis sinartisis
p1 = -sinh(a)*cosd(y1) + 1i*cosh(a)*sind(y1);
p2 = -sinh(a)*cosd(y2) + 1i*cosh(a)*sind(y2);
p3 = -sinh(a)*cosd(y3) + 1i*cosh(a)*sind(y3);
p4 = -sinh(a)*cosd(y4) + 1i*cosh(a)*sind(y4);
p5 = -sinh(a)*cosd(y5) + 1i*cosh(a)*sind(y5);
%% W, Q polwn Chebyshev
W1 = 0;
Q1 = 0.5;
w0_23 = sqrt((real(p2))^2+(imag(p2))^2);
Q_23 = w0_23 /(2 * abs(real(p2)));

w0_45 = sqrt((real(p4))^2+(imag(p4))^2);
Q_45 = w0_45/(2*abs(real(p4)));
%% Antistrofi polwn
whp = wp/Whp;

s1 = wp/abs(p1);
w23 = wp/w0_23;
w45 = wp/w0_45;
%% Piknotis
C = 10^(-8);
%% Proti monada
R11 = 1;
C11 = 1;
p_1 = 1/(R11*C11);
%% Klimakopoihsh
kf1 = s1;
C11new = C;
km1 = C11/(C11new*kf1);
R11new = R11*km1;
%% Deuteri monada
C21 = 1;
C22 = 1;
R21 = 1;
R22 = 1;
r21 = 1;
r22 = 2 - (1/Q_23);
k2 = 3 - (1/Q_23);
%% Klimakopoihsh
kf2 = w23;
C21new = C;
C22new = C;
km2 = C21 /(C21new*kf2);
R21new = R21*km2;
R22new = R22*km2;
r21new = r21*km2;
r22new = r22*km2;
%% Triti monada
C31 = 1;
C32 = 1;
R31 = 1;
R32 = 1;
r31 = 1;
r32 = 2 - (1/Q_45);
k3 = 3 - (1/Q_45);
%% Klhmakopoihsh
kf3 = w45;
C31new = C;
C32new = C;
km3 = C31 /(C31new*kf3);
R31new = R31*km3;
R32new = R32*km3;
r31new = r31*km3;
r32new = r32*km3;
%% Rythmisi kerdous
ktotal = k2*k3;
a_k = (10^(0.5))/ktotal;
% a<1 ara pathitiki
%% Synartiseis metaforas

T1 = tf([1 0],[1 s1]);
T2 = tf([k2 0 0],[1 w23/Q_23 w23^2]);
T3 = tf([k3 0 0],[1 w45/Q_45 w45^2]);

T_tot = a_k*T1*T2*T3;

plot_transfer_function(T1, [fp fs]);
plot_transfer_function(T2, [fp fs] );
plot_transfer_function(T3, [fp fs]);
plot_transfer_function(T_tot, [90000 fp fs]);

Invsys=inv(T_tot);
plot_transfer_function(Invsys,[90000 fp fs])

ltiview({'bodemag'}, T1, T2, T3, T_tot)
%% Rythmisi kerdous sta 0dB
a_0 = 1/ktotal;
T = a_0*T1 * T2 *T3;
Invsys2=inv(T);
plot_transfer_function(Invsys2,[90000 fp fs])
%% Fourier Analysis

Time = (1/100);
Fs = 1000000; % sample frequency
dt = 1/Fs;
t = 0:dt:Time-dt;

fin_1=(0.4*ws)/(2*pi);
fin_2=(0.9*ws)/(2*pi);
fin_3=(1.4*wp)/(2*pi);
fin_4=(2.4*wp)/(2*pi);
fin_5=(4.5*wp)/(2*pi);

% input signal
x =cos(2*pi*fin_1*t)+0.5*cos(2*pi*fin_2*t)+cos(2*pi*fin_3*t)+0.7*cos(2*pi*fin_4*t)+0.5*cos(2*pi*fin_5*t);

figure;
plot(t,x);
title('Input signal');

y=lsim(T_tot,x,t);
figure;
plot(t,y);
title('Output signal');

figure;
plot(t,x);
hold on;
plot(t,y);
hold off;
title('Input and Output signals');

Xf=fft(x);
L=length(x);
P2 = abs(Xf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.001 20000 0 inf]);
title('Single-Sided Amplitude Spectrum of of Input signal');

Yf=fft(y);
L=length(y);
P2 = abs(Yf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.001 20000 0 inf]);
title('Single-Sided Amplitude Spectrum of Output signal');