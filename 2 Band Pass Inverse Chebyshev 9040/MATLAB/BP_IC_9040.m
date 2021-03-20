%% Anastasios Mouratidis 9040
%% Band Pass Inverse Chebyshev
clear;
close all;
% AEM
a1 = 9;
a2 = 0;
a3 = 4;
a4 = 0;
%% Prodiagrafes provlimatos
f0 = 1000;
f1 = 650 + (25*a4);
f2 = (f0^2)/f1;
D  = 2.1*((f0^2)-(f1^2))/f1;
f3 = (-D + sqrt(D^2 + 4*f0^2))/2;
f4 = (f0^2)/f3;

amin = 35 - a3;
amax = 0.4 + (a4/36);

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

Wp = 1;
Ws = (w4 - w3)/(w2 - w1);

bw = w2-w1;
qc = w0/bw;
n = acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Ws);
n = ceil(n);

e = 1/sqrt(10^(amin/10)-1);
a = (1/n)*( asinh(1/e));

%% sixnotita imiseias isxuos
whp = 1/cosh(acosh(1/e)/n);
%% Gwnies
y1 = 22.5;
y2 = -22.5;
y3 = 67.5;
y4 = -67.5;
%% Poloi
p1 = -sinh(a) * cosd(y1) + 1i*cosh(a) * sind(y1);
p2 = -sinh(a) * cosd(y2) + 1i*cosh(a) * sind(y2);
p3 = -sinh(a) * cosd(y3) + 1i*cosh(a) * sind(y3);
p4 = -sinh(a) * cosd(y4) + 1i*cosh(a) * sind(y4);

W_12 = sqrt((real(p1))^2+(imag(p1))^2);
W_34 = sqrt((real(p3))^2+(imag(p3))^2);

Q_12 = W_12/(2*abs(real(p1)));
Q_34 = W_34/(2*abs(real(p3)));

%% antistrofi polon
InvW1 = 1 / W_12;
InvW2 = 1 / W_34;

%% klimakopoihsh sixnotitas
InvW1 = InvW1 * Ws;
InvW2 = InvW2 * Ws;
%% neoi poloi
S12 = InvW1 / ( 2 * Q_12 );
S34 = InvW2 / ( 2 * Q_34);
W12 = sqrt( InvW1^2 - S12^2 );
W34 = sqrt( InvW2^2 - S34^2 );
%% midenika gia k=1 kai k=3
Z1 = sec(pi/8);
Z2 = sec(3*pi/8);

% klimakopoihsh mhdenikwn
Z1new = Z1 * Ws;
Z2new = Z2 * Ws;
%% Metasximatismos migadidkou polou
C1 = S12^2 + W12^2;
D1 = (2*S12) / qc;
E1 = 4 + C1/ (qc^2);
G1 = (E1^2 - 4 * D1^2)^(1/2);
Q1 = 1/D1 * sqrt (1/2*(E1+ G1));
y1 = acosd(1/(2*Q1));
K1 = (S12*Q1) / qc;
W1 = K1 + sqrt( K1^2 -1);
wo1 = w0 / W1;
wo2 = W1 * w0;
%% Metasximatismos migadidkou polou
C2 = S34^2 + W34^2;
D2 = (2*S34) / qc;
E2 = 4 + C2/ (qc^2);
G2 = (E2^2 - 4 * D2^2)^(1/2);
Q2 = 1/D2 * sqrt (1/2*(E2+ G2));
y2 = acosd(1/(2*Q2));
K2 = (S34*Q2) / qc;
W2 = K2 + sqrt( K2^2 -1);
wo3 = w0 / W2;
wo4 = W2 * w0;
%% metasximatismos 1ou fantastikou midenikou
Kzero1 = 2 + (Z1new^2) / (qc^2);
x1 = ( Kzero1 + sqrt( Kzero1^2 - 4 ) ) / 2;
wz1 = w0 * ( sqrt(x1) );
wz2 = w0 / ( sqrt(x1) );
%% metasximatismos 2ou fantastikou midenikou 
Kzero2 = 2 + (Z2new^2) / (qc^2);
x2 = ( Kzero2+ sqrt( Kzero2^2 - 4 ) ) / 2;
wz3 = w0 * ( sqrt(x2) );
wz4 = w0 / ( sqrt(x2) );
%% Piknotis
C= 10^(-8);
%% Proti monada - LPN
wzo1 = wz1/wo1;
R11 = 1;
R12 = 4*Q1^2;
R13 = (wzo1^2)/(2*Q1^2);
R14 = 1;
R15 = (4*Q1^2)/(wzo1^2-1);
C11 = 1/(2*Q1);
k1  = 1/(1+(wzo1^2)/(2*Q1^2)); 
% Klimakopoihsh
kf1 = wo1;
km1 = C11/(kf1*C);
C11new = C;
R11new = R11*km1;
R12new = R12*km1;
R13new = R13*km1;
R14new = R14*km1;
R15new = R15*km1;
%% Deuteri monada - HPN
wzo2 = wz2/wo2;
k21 = (1/wzo2)^2 -1;
k22 = ((2+k21)*Q1^2)/((2+k21)*Q1^2 +1 );
k2  = k22*(1/wzo2)^2;
R21 = 1;
R22 = (k21+2)^2*Q1^2;
R23 = 1;
R24 = (k21+2)*Q1^2;
C_2 = 1/(Q1*(2+k21));
C21 = k21*C_2;
% Klimakopoihsh
kf2 = wo2;
km2 = C_2 / (kf2*C);
C2new  = C;
C21new = C21/(km2 * kf2);
R21new = R21 * km2;
R22new = R22 * km2;
R23new = R23 * km2;
R24new = R24 * km2;
%% Triti monada - LPN
wzo3 = wz3/wo3;
R31 = 1;
R32 = 4*Q2^2;
R33 = (wzo3^2)/(2*Q2^2);
R34 = 1;
R35 = (4*Q2^2)/(wzo3^2-1);
C31 = 1/(2*Q2);
k3  = 1/(1+(wzo3^2)/(2*Q2^2)); 
% Klimakopoihsh
kf3 = wo3;
km3 = C31/(kf3*C);
C31new = C;
R31new = R31*km3;
R32new = R32*km3;
R33new = R33*km3;
R34new = R34*km3;
R35new = R35*km3;
%% Monada 4 - HPN
wzo4 = wz4/wo4;
k41 = (1/wzo4)^2 -1;
k42 = (2+k41)*Q2^2/((2+k41)*Q2^2 +1);
k4  = k42*(1/wzo4)^2;
R41 = 1;
R42 = (k41+2)^2*Q2^2;
R43 = 1;
R44 = (k41+2)*Q2^2;
C_4 = 1/(Q2*(2+k41));
C41 = k41*C_4;
% Klimakopoihsh
kf4 = wo4;
km4 = C_4 / (kf4*C);
C4new  = C;
C41new = C41/(km4 * kf4);
R41new = R41 * km4;
R42new = R42 * km4;
R43new = R43 * km4;
R44new = R44 * km4;
%% synartiseis metaforas
T1 = tf( [k1 0 ( k1 * wz1^2 ) ], [ 1 ( wo1 / Q1 ) wo1^2 ] );
T2 = tf( [k2 0 ( k2 * wz2^2 ) ], [ 1 ( wo2 / Q1 ) wo2^2 ] );
T3 = tf( [k3 0 ( k3 * wz3^2 ) ], [ 1 ( wo3 / Q2 ) wo3^2 ] );
T4 = tf( [k4 0 ( k4 * wz4^2 ) ], [ 1 ( wo4 / Q2 ) wo4^2 ] );

%% Rythmisi kerdous
K_total = k1*k2*k3*k4;
T_tot = T1*T2*T3*T4;
gain = abs(evalfr(T_tot, w0 * 1i));

a_k = (10^(0.5))/gain;

% a_k<1
r1 = 10000;
r2 = r1 * a_k;

%% Plots
T_tot = a_k * T_tot;

plot_transfer_function(T1, [f1 f2 f3 f4])

plot_transfer_function(T2, [f1 f2 f3 f4])

plot_transfer_function(T3, [f1 f2 f3 f4])

plot_transfer_function(T4, [f1 f2 f3 f4])

ltiview({'bodemag'}, T1, T2, T3, T4, T_tot)

plot_transfer_function(T_tot, [f1 f2 f0 f3 f4])

InvSys_new = inv (T_tot);
plot_transfer_function(InvSys_new, [f1 f2 f0 f3 f4])
%% Rythmisi sta 0 dB gia tin deyteri sinartisi aposvesis
a_0 = 1/gain;
T = a_0*T1 * T2 *T3 * T4;
Invsys2=inv(T);
plot_transfer_function(Invsys2,[f1 f2 f0 f3 f4])
%% Fourier Analysis
Time = (1/100);
Fs = 1000000;
dt = 1/Fs;
t = 0:dt:Time-dt;

fin_1 = (w0 - ((w0-w1)/2)) / (2*pi);
fin_2 = (w0 + ((w0+w1)/3)) / (2*pi);
fin_3 = (0.4 * w3) / (2*pi);
fin_4 = (2.5 * w4) / (2*pi);
fin_5 = (3 * w4) / (2*pi);
% input signal
x = cos(2*pi*fin_1*t)+0.8*cos(2*pi*fin_2*t)+0.8*cos(2*pi*fin_3*t) ...
    +0.6*cos(2*pi*fin_4*t)+0.5*cos(2*pi*fin_5*t);

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
axis([0.01 20000 0 inf]);
title('Single-Sided Amplitude Spectrum of Input signal');

Yf=fft(y);
L=length(y);
P2 = abs(Yf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('Single-Sided Amplitude Spectrum of Output signal');


