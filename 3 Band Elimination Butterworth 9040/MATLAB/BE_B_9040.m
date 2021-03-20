%% Anastasios Mouratidis 9040
%% Band Elimination Butterworth
clear;
close all;

a1 = 9;
a2 = 0;
a3 = 4;
a4 = 0;

f0 = 1.75*1000;
f1 = 1000 + (50*a4);
f2 = (f0^2)/f1;
D  = (1/(3.5*f1))*((f0^2)-(f1^2));
f3 = 0.5*(-D + sqrt((D^2) + (4*f0^2)));
f4 = (f0^2)/f3;

amin = 32  + (5/9)*a3;
amax = 0.4 + (0.25/9)*a4;

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

bw = w2 - w1;
Wp = 1;
Ws = (w2-w1)/(w4-w3);

n = (log((10^(amin/10)-1)/((10^(amax/10)-1))))/(2*log(Ws/Wp));
n = ceil(n);

Wo = Wp /(((10^(amax/10)-1))^(1/(2*n)));

wz = w0;
%% Poloi ths LP kai ths HP
p1 = -1;
p2 = -0.8090 + 1i * 0.5877;
p3 = -0.8090 - 1i * 0.5877;
p4 = -0.3090 + 1i * 0.9510;
p5 = -0.3090 - 1i * 0.9510;
%% Kanonikopoihsh
Wo_hat = 1 / Wo;

p1_hat = p1 * Wo_hat;
p2_hat = p2 * Wo_hat;
p3_hat = p3 * Wo_hat;
p4_hat = p4 * Wo_hat;
p5_hat = p5 * Wo_hat;
%% Metasximatismos pragmatikou polou -0.7915
S1 = abs(p1_hat);
qc = w0/bw;
Q1 = qc/S1;
y1 = acosd(1/(2*Q1));
wo1 = w0;

%% Metasximatismos migadidkou polou  -0.6403 +- 0.4652i
S2 = abs(real(p2_hat));
W2 = abs(imag(p2_hat));
C2 = S2^2 + W2^2;
D2 = (2*S2) / qc;
E2 = 4 + C2/ (qc^2);
G2 = (E2^2 - 4 * D2^2)^(1/2);
Q2 = 1/D2 * sqrt (1/2*(E2+ G2));
y2 = acosd(1/(2*Q2));
K2 = (S2*Q2) / qc;
W_2 = K2 + sqrt( K2^2 -1);
wo2 = W_2 * w0;
wo3 = w0 / W_2;

%% Metasximatismos migadikou polou  -0.2446 +- 0.7527i
S3 = abs(real(p4_hat));
W3 = abs(imag(p4_hat));
C3 = S3^2 + W3^2;
D3 = (2*S3) / qc;
E3 = 4 + C3/ (qc^2);
G3 = (E3^2 - 4 * D3^2)^(1/2);
Q3 = 1/D3 * sqrt (1/2*(E3 + G3));
y3 = acosd(1/(2*Q3));
K3 = (S3*Q3) / qc;
W_3 = K3 + sqrt( K3^2 -1);
wo4 = W_3 * w0;
wo5 = w0 / W_3;
%% Piknotis
C = 10^(-8);
%% Monada 1 - Boctor High Pass
% wzo2 =  wz /wo2;
BoctorHighPass(wz, wo2, Q2,100,10^(-8)) % timi eisodou C=10^(-8)
R11 = circuit.R_1; 
R12 = circuit.R_2; 
R13 = circuit.R_3; 
R14 = circuit.R_4; 
R15 = circuit.R_5; 
R16 = circuit.R_6;
C11 = circuit.C_1;
C12 = circuit.C_2;
k1  = circuit.H;
%% Monada 2 - Notch
% wo = 1 ara wz = wo1 = 1
k21 = 0;
k22 = ((2+k21)*Q1^2)/((2+k21)*Q1^2 +1 );
k2  = k22;
R21 = 1;
R22 = 4*Q1^2;
R23 = 1;
R24 = 2*Q1^2;
C_2  = 1/(Q1*(2+k21));
%% Klimakopoihsh
kf2 = wo1;
km2 = C_2 / (kf2*C);
C2new  = C;
R21new = R21 * km2;
R22new = R22 * km2;
R23new = R23 * km2;
R24new = R24 * km2;

%% Triti monada
%  0.5599 < k31 < 1 kai epilegw k31 = 0.9
wzo3 =  wz /wo3;
k31 = 0.9;
R31 = 2/((k31*wzo3^2)-1);
R32 = 1/(1-k31);
R33 = (1/2)*((k31/(Q2^2))+(k31*wzo3^2)-1);
R34 = 1/k31;
R35 = 1;
R36 = 1;
C31 = k31/(2*Q2);
C32 = 2*Q2;
k3 = 1/(0.5*(k31/Q2^2+k31*wzo3^2+1));
%% Klimakopoihsh
kf3 = wo3;
km3 = C31 / (kf3*C);
C31new = C;
C32new = C32/(km3 * kf3);
R31new = R31 * km3;
R32new = R32 * km3;
R33new = R33 * km3;
R34new = R34 * km3;
R35new = R35 * km3;
R36new = R36 * km3;
%% Monada 4 - HPN
% kanonika ithele Boctor alla den isxiei i sinthiki
wzo4 = wz/wo4;
k41 = (1/wzo4)^2 -1;
k42 = (2+k41)*Q3^2/((2+k41)*Q3^2 +1 );
k4  = k42*(1/wzo4)^2;
R41 = 1;
R42 = (k41+2)^2*Q3^2;
R43 = 1;
R44 = (k41+2)*Q3^2;
C_4 = 1/(Q3*(2+k41));
C41 = k41*C_4;
%% Klimakopoihsh
kf4 = wo4;
km4 = C_4 / (kf4*C);
C4new  = C;
C41new = C41/(km4 * kf4);
R41new = R41 * km4;
R42new = R42 * km4;
R43new = R43 * km4;
R44new = R44 * km4;
%% Pempti monada
% 0.42 < k51 < 1 kai epilegw k51 = 0.9
wzo5 =  wz /wo5;
k51 = 0.9;
R51 = 2/((k51*wzo5^2)-1);
R52 = 1/(1-k51);
R53 = (1/2)*((k51/(Q3^2))+(k51*wzo5^2)-1);
R54 = 1/k51;
R55 = 1;
R56 = 1;
C51 = k51/(2*Q3);
C52 = 2*Q3;
k5 = 1/(0.5*(k51/Q3^2+k51*wzo5^2+1));
%% Klimakopoihsh
kf5 = wo5;
km5 = C51 / (kf5*C);
C51new = C;
C52new = C52/(km5 * kf5);
R51new = R51 * km5;
R52new = R52 * km5;
R53new = R53 * km5;
R54new = R54 * km5;
R55new = R55 * km5;
R56new = R56 * km5;
%% Rythmisi kerdous sta 10dB
ktotal = k1*k2*k3*k4*k5;
a_k = (10^(0.5))/ktotal;
%% Synartiseis Metaforas

T1 = tf([k1 0 k1*(wz^2)],[1 wo2/Q2 wo2^2]);
T2 = tf([k2 0 k2*(wz^2)],[1 wo1/Q1 wo1^2]);
T3 = tf([k3 0 k3*(wz^2)],[1 wo3/Q2 wo3^2]);
T4 = tf([k4 0 k4*(wz^2)],[1 wo4/Q3 wo4^2]);
T5 = tf([k5 0 k5*(wz^2)],[1 wo5/Q3 wo5^2]);


T_tot = a_k*T1*T2*T3*T4*T5;

plot_transfer_function(T1, [f1 f2 f3 f4])
plot_transfer_function(T2, [f1 f2 f3 f4])
plot_transfer_function(T3, [f1 f2 f3 f4])
plot_transfer_function(T4, [f1 f2 f3 f4])
plot_transfer_function(T5, [f1 f2 f3 f4])
plot_transfer_function(T_tot, [100 f1 f2 f3 f4])

ltiview({'bodemag'}, T1, T2, T3, T4,T5, T_tot)

InvSys_new = inv (T_tot);
plot_transfer_function(InvSys_new, [100 f1 f2 f3 f4])
%% Rythmisi sta 0 dB gia tin deyteri sinartisi aposvesis
a_0 = 1/ktotal;
T = a_0*T1*T2*T3*T4*T5;
Invsys2=inv(T);
plot_transfer_function(Invsys2,[100 f1 f2 f3 f4])
%% Fourier Analysis
Time = (1/100);
Fs = 1000000; % sample frequency
dt = 1/Fs;
t = 0:dt:Time-dt;

fin_1 = (w0 - ((w0-w3)/3)) / (2*pi);
fin_2 = (w0 + ((w0+w3)/4)) / (2*pi);
fin_3 = (0.4 * w1) / (2*pi);
fin_4 = (2.5 * w2) / (2*pi);
fin_5 = (3 * w2) / (2*pi);
% input signal
x = 0.5*cos(2*pi*fin_1*t)+0.8*cos(2*pi*fin_2*t)+0.8*cos(2*pi*fin_3*t)+0.6*cos(2*pi*fin_4*t)+1.2*cos(2*pi*fin_5*t);

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