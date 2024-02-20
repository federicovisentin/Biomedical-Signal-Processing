% Elaborazione di Segnali Biomedici
% Ottavo laboratorio: trasformate tempo - frequenza

clear all
close all
pack
clc

%segnale sinusoidale 20 Hz tra 0 e 0.4 e 50 Hz tra 0.6 e 1


fc = 256;		% sampling rate
T = 1;	% window length
asset1 = [0:1/fc:0.4];	% x axis
asset1_1 = [0.4+1/fc:1/fc:0.6];
asset2 = [0.6+1/fc:1/fc:1];
asset = [0:1/fc:1-1/fc];
fsin1 = 20;		% frequency of the sine wave
fsin2 = 50;
x1=sin(2*pi*fsin1*asset1);
x1_1=zeros(1,length(asset1_1));
x2=sin(2*pi*fsin2*asset2);
x =[x1 x1_1 x2] ;% sine wave (not analytical)
plot(x);

x=x-mean(x);

x = hilbert(x);	% use Hilbert transform to obtain analytic signal

ntlag = max(size(x))/2;		% number of time lags
sigma = 1;					% kernel parameter
w_cw = dcw(x,fc,ntlag,sigma,T);  % Time-Frequency transform (see functions dcw.m and wig.m)
w_wig = wig(x,fc,ntlag);

[roww_cw, colw_cw] = size(w_cw);
[roww_wig, colw_wig] = size(w_wig);
upper = fc/4;	% upper frequency since non alias-free implementation
assef = (0:roww_cw-1)/roww_cw*upper;			% frequency axis

figure(1)

subplot(2,1,1)
mesh(asset,assef,real(w_wig)),xlabel('time (s)'),ylabel('frequency (Hz)'),title(' WIG transform '),
axis([min(asset) max(asset) min(assef) max(assef) 2*min(min(real(w_wig))) 2*max(max(real(w_wig)))]);

subplot(2,1,2)
mesh(asset,assef,real(w_cw)),xlabel('time (s)'),ylabel('frequency (Hz)'),title(' CW transform '),
axis([min(asset) max(asset) min(assef) max(assef) 2*min(min(real(w_cw))) 2*max(max(real(w_cw)))]);

figure(2)

subplot(2,1,1)
contour(asset,assef,real(w_wig)/max(max(real(w_wig))),.1:.1:1);
xlabel('time (s)');
ylabel('frequency (Hz)');
title( 'WIG transform');

subplot(2,1,2)
contour(asset,assef,real(w_cw)/max(max(real(w_cw))),.1:.1:1);
xlabel('time (s)');
ylabel('frequency (Hz)');
title( 'CW transform');

%come plottare la funzione di ambiguità
amb=daf(x,ntlag,fc,T);
assesig=[-127:1:128];
assetau=[0:1:127];
figure(3)
contour(assesig,assetau,real(amb)/max(max(real(amb))));








