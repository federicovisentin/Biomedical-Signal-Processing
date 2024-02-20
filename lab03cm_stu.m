% Terzo Laboratorio di Elaborazione di Segnali Biomedici.
% Stima della ACS e stima spettrale basata sul correlogramma.

clear all
close all
pack
clc

% generazione di una sequenza numerica di prova
%----------------------------------------------
x = [1:10];
ntlag =5;
bias =0;
xacs_mine = acs(x,ntlag,bias);
xacs_mat = xcorr(x,ntlag,'unbiased');       %cambiare in 'biased' se necessario

figure(1),plot(xacs_mine,'b'),hold on,plot(xacs_mat,'r')
err = sum(xacs_mine - xacs_mat);
str = sprintf('L''errore tra le due stime è pari a %2.2f',err);
title(str);

% determinazione della ACS di un segnale deterministico
%------------------------------------------------------
fc = 100;
T = 1;
t = [0:1/fc:T];
f0 = 15;
x = sin(2*pi*f0*t); %ottenere qui la sinusoide richiesta
ntlag = round(length(x)/10);                       %cambiare se necessario
xacs_mine = acs(x,ntlag,bias);
xacs_mat = xcorr(x,ntlag,'unbiased');       %cambiare in 'biased' se necessario
figure(2),
subplot(211),plot(xacs_mine),title('My routine')
subplot(212),plot(xacs_mat),title('Matlab routine')

% determinazione della ACS di un segnale stocastico
%--------------------------------------------------
x = randn(1,1024);%ottenere qui un vettore di dati casuali
x=x-mean(x);
ntlag = round(length(x)/10);                %cambiare se necessario
xacs_mine = acs(x,ntlag,bias);
xacs_mat = xcorr(x,ntlag,'unbiased');       %cambiare in 'biased' se necessario
figure(3),
subplot(211),plot(xacs_mine),title('My routine')
subplot(212),plot(xacs_mat),title('Matlab routine')

% stima spettrale mediante ACS
%-----------------------------
NFFT = 512;                             %numero di punti su cui calcolare la FFT
XX = fft(xacs_mine,NFFT);                    %ottenere la FFT de segnale
XX = XX(1:length(XX)/2+1);              %elimina la replica spettrale
f = (0:length(XX)-1)/length(XX)*(fc/2); %creazione dell'asse delle frequenze
figure(4),plot(f,abs(XX))               %plot dei risultati
xlabel('frequency (Hz)'),ylabel('FFT absolute value (a.u.)')
axis([0 50 0 50]);
XXX=fft(x,NFFT);
YYY=XXX.^2;
YYY = YYY(1:length(YYY)/2+1)
hold on 
plot(f,abs(YYY),'r')
