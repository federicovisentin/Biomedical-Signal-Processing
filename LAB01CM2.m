% Programma per la verifica delle prestazioni dell'averaging su 100 epoche di segnale
% con e senza riallineamento delle forme d'onda in presenza di jitter distribuito
% uniformemente.

clear all
close all
pack
clc

load signorm5



xor=[s1 zeros(1,68)]/max(s1);			%forma d'onda originale

% visualizzazione di 300 campioni di segnale
%-------------------------------------------
plot(x(1:300)/max(x(1:300))), grid on,title('Primi 300 campioni di segnale'),xlabel('campioni'),ylabel('ampiezza (a.u.)')
clc,fprintf('Press a key to proceed with the spectral matching alignment\n'),pause

% spectral matching
%------------------
for i=0:98
   xact=x(101+100*i:200+100*i);
   xc=xcorr(xact,xor);					%cross-correlazione del segnale originale e della epoca i-esima
   M=max(xc);
   indM=find(xc==M)-100;
   % rifinisce allineamento
   d(i+2)=delay(real(fft(xact)),imag(fft(xact)),real(fft(xor)),imag(fft(xor)),indM); %stima del ritardo in campioni
end
d=round(abs(d));

% creazione della matrice di epoche riallineate
%----------------------------------------------
xrial(1,:)=x(1:100);
for i=0:98
   segment=x(i*100+101:i*100+200);
   if d(i+2)==0
      xrial(i+2,:)=segment;
   else
      xrial(i+2,:)=[segment(1+d(i+2):100) segment(1:d(i+2))];
  end
end

% plot delle 100 epoche sovrapposte
%----------------------------------
for j=0:99
   matx(j+1,:)=x(j*100+1:j*100+100);
end
figure,plot(matx'),grid on,xlabel('campioni'),ylabel('ampiezza'),
title('Sovrapposizione delle 100 epoche di segnale senza riallineamento'),
clc,fprintf('Press any key to proceed to the averaging\n'),pause

figure,plot(xrial'),grid on,xlabel('campioni'),ylabel('ampiezza'),
title('Sovrapposizione delle 100 epoche di segnale con riallineamento'),
clc,fprintf('Press any key to proceed to the averaging\n'),pause

figure,
subplot(2,1,1),plot(xor/max(xor)),title('Forma d''onda originale')
subplot(2,1,2),plot(mean(matx)),title('Risultante dopo averaging di 100 epoche senza riallineamento'),
clc,pause

figure,
subplot(2,1,1),plot(xor/max(xor)),title('Forma d''onda originale')
subplot(2,1,2),plot(mean(xrial)),title('Risultante dopo averaging di 100 epoche con riallineamento')


