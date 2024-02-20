% LAB01CM_stu.M
% Primo laboratorio del corso di Elaborazione Segnali Biomedici

 close all
 clear;
 pack;
 clc;

 fprintf('\n Programma per l''elaborazione di potenziali evocati - Averaging \n\n');
 fprintf(' Press a key to choose an evpotxxx.bin file \n\n');
 pause

%-------------------
% Apertura file dati
%-------------------
 [file,path]=uigetfile('*.bin','Load');
 filename=sprintf('%s%s',path,file);
 h=fopen(filename,'r');
 x=fread(h,inf,'float');
 fclose(h);
 
%---------------------------------------------------
% Apertura file pot5100.bin (segnale di riferimento)
%---------------------------------------------------
 h=fopen('pot5100.bin','r');
 xor=fread(h,inf,'float');
 fclose(h);
 %Dal workspace di MATLAB, prendere nota delle dimensioni dei vettori x e
 %xor 

 subplot(2,1,1)
 plot(xor(1:300),'r'); %Plottare tre epoche del segnale originale, dal 1° al 300° campione
 axis([0 300 -0.3 0.3]);
 title('Synthesized potential')

 %--------------------
 % Dal segnale di riferimento, determinare i range (in campioni)
 % entro i quali calcolare:
 % - il segnale picco-picco 
 % - il rumore
 %--------------------

 range_pp = [1:50];
 range_rumore = [32:100];
 
 signal_pp=max(xor(range_pp))-min(xor(range_pp));
 noise_stdev=std(x(range_rumore));
 
 %--------------------
 % Calcolare il rapporto segnale-rumore
 %--------------------
 
 SNR = 20*log10(signal_pp/(4*noise_stdev));
 
 subplot(2,1,2)
 plot(x(1:300),'g');
 axis([0 300 -0.3 0.3]);
 title(['Corrupted potential - SNR = ', num2str(SNR,'%4.1f')])

 
 
 clc
 fprintf('\n press a key to continue with averaging \n\n');
 pause

%--------------------
% Parametri averaging
%--------------------

 %Determinare il numero di campioni di ogni epoca e il numero di epoche
 N =  100;
 n_epoch = fix(length(x)/(N));

 FR = n_epoch;		% calcola il numero delle epoche

 %Initializzazione del vettore contentente il segnale mediato. NON
 %MODIFICARE
 xav = zeros(N,1);

 %clc
 subplot(1,1,1)

%--------------------
% Calcolo della media
%--------------------
  for jj=1:FR
    
    %---------------------------------------------------------------------------------------------------------------- 
    % inserire qui sotto il calcolo della media dei brani di lunghezza 100
    % elementi 
    %
    xav = xav+x(100*jj-99:100*jj);  
    xav = xav./jj;
    %
    %---------------------------------------------------------------------------------------------------------------- 
 
    h=gcf;
    clf
    
    noise_stdev_av = std(xav(range_rumore));
    SNR = 20*log10(signal_pp/(4*noise_stdev_av));
    
    %----------
    % Aggiornare due vettori che conterranno i valori sperimentali e teorici in
    % funzione del numero di epoche considerate
    %----------
    SNR_teor(jj,1)=20*log10(signal_pp/(4*noise_stdev)*sqrt(jj));
    SNR_sper(jj,1)=SNR;
    
    
    subplot(2,1,1)
    plot(xor(1:100),'r')
    axis([0 100 -.3 .3]) 
    title('Synthesized potential')
    subplot(2,1,2)
    plot(xav(1:100),'g')
    axis([0 100 -0.3 0.3]) 
    title(['Averaged potential - Epoch ', int2str(jj), '    SNR = ', num2str(SNR,'%4.1f')])
    figure(h)    
    home
    %pause(0.1)
    drawnow
    t1=clock;
    t2=t1;
    e=etime(t2,t1);
    while e <= 0.001
    t2=clock;	
    e=etime(t2,t1);
    end
    %
    % ---------------------------------------------------------------------
    % inserire denormalizzazione dell'ampiezza del potenziale mediato
    % a cui sommare l'epoca successiva alla prossima iterazione del ciclo
    xav= xav.*jj;
    % ---------------------------------------------------------------------
    
  end
  
    % ---------------------------------------------------------------------
    % Inserire qui sotto denormalizzazione finale per ricostruire
    % l'ampiezza
    xav= xav./130;   
    % ---------------------------------------------------------------------
    

    % ---------------------------------------------------------------------
    % Inserire qui sotto il codice necessario per rappresentare
    % l'andamento di: 1) rapporto segnale-rumore teorico in funzione
    % delle epoche mediate, 2) rapporto segnale-rumore misurato
    % sperimentalmente
    % ---------------------------------------------------------------------
figure
plot(SNR_teor,'r');
hold on
plot(SNR_sper,'b');
title('SNR teorico(r) SNR misurato(b)')


    
    % ---------------------------------------------------------------------
    % Inserire qui sotto il codice necessario per creare una matrice che
    % avrà un numero di righe pari al numero di campioni di ogni epoca, e 
    % un numero di colonne pari al numero di epoche. Leggere l'help
    % della funzione reshape.
    % 
    % Fare la trasposta della matrice ottenuta in modo da ottenere le
    % epoche sulle righe.
    % 
    % Calcolare il segnale mediato totale con la funzione mean e
    % confrontare il SNR ottenuto con quello ottenuto all'ultima iterazione
    % del ciclo precedente
    % ---------------------------------------------------------------------
mat_av=reshape(x,N,n_epoch);
mat_av=mat_av';
xav_mat=mean(mat_av);
SNR_mat=20*log10(signal_pp/(4*std(xav_mat(range_rumore))));



