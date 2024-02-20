% Quarto laboratorio del corso di ESB
% Stima spettrale non parametrica applicata al segnale EEG.

clear all
close all
pack
clc

% caricamento dei dati
%---------------------
h = fopen("EEG2501.bin");

% Leggere l'help della funzione fread, in particolare su come si può
% definire il parametro di ingresso 'sizeA' per avere in uscita una matrice
eeg = fread(h,[25 inf,],'float');                                 
fclose(h);

% Andare a vedere nel Workspace di MATLAB le dimensioni della matrice 'eeg'
% e fare in modo che ogni colonna della matrice corrisponde ad un canale 


fc =  512;                    %definire la frequenza di campionamento
nest = 1;                 %contatore per marcare il rilascio della figura

% visualizzazione segnali nel dominio del tempo
%----------------------------------------------
chan7=eeg(7,:);
chan9=eeg(9,:);
figure(1),
subplot(2,1,1),plot(chan7(1,9*fc:10*fc)),title('Canale 7')
subplot(2,1,2),plot(chan9(9*fc:10*fc)),title('Canale 9')
xlabel('campioni'), pause
clf

% ciclo sui dati
%---------------
for jj = 7:2:9
    
    x = eeg(jj,:);
    %Rimuovere il valor medio
    x = x-mean(x);
    
    ax(1) = subplot(3,1,1); ax(2) = subplot(3,1,2); ax(3) = subplot(3,1,3); %  tre assi: uno per ogni plot (metodo di stima)
    % Calcolo periodogramma tramite periodogramma classico
    %----------------------------------------------------

    if nest==4
        hold off
    end

    %Trovare il valore NFFT sapendo che si vuole ottenere una risoluzione
    %apparente pari a 0.025Hz:
    ncamp=7680;
    NFFT = fc/0.025;
    w=boxcar(ncamp);
    [Pxx,f]=pwelch(x,w,0,NFFT,fc);
    plot(ax(1),f,Pxx/max(Pxx)); 

    title('Periodogramma classico (ris. 0.025 Hz) ') 
    axis(ax(1),[0 127 0 1]); 
    nest=nest+1;
    hold on
    
    % Suddivisione potenza media tra le 5 bande di interesse
    %-------------------------------------------------------
    %Leggere l'help della funzione FIND e trovare i range di indici per la
    %ripartizione della potenza media del segnale nelle bande di interesse
    d=find(f>0.5 & f<3.5);
    t=find(f>3.5 & f<7);
    a=find(f>7 & f<14);
    b1=find(f>14 & f<21);
    b2=find(f>21 & f<30);
    Ptot=sum(Pxx);
    P=sum(Pxx(d));
    Pd=P/Ptot;
    P=sum(Pxx(t));
    Pt=P/Ptot;
    P=sum(Pxx(a));
    Pa=P/Ptot;
    P=sum(Pxx(b1));
    Pb1=P/Ptot;
    P=sum(Pxx(b2));
    Pb2=P/Ptot;
    fprintf('\n Canale %i',jj); 
    fprintf('\n Periodogramma classico - risoluzione 0.025 Hz'); 
    fprintf('\n Pd = %f Pt = %f Pa = %f Pb1 = %f Pb2 = %f \n',Pd,Pt,Pa,Pb1,Pb2);
 
    % Ripeto conti con Welch
    %-----------------------
    
    %Trovare il valore NFFT sapendo che si vuole ottenere una risoluzione
    %apparente pari a quella teorica:
    w=hamming(1024);
    NFFT = 512/0.5; % Nota: RisSp = df := 1/T = 1/(Ndt) = fc/N => N = fc/df;  
    
    [Pxx,f]=pwelch(x,w,512,NFFT,fc);
    plot(ax(2),f,Pxx/max(Pxx),'r');
    title('Metodo di Welch (Ris. 0.5 Hz) ')
    axis(ax(2),[0 127 0 1]); 
    nest=nest+1;
    
    % Suddivisione potenza media tra le 5 bande di interesse
    %-------------------------------------------------------
    %Leggere l'help della funzione FIND e trovare i range di indici per la
    %ripartizione della potenza media del segnale nelle bande di interesse
    d=find(f>0.5 & f<3.5);
    t=find(f>3.5 & f<7);
    a=find(f>7 & f<14);
    b1=find(f>14 & f<21);
    b2=find(f>21 & f<30);
    Ptot=sum(Pxx);
    P=sum(Pxx(d));
    Pd=P/Ptot;
    P=sum(Pxx(t));
    Pt=P/Ptot;
    P=sum(Pxx(a));
    Pa=P/Ptot;
    P=sum(Pxx(b1));
    Pb1=P/Ptot;
    P=sum(Pxx(b2));
    Pb2=P/Ptot;
    fprintf('\n Metodo di Welch - risoluzione 0.5 Hz') 
    fprintf('\n Pd = %f Pt = %f Pa = %f Pb1 = %f Pb2 = %f \n',Pd,Pt,Pa,Pb1,Pb2);
    % Ripeto conti con correlogramma
    %-------------------------------
    % Trovare il numero di ritardi da considerare per ottenere una sequenza
    % di autocorrelazione di almeno 512 valori
    ntlag = 256;
    
    acs= xcorr(x,ntlag,'biased');   
   
    w1=hamming(512); 
    w1=w1';
    %Prendere il pezzo centrale della sequenza di autocorrelazione della
    %lunghezza della finestra considerata 
    acs1 = acs(1:512);
    %calcolare la ACS finestrata
    acsw=acs1.*w1;
   
    %Calcolare la risoluzione teorica:
    ris = 1; %perchè la lunghezza della finestra è 512 camp
    
    %Trovare il valore NFFT sapendo che si vuole ottenere una risoluzione
    %apparente pari a quella teorica:
    NFFT = 512;    
    
    %stima spettrale
    Pxxc= fft(acsw,NFFT);
    %eliminare la replica spettrale
    Pxxc = Pxxc(1:length(Pxxc)/2+1)
    
    f=512*(0:256)/512;
    plot(ax(3),f,abs(Pxxc/max(abs(Pxxc))),'g'); 
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);
    axis(ax(3),[0 127 0 1]); 
    nest=nest+1;
    
    % Suddivisione potenza media tra le 5 bande di interesse
    %-------------------------------------------------------
    %Leggere l'help della funzione FIND e trovare i range di indici per la
    %ripartizione della potenza media del segnale nelle bande di interesse
        
    d=find(f>0.5 & f<3.5);
    t=find(f>3.5 & f<7);
    a=find(f>7 & f<14);
    b1=find(f>14 & f<21);
    b2=find(f>21 & f<30);
    Ptot=sum(abs(Pxxc));
    P=sum(abs(Pxxc(d)));
    Pd=P/Ptot;
    P=sum(abs(Pxxc(t)));
    Pt=P/Ptot;
    P=sum(abs(Pxxc(a)));
    Pa=P/Ptot;
    P=sum(abs(Pxxc(b1)));
    Pb1=P/Ptot;
    P=sum(abs(Pxxc(b2)));
    Pb2=P/Ptot;
    fprintf('\n Correlogramma - risoluzione 1 Hz') 
    fprintf('\n Pd = %f Pt = %f Pa = %f Pb1 = %f Pb2 = %f \n',Pd,Pt,Pa,Pb1,Pb2); linkaxes(ax,'x'), pause, clf % Studiare sull'help la funzione linkaxes
end
hold off

% Mappa valore efficace
%----------------------
% Costruzione della matrice
Z=[std(eeg(:,1:5)); std(eeg(:,6:10)); std(eeg(:,11:15)); std(eeg(:,16:20)); std(eeg(:,21:25))]';

% Interpolazione
mappa(Z,'Mappa distribuzione valore efficace segnale EEG');
pause;

% Mappe potenze relative ritmi diversi
for i=1:5
    for j=1:5
        x2=eeg(:,j+(i-1)*5);
        [Zd(i,j),Zt(i,j),Za(i,j),Zb1(i,j),Zb2(i,j),Zu(i,j)]=rel_pot(x2);
    end
end

mappa(Zd,'Mappa ritmo delta');
pause;
mappa(Zt,'Mappa ritmo teta');
pause;
mappa(Za,'Mappa ritmo alfa');
pause;
mappa(Zb1,'Mappa ritmo beta 1');
pause;
mappa(Zb2,'Mappa ritmo beta 2');
