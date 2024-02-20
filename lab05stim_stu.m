% Quinto laboratorio del corso di Elaborazione Segnali Biomedici

 clear all
 pack
 clc
 close all

%-----------------------------------------------------------------------
% Selezione del file sul quale lavorare
%-----------------------------------------------------------------------
 [file,path]=uigetfile('*.sig','Load');
 filename=sprintf('%s%s',path,file);

%-----------------------------------------------------------------------
% Apertura file e lettura dati nella matrice XX
%-----------------------------------------------------------------------
 h = fopen(file,'r');
 XX = fread(h,[4 inf],'short');                  %I dati nei file sono di formato "short"
 fclose(h);
 XX=XX';

%**********     la matrice XX contiene i dati demultiplexati     *******
%********** col 1 = SD; col 2 = DD1; col 3 = DD2, col 4 = coppia *******

%-----------------------------------------------------------------------
% Dati aggiuntivi
%-----------------------------------------------------------------------
 fsamp = 1024;
 fstim = input('Type the stimulation frequency (Hz):   ');
 nsamp = fix(fsamp/fstim);
 int_dist = 10e-3;

 
%-----------------------------------------------------------------------
% Visualizza l'andamento dei segnali su 0.25 secondi
%-----------------------------------------------------------------------
figure(1),
 set(gcf,'Color',[1 1 1]);
 subplot(4,1,1), plot(XX((nsamp*fstim):(fix(nsamp*fstim*1.5)),1),'b'); 
 title('Single differential'); axis([1 (fix(nsamp*fstim*0.25)) -25000 25000]);axis('off'); 
 subplot(4,1,2), plot(XX((nsamp*fstim):(fix(nsamp*fstim*1.5)),2),'g'); 
 title('Double differential 1'); axis([1 (fix(nsamp*fstim*0.25)) -25000 25000]); axis('off'); 
 subplot(4,1,3), plot(XX((nsamp*fstim):(fix(nsamp*fstim*1.5)),3),'r'); 
 title('Double differential 2'); axis([1 (fix(nsamp*fstim*0.25)) -25000 25000]); axis('off'); 
 subplot(4,1,4), plot(XX((nsamp*fstim):(fix(nsamp*fstim*1.5)),4),'c'); 
 title('Torque'); axis([1 (fix(nsamp*fstim*0.25)) 0 30000]); 
 hold on,axis('off');h=axis;
 axes;
 set(gca,'Box','off','Xlim',[h(1) h(2)],'Ylim',[h(3) h(4)],'XGrid','on','YTick',[],'Color','none');
 xlabel('sample'),ylabel(file)
 
  fprintf('\n press a key to continue \n\n');
 pause

%------------------------------------------------------------------------
% Calcola il numero delle epoche
%------------------------------------------------------------------------
 epoch_len=1;       % Lunghezza dell'epoca IN SECONDI

 %Trovare il numero di epoche (Ricordare che deve essere un numero intero!)
 n_epoch = fix(length(XX)/fsamp);         
 
 %Inizializzare a vettori di zeri con n_epoch righe e 1 colonna
 fmeanv = zeros(epoch_len,1);
 rmsv = zeros(epoch_len,1);
 cv = zeros(epoch_len,1);

 start=2.5;

 fprintf('\n');
 for i=1:n_epoch
     
   fprintf('\n processing epoch # %d',i);
   xm = zeros(nsamp,1);
   x = zeros(1024,1);
   xd1m = zeros(nsamp,1);
   xd1 = zeros(1024,1);
   xd2m = zeros(nsamp,1);
   xd2 = zeros(1024,1);
   %--------------------------
   % calcolo dell'onda M media
   %--------------------------
   for jj=1:fstim
     xm = xm+XX((( (fstim*nsamp*(i-1))+(jj-1)*nsamp+1):((fstim*nsamp*(i-1))+jj*nsamp)),1);   
     xd1m = xd1m+XX((( (fstim*nsamp*(i-1))+(jj-1)*nsamp+1):((fstim*nsamp*(i-1))+jj*nsamp)),2);     
     xd2m = xd2m+XX((( (fstim*nsamp*(i-1))+(jj-1)*nsamp+1):((fstim*nsamp*(i-1))+jj*nsamp)),3);
   end
   
   xm = xm/fstim;
   xd1m = xd1m/fstim;
   xd2m = xd2m/fstim;
   
   x(1:nsamp) = xm;
   xd1(1:nsamp) = xd1m;
   xd2(1:nsamp) = xd2m;

   %---------------------------------------------------
   % ATTENZIONE!!! COSTRUIRE LA FUNZIONE PER IL CALCOLO 
   %               DELLA FREQUENZA MEDIA
   %---------------------------------------------------
   fmeanv(i) = fmean(x,fsamp,epoch_len);

   rmsv(i) = 5/32767*sqrt((xm'*xm)/nsamp);
   fft1r=real(fft(xd1));
   fft1i=imag(fft(xd1));
   fft2r=real(fft(xd2));
   fft2i=imag(fft(xd2));
   
   %Richiamare la funzione delay per trovare il tempo (IN CAMPIONI)
   tmp = delay(fft1r,fft1i,fft2r,fft2i,start);
   
   %ottenere qui la velocità di conduzione (m/s)
   cv(i) = int_dist/(tmp/fsamp);
   
   if (i==1)
     input_stft=x; 
   else
     input_stft=[input_stft; x];
   end
   
 end

 figure(2),
 ticks=[0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9];
 axis('off');
 axes('Xlim',[0 n_epoch],'Ylim',[0.5 2.0],'YTick',ticks,'YGrid','on');hold on
 plot(fmeanv/(fmeanv(1)),'r-');
 axis([0 30 0.5 2.0])
 hold on
 plot(rmsv/rmsv(1),'g--');
 hold on
 plot(cv/cv(1),'b-.');
 hold off
 axis([0 n_epoch 0.5 2.0]); 
 label=sprintf('Fatigue plot %s - rms = green, mnf = red, cv = blue',file);
 title(label),xlabel('time (s)');
 
%---------------------------------------------------------------------
% NON TOCCARE LE RIGHE SOTTOSTANTI
% Viene effettuata qui una stima spettrale mediante tempo-frequenza
% per osservare le modificazioni dello spettro del segnale in funzione
% del tempo
%---------------------------------------------------------------------
 fprintf('\n');
 fprintf('\n press a key to continue \n\n');
 pause;
 fprintf('\n');
 fprintf('\n STFT estimation \n');
 nfft=1024;
 overlap=0;
 [STFT,f,t]=specgram(input_stft,nfft,fsamp,boxcar(nfft),overlap);

 fprintf('\n');
 fprintf('\n SPEC estimation \n');
 SPEC=STFT.*conj(STFT);
 SPEC=SPEC/max(max(SPEC));
 [rSPEC,cSPEC]=size(SPEC);
 for i=1:cSPEC
   modSPEC(:,i)=SPEC(:,i)/max(SPEC(:,i));
 end
 figure(3),
 [xi,yi]=meshgrid(1:1:30,0:5:300);
 zi=interp2(t,f,SPEC,xi,yi,'nearest');
 pcolor(xi,yi,zi);
 shading interp
 colormap(hot(32));
 axis([min(t) max(t) 0 300]),
 xlabel('time (s)'),ylabel('frequency (Hz)');
 label=sprintf('PSD %s',file);
 title(label);

 figure(4),
 surfl(xi,yi,zi,[300,20,2])
 colormap(bone(32)),brighten(0.4),
 axis([min(t) max(t) 0 300 0 1.2]),view(70,-55),
 xlabel('time (s)'),ylabel('frequency (Hz)'),
 zlabel('relative PSD'),title(file)
