%% RENAME DATA SHORTER %%

SOUT=inputSignalData;
SIN=sigIN';
n=length(T);

%% PROESSING - FILTERS - FFT
ACCGAIN=2.6;


low=0.1/sr;
high=freqIN/sr*1.1;

citylow=49.9/sr;
cityhigh=50.01/sr;

SOUT=ACCGAIN*SOUT;

[b,a]=butter(1,[low high ],'bandpass');
% [c,d]=butter(1,[citylow,cityhigh],'stop');

% H=hann(length(T));
% SIN=SIN.*H;
% SOUT=SOUT.*H;

SIN=filter(b,a,SIN);
SOUT=filter(b,a,SOUT);

% SIN=filter(c,d,SIN);
% SOUT=filter(c,d,SOUT);

[F_SIN,        FFT_SIN,     PHASE_SIN,     RAW_SIN]=gg_fft(T,SIN);
[F_SOUT,    FFT_SOUT, PHASE_SOUT, RAW_SOUT]=gg_fft(T,SOUT);

%% PLOTS FOR  RAW VDC IN-OUT

figure;
hold;
plot(T,SIN);
plot(T,SOUT,'--r')
hold;
title({['SIGNALS IN TIME DOMAIN'],[titlesforPlots]});                                  % ADDED Titles for Plots
xlabel('Time (sec)')
ylabel('Amp (VDC)')
legend('SIGNAL SENT', 'SIGNAL RECIEVED')
% set(gca,'XLim',[10*1/freqIN,20*1/freqIN]);
gg_ZoomAndGrid2D;
pause();

figure;
p3=FFT_SIN(1:length(F_SIN));
p4=FFT_SOUT(1:length(F_SOUT));
loglog(F_SIN,p3');
hold;
loglog(F_SIN,p4','--r');
hold;
% set(gca,'XLim',[0.5*freqIN,1.5*freqIN]);
title({['SIGNALS IN FREQUENCY DOMAIN FILTERED'],[titlesforPlots]});  % ADDED Titles for Plots
xlabel('Freq (Hz)')
ylabel('Amp |G(w)|')
legend('SIGNAL SENT', 'SIGNAL RECIEVED')
% gg_ZoomAndGrid2D;
pause();
