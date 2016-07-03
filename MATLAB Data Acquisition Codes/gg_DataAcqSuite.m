%% DATA ACQUISITION SESSION CODES

% Send And Read Data from NIDAQ Analog Channels- AI-AO
% set 
addpath('C:\Documents and Settings\ceylana\Desktop\DASessions\MASTER CODES')

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS OF DATA ACQ SESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sessionName='NIDAQOUT-Shaker-Accel-NIDAQIN-withPleksiBeam';

duration=3;
stopsending=0;

samplingRate=10000;
sr=samplingRate;
AMP=1;                                        %VDC 
freqIN=100;                                  %Hz gets in sine as rad/sec

T=linspace(0,duration,duration*samplingRate);
sigIN=AMP*sin(freqIN*2*pi*T);

ACCGAIN=2.6;

low=0.1/sr;
high=freqIN/sr*1.1;

citylow=49.9/sr;
cityhigh=50.01/sr;

% ADDED Titles information for Plots 
titlesforPlots=(['Inputs to System ---> Amplitude: ',num2str(AMP),', Frequency: ', num2str(freqIN), ', Duration: ', num2str(duration), ' seconds']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tempstr1=pwd;
if(~strcmp(tempstr1,'C:\Documents and Settings\ceylana\Desktop\DASessions'))
cd 'C:\Documents and Settings\ceylana\Desktop\DASessions'
end

tempstrtosave=['ma_',sessionName,' ',strrep(datestr(now),':','-'),'.mat'];
tempstr=['ma_',sessionName,' ',strrep(datestr(now),':','-')];

%% Standart NIDAQ Analog Input Codes

AI= analoginput('nidaq','Dev1');   % Analog Input Activation
chan = addchannel(AI,0,{'Laser'}); % Prepare Channel

inputActualRate = setverify(AI,'SampleRate',samplingRate); % Sample Rate Setting and Verify
samplePerTrigger=inputActualRate*duration;  % Trigger Setting
set(AI,'InputType','SingleEnded')
set(AI,'SampleRate',inputActualRate);
set(AI,'SamplesPerTrigger',samplePerTrigger);
set(AI,'TransferMode','Interrupts');

%% Standart NIDAQ Analog Output Codes

AO0 = analogoutput('nidaq','Dev1');
addchannel(AO0,0);

set(AO0,'SampleRate',samplingRate);
set(AO0,'TriggerType','Immediate');
outputActualRate = get(AO0,'SampleRate');

%% Data Acquisiton

putdata(AO0,sigIN');
start([AI AO0]);

% trigger(AO0);
% wait(AI, duration+1);

inputSignalData=getdata(AI);

save( tempstrtosave);

delete(AO0);
clear AO0;

delete(AI);
clear AI;

%% RENAME DATA SHORTER %%

SOUT=inputSignalData;
SIN=sigIN;
n=length(T);

%% PROESSING - FILTERS - FFT
SOUT=ACCGAIN*SOUT;


[b,a]=butter(1,[low high ],'bandpass');
[c,d]=butter(1,[citylow,cityhigh],'stop');

% SIN=filter(c,d,[filter(b,a,SIN)]);
% SOUT=filter(c,d,[filter(b,a,SOUT)]);

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


