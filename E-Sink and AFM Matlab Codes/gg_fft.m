function [f,amp,phase,rawFFT] = gg_fft(T,X,shouldTrans)
% function [f,amp,phase,rawFFT] = gg_fft(T,X,shouldTrans)
% Returns FFT in hertz and normalized form
% This function is designed to take fft of multiple signal at once
% If row and/or column formation is reversed than send shouldTrans
% parameter as 1 aka 'TRUE'
% [f,amp,phase,rawFFT] = gg_fft(T,X)
% f     - Frequency vec
% amp   - Relative Amp |G(w)|
% phase - Phi (rad)
% rawFFT- Raw data a+bi form FFT;

% clean Vectors from NaNs;

if nargin<3,shouldTrans=0;end

for i=1:length(X);
if(isnan(X(i)));
X(i)=0;
end
end

for i=1:length(T);
if(isnan(T(i)));
X(i)=0;
end
end

if (shouldTrans);
    T=transpose(T);
    X=transpose(X);
end;

% clean Vectors from NaNs;
 
MAX=max(T);
N=2^nextpow2(length(T));
dt=MAX/N;
df=1/MAX;

tempf=0:df:(1/dt)-df;
f=transpose(tempf(1:length(tempf)/2));

tempt=0:dt:MAX;

vecsize=size(X);

for i=1:vecsize(1);
tempx=interp1(T,X(i,:),tempt);
rawFFT(i,:)=fft(tempx,N)./N*2;
tempraw=rawFFT(i,:);
amp(i,:)=abs(tempraw(1:length(f)));
phase(i,:)=unwrap( angle( tempraw(1:length(f)) ) );
end;


 