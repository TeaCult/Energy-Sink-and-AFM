% Satelitte to Satelitte / Satelitte to Satelitte  / Primary to Primary

% Adapting old simulations to new variables 
 if(~exist('simNo','var'))
simNo='';
 end

tForce=T;
extForce=subs(ForceSym,T);

%% EigenValue Calculation 
[U,W]=eig(K,M,'chol');  % cholesky gives alreadry mass normalized
%% Calculate Energies
XX(length(X),n)=0;
VV(length(X),n)=0;

for i=1:n;
XX(:,i)=X(:,2*i-1);
VV(:,i)=X(:,2*i);
end

% Modal Energy Calculation
eta=transpose(inv(U)*transpose(XX)); % eta=XX/U;
etaV=transpose(inv(U)*transpose(VV)); % etaV=VV/U;

E(length(XX),n)=0;
Em(length(eta),n)=0;

for i=2:n;

E(:,i)= 0.5*m(i)*VV(:,i).^2  + 0.5*k(i)*(XX(:,i)-XX(:,1)).^2;
Em(:,i)= 0.5*(W(i,i)*eta(:,i).^2  + etaV(:,i).^2);
end
i=1;
E(:,i)= 0.5*m(i)*VV(:,i).^2  + 0.5*k(i)*XX(:,i).^2;
Em(:,i)= 0.5*(W(i,i)*eta(:,i).^2  + etaV(:,i).^2);

% Total Energy Of System Over Time
Es=E(:,2:n);
Esys=sum(E,2);

% Imparted Energy From Initial Conditions
EIC=sum(E(1,:));

% Imparted Energy From External Excitation
EimpF=cumtrapz(X(:,1),extForce);

clear i VV eta etaV;

% gg_OldDataAdaptToGraphs;


%% Total Energy - Calculation Check 

g01=figure;
plot( T, EIC+EimpF-Esys);

title(['Imparted Energy - Total Energy (Error Check)',simNo])
xlabel('Time (sec)')
ylabel('Energy (j)')

gg_ZoomAndGrid2D;
set(gcf,'Position' , [0,50 ,800,600])
set(gcf, 'Renderer', 'ZBuffer')

saveas(g01,[simname,'-0-CalcCheck.jpg'])
close(g01);
clear g01;


%% Displacement of Primary 

g1=figure;

plot(T,X(:,1))
title(['Displacement of Primary ',simNo])
xlabel('Time (sec)')
ylabel('Displacement (m)')

gg_ZoomAndGrid2D;

set(gcf,'Position' , [0,50 ,800,600])
set(gcf, 'Renderer', 'ZBuffer')

saveas(g1,[simname,'-1-PrimaryDisp.jpg'])
close(g1);
clear g1;
%% FFT of Primary
MAXFREQ=wf;
if sqrt(max(max(W))) > MAXFREQ
    MAXFREQ=sqrt(max(max(W)));
end 

fftPlotRange=[-MAXFREQ/20 1.1*MAXFREQ];     % sets maximum fft frequency to either highest mode 
                                            % or it takes the forcing frequency

[fx,fy]=gg_fft(T,X(:,1),1);

g2=figure;
plot(fx*2*pi,20*log10(fy));
title(['Displacement of Primary in Frequency Domain',simNo]);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (dB)');
grid;
 set(gca,'XLim',fftPlotRange);
 set(gcf,'Position' , [0,50,800,600])
 set(gcf, 'Renderer', 'ZBuffer')

saveas(g2,[simname,'-2-PrimaryFFT.jpg'])

close(g2);
clear g2 fx fy;

%% Energy Of Primary
dispPlotRange=[-simtime/20,simtime];

g00=figure;

plot(T,E(:,1));
title(['Enrgy of Primary',simNo])
xlabel('Time (sec)')
ylabel('Energy (j)')

gg_ZoomAndGrid2D;

saveas(g00,[simname,'-3-PrimaryEn.jpg'])

close(g00);
clear g00;


if (n>2);
%% Plotting Energy Dist Of Satelittes;

g3=figure;

index(length(T),n-1)=0;
for i=1:n-1;index(:,i)=ones(length(T),1)*i+1; 
end

surf(index,T,Es,Es)
shading flat;

% hold;
% plot3(index,T*ones(1,99),Es,'k');

title(['Energy Distribution of Satelittes ',simNo]);
ylabel('Time (sec)');
xlabel('Satellite Index (i)');
zlabel('Energy (J)');

set(gcf,'Position' , [50,50,800,600])
set(gcf, 'Renderer', 'ZBuffer')

view(-60,60);
colormap jet;
colorbar;

saveas(g3,[simname,'-4-SatelittesEn.jpg']);
close(g3);
clear Es g3 index;


%% Plotting Modal Energy Dist

g4=figure;

index(length(T),n)=0;
for i=1:n; index(:,i)=ones(length(T),1)*i; 
end

surf(index,T,Em,Em)
shading flat;

% hold;
% plot3(index,T*ones(1,100),Em,'--k');

title(['Modal Energy Distribution ',simNo]);
ylabel('Time (sec)');
xlabel('Mode Index (i)');
zlabel('Energy (J)');

set(gcf,'Position' , [50,50,800,600])
set(gcf, 'Renderer', 'ZBuffer')

view(-60,60);
colormap jet;
colorbar;

saveas(g4,[simname,'-5-ModsEn.jpg'])
close(g4)
clear g4 Em index;


%% Plotting Satelitte Displacements

index(length(T),n)=0;
for i=1:n-1;
    index(:,i)=ones(length(T),1)*i+1; 
end


g33=figure;

surf(index,T,XX,XX)
shading flat;

title(['Satelitte Movement as a Surface ',simNo]);
ylabel('Time (sec)');
xlabel('Satelitte Index (i)');
zlabel('Displacement (m)');
set(gcf,'Position' , [50,50,800,600])
set(gcf, 'Renderer', 'ZBuffer')

view(60,60);
colormap jet;
colorbar;

saveas(g33,[simname,'-6-SatelittesDisp.jpg']);
close(g33);
clear g33 XX index;

%% Plotting FFT of Satelitte Displacements

for i=1:n-1; 
    [ff,XXF(:,i)]=gg_fft(T,X(:,2*(i+1)-1),1);
    index(:,i)=ones(length(ff),1)*i+1; 
end


g44=figure;

surf(index,2*pi*ff,20*log10(XXF),20*log10(XXF))
shading flat;

alpha(.5);

% plot(2*pi*ff,20*log10(XXF));

hold; 
plot3(index,2*pi*ff,20*log10(XXF),'k');

title(['FFT of Satelitte Movement as a Surface ',simNo]);
ylabel('Frequency (rad/sec)');
xlabel('Satelitte Index (i)');
zlabel('Amplitude (dB)');
set(gca,'YLim',fftPlotRange);
set(gcf,'Position' , [50,50,800,600])
set(gcf, 'Renderer', 'ZBuffer')


view(60,15);
colormap jet;
colorbar;

saveas(g44,[simname,'-7-SatelittesFFT.jpg']);
close(g44)
clear g44 XXF index;


end;

%% Plotting Independent Natural Frequency Dist Of satellites;
g11=figure;
stem(w,'-r*');

WW(n)=0;
for i=1:n;
    WW(i)=W(i,i);
end
hold;
stem(sqrt(WW));
hold;

title(['Independent vs Eigen Frequency Distribution ',simNo]);
xlabel('Satelitte and Mode Index');
ylabel('Frequency (rad/sn)');
set(gcf, 'Renderer', 'ZBuffer') ;
set(gcf,'Position' , [50,50,800,600])
grid

legend('Independent Freqs','Eigen Freqs');
saveas(g11,[simname,'-8-Freqs.jpg']);
close(g11);
clear g11 WW;

%% Input Signal in Time Domain

if (exist('extForce','var')) 

g0=figure;

plot(tForce,extForce);
title(['Input Signal',simNo]);
xlabel('Time (sec)');
ylabel('Amplitude (N)');
grid;
set(gca,'XLim',[0 3*2*pi/wf]);
set(gcf, 'Renderer', 'ZBuffer');

saveas(g0,[simname,'-9-Input.jpg'])
close(g0);
clear g0

%% Input Signal in Frequency Domain 

g01=figure;

[ffx,ffy]=gg_fft(tForce,extForce,1);
plot(ffx*2*pi,20*log10(ffy));
title(['Input Signal in Frequency Domain ',simNo]);
xlabel('Frequency (rad/sn)');
ylabel('Amplitude (dB)');
grid;
set(gca,'XLim',fftPlotRange);
set(gcf, 'Renderer', 'ZBuffer') ;

saveas(g01,[simname,'-10-InputFFT.jpg'])

close(g01);
clear g01 ffx ffy

end
