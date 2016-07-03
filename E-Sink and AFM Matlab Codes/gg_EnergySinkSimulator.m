% This script is designed to run mdof simulations for energy sink system.
% Important details are included and commented in the script file.
% Written by Gediz GURSU 8/12/2010

clc;
clear all;
close all;

load('temp.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Variables
% 0.5 1 1.2 1.4 2 5 
simname=['ES-LFD-HE_Batch',num2str(pubCount)];       % Simulation Name
simtime=6000;                % Simulation Time

n=99; % Number of sattalites
                   
MB=1; 
KB=1;
CB=0;
msat=MB/10/n;                                      

% Linear Frequency Distribution :
wtemp=linspace(2/n,2,n); 

% Optimum Frequency Distribution :
%wtemp=[1.25046342663753,  1.21793947218943,  1.00287175105722, 0.829503204019052,  0.994275840616342,  0.990037114421945,  0.878955118528883,  0.976658557288460, 0.883337610340081,  1.08333624207787,  0.874491508119232,  1.17872996443854,     0.887746808658911,  1.15180402512120,  1.07296422457843,  0.751537039383106,      0.790483820637464,  0.905795442762578,  1.05802431446775,  0.914842187442993,      0.985842526044666,  1.23903629674719,  0.952109652092449,  0.843099352145036,      1.02954081273606,  0.966679060643153,  0.847590209740284,  0.856578284420052,      0.869989274736065,  0.852110555888363,  1.20795613516513,  0.981874031033416,      1.11046271067682,  0.942585186089950,  1.13954745891306,  0.924035952578139,      0.947320681238165,  0.971617582187494,  0.892150650378390,  0.919421559110920,      1.10488796347336,  0.937913494599069,  1.17190083920356,  0.928655976535449,      1.01601014968708,  0.933260740696171,  1.12768058061995,  0.901311267952485,      1.04827523124165,  0.824855274584758,  1.11605182499893,  0.834029193653786,      1.03414623167133,  0.961781657585011,  1.32834080099387,  1.20076116539589,     1.29060001098521,  1.19359174161462,  0.896774961110915,  1.18573665036900,      1.30754500481979,  1.12177337509548,  1.05313985487859,  0.838572997746105,      1.09938036855499,  0.800826411766995,  1.07809827521806,  0.820165546813477,     1.03880637115979,  0.861015501126513,  1.14561018590893,  0.956925830720000,      1.04346214472772,  0.773672862166070,  1.01155225778866,  0.810676785734697,      0.910274963443276,  1.16516475720235,  1.02497996083436,  1.13361881393468,     0.865466552066545,  1.22826076345449,  0.785130028902354,  0.795712997161774,      0.767286440816674,  0.760073282337135,  0.998525806992197,  1.15797528731106,      0.779591037255434,  1.00718180579661,  1.27597187828249,  0.805836277584847,      0.815471862864116,  1.02043518068529,  1.26282207192345,  1.08861745881931,      1.06791401873149,  1.09398539432362,  1.06294125793547];

w=[1,sort(wtemp)];

% Initial Conditions
n=n+1;  % changing number of sattalite variables to mass index var;

IC(2*n)=0;
IC(2)=0;

% Force Applied to BASE as Fp=sin(t*wf+phi)


wf=MODIFIER; % rad/sec
% phi=0;
% dt=1/wf/20;
% srate=simtime/dt;
ForceSym=['sin(',num2str(MODIFIER),'*t)'];
% tForce=linspace(0,simtime,srate);
% extForce=subs(ForceSym,'t',tForce);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simNo='';
simDescription=0;
fmaxOfWN=2;                    % max frequency for FFT
fminOfWN=0/(2*pi);             % min frequency for FFT

n=n-1; % changing back to sattalite index

%% Creation of Mass And Stiffness Matrices
m(n+1)=0;
c(n+1)=0;
k(n+1)=0;

m(:)=msat;
c(:)=0;

for i=1:n+1;    
    k(i)=w(i)^2*m(i);
end;

m(1)=MB;
c(1)=CB;
k(1)=KB;

M=zeros(length(k));
K=M;
C=M;

for i=1:length(m);
M(i,i)=m(i);
K(i,i)=k(i);
C(i,i)=c(i);
end;

K(1,:)=-k(:);
K(:,1)=-k(:);
K(1,1)=sum(k);

C(1,:)=-c(:);
C(:,1)=-c(:);
C(1,1)=sum(c);
n=n+1; % changing to normal index again; 

%% Solution
h=waitbar(0,'ODE Solution Progress');
options = odeset('RelTol',1e-4,'AbsTol',1e-6);%,'InitialStep',1e-24);%,'Refine',1,'InitialStep',1e-24);

[T,X] = ode23t(@(t,x) gg_EnergySinkSolver(t,x,M,K,C,simtime,h),[0 simtime],IC,options); % Checked if stiff and ode113 is found to be the best precision and best performance
delete(h);

save(simname); % exporting data using simulation name