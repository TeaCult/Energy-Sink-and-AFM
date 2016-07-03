function dx=gg_EnergySinkSolver(t,x,M,K,C,simtime,waitbarptr)
% Designed to solve MDOF system which consist of parallel oscillators
% attached to a primary mass and spring.
% Written by Gediz GURSU 8/12/2010.
% t - time
% x - position
% M - Mass matrix
% K - Stiffness Matrix
% C - Damping Matrix
% Forcing will be added to file by fopen fwrite instead of interp1

load('temp.mat');

Meq=mc*0.178;
Keq=1.07*kc;

%% Solution For MDOF System
n=2*rank(M);
dx=zeros(n,1);

A(n/2,1)=0;
V(n/2,1)=0;
X(n/2,1)=0;
F(n/2,1)=0;

for i=2:2:n;
    j=i/2;
    
    A(j,1)=dx(i);  % Acceleration Vector
    V(j,1)=x(i);   % Velocity Vector
    X(j,1)=x(i-1); % Position Vector
    F(j,1)=0;      % Force Applied to sattalites 
    
 dx(i-1)=V(j,1);   % Assign Position Vector Back for solver iterration
end;
F(1,1)=sin(MODIFIER*t);   % Force Applied to Primary 

A=inv(M)*(-K*X-C*V+F); % Calculation for itteration
                                                        
% Assign Acceleration Vector Back for solver iterration

for i=2:2:n;
    j=i/2;
    dx(i)=A(j,1); 
end;

% Progress indicator updater 
if( (simtime/t-round(simtime/t)) < 0.01); 
waitbar(t/simtime,waitbarptr); 
end;
