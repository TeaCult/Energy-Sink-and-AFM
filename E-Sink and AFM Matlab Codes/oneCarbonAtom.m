function dx = oneCarbonAtom(t,x)

% ONE CARBON ATOM

avagadro=6.022142*10^23; % avagadro number
Ebond=280/avagadro; % Bond Energy per kj/mol/avagadro
kd=2.55*10^10; %decay factor
ze = 2.42*10^-10;  % bond distance 2.42 A in m
mc = 12/avagadro;

%Vpott=-Ebond*(-2*exp(-kd*(z-ze))+exp(-2*kd*(z-ze))); %therefore
%Fts=Ebond*((2*kd)/exp(kd*(z - ze)) - (2*kd)/exp(2*kd*(z - ze)));

% Solution of Equations of Motion
dx = zeros(2,1);    % a column vector
dx(1)=x(2);
dx(2)= -Ebond*((2*kd)/exp(kd*(x(1) - ze)) - (2*kd)/exp(2*kd*(x(1) - ze)))/mc;