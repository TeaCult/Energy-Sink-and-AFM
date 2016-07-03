clear all;
close all;
clc

% kc=42; % cantilever stiffness
% wc = 300*10^3*2*pi; % cantilever natural frequency
% mc = kc/wc^2 ; % cantilver mass
% LinPot=2*D*b^2; % Potential is Linearized from Akays Paper

Angst=1e-10;                % Angstrom;
Nano=1e-9;                  % Nano Scale Factor
PM=Angst/100;               % PicoMeter
Avagadro=6.022142*10^23;    % Avagadro number
AH=1e-18;                     % Hamaker constant in eV;
kN2N=1e3;                   % Kilo Newton to Newton


b = 2.5e10;         % decay factor
a = 74*PM;          % bond distance 2.42 A in m
D = 436;            % Bond Energy per kj/mol

z=linspace(0,20*Angst,1000);
z=z';

VMorse=-D*(2*exp(-b*(z - a))-exp(-2*b*(z - a))); % kJ/mol
FMorse=(2*b*exp(b*(a - z)) - 2*b*exp(2*b*(a - z)))*D; % kNewton /mol

FMorsePerMolecule=FMorse/Avagadro*kN2N; % Newton per molecule


% Decay Factor Confirmation From a Web Source using H-H Bond
% http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html

plot(z/PM,VMorse);
set(gca,'XLim',[30 250]);
set(gca,'YLim',[-490 100]);
grid;

title('H-H Van Der Waals Forces');
xlabel('Distance z(pm)')
ylabel('Potential Energy (kj/Mol)');



% AFM Moleculer Force Confirmation : Bushan Springers HandB. of Nano.
clear VMorse FMorse FMorsePerMolecule;

a = 296*PM;          % bond distance 2.42 A in m
D = 130;            % Bond Energy per kj/mol
b = 2.5e10;         % decay factor

VvdW=-AH./(6*z);
FvdW=AH./(6.*z.^2);

VMorse=-D*(2*exp(-b*(z - a))-exp(-2*b*(z - a))); % kJ/mol
FMorse=(2*b*exp(b*(a - z)) - 2*b*exp(2*b*(a - z)))*D; % kNewton /mol
FMorsePerMolecule=FMorse/Avagadro*kN2N; % Newton per molecule

FTotal=FMorsePerMolecule/Nano+FvdW;

figure;
hold;

plot(z/Angst,-FMorsePerMolecule/Nano,'-b');
plot(z/Angst,-FvdW,'-g');
plot(z/Angst,-FTotal,'-r');

set(gca,'XLim',[0 20]);
set(gca,'YLim',[-5 5]);
grid;
title('Molecular Forces of a typical Tip Sample Interraction Fts');
xlabel('Distance z(Angst)')
ylabel('Fts (nN)');

legend('Short-Range Morse Potential', 'Long-Range Van der Waals Force','Total Force');


% 
% 
% plot(x,Vpott)
% figure;
% plot(x,Fts)
