% *************************************************************************
% LICENSE
% *************************************************************************
% Code by Filip Donchevski
% For Paper, "Effects of Inductor Parasitics on Loop Gain in Single-Ended
% Common-Drain Colpitts Oscillator"
% by F. Donchevski, Z. Kokolanski, and M. Stankovski
%
% 1. Grant of License
% 
% You are free to:
% 
%     Use this software for personal, educational, and research purposes.
%     Modify and distribute the software, provided you retain this license
%     notice. Share the software with others under these same terms.
% 
% 2. Commercial Use Restriction
% 
% This software cannot be used for commercial purposes without prior
% written permission from the original author(s). Commercial use includes,
% but is not limited to:
% 
%     Selling, leasing, or licensing the software. Using the software in a
%     product or service that generates revenue. Integrating the software
%     into proprietary software or commercial applications.
% 
% 3. Commercial Use Licensing
% 
% If you wish to use this software for commercial purposes, you must obtain
% a separate commercial license. Please contact Filip Donchevski
% (fdoncevski@feit.ukim.edu.mk) for licensing inquiries.
% 
% 4. Disclaimer
% 
% THIS SOFTWARE IS PROVIDED "AS IS," WITHOUT WARRANTIES OF ANY KIND. THE
% AUTHORS ARE NOT LIABLE FOR ANY DAMAGES ARISING FROM THE USE OF THIS
% SOFTWARE.
%
%
%
% *************************************************************************
% CODE INFORMATION
% *************************************************************************
% Relevant information:
%   This script generates the root locus plots for an ideal and non-ideal
%   single-ended common-drain Colpitts oscillator. The non-ideal oscillator
%   considers non-ideal inductor with equivalent series resistance and
%   equivalent parallel capacitance. The script is organized as following:
%
%       1. Initialization of circuit parameters
%       2. Calculating oscillation condition for an ideal oscillator
%       3. Root locus plot of an ideal oscillator
%       4. Root locus plot of an ideal oscillator (magnified)
%       5. Finding imaginary poles location for an ideal oscillator
%       6. Calculating oscillation condition for a non-ideal oscillator
%       7. Root locus plot of a non-ideal oscillator
%       8. Root locus plot of a non-ideal oscillator (magnified)
%       9. Finding imaginary poles location for a non-ideal oscillator
%
% Relevant variables:
%   L - LC resonant circuit inductance
%   C1/C2 - LC resonant circuit capacitances
%   R_esr - LC inductor equivalent series parasitic resistance
%   ro - MOSFET small signal output resistance
%   R_L - Load resistance
%   C_epc - LC inductor equivalent parallel parasitic capacitance
%   R - Total circuit resistance
%   w0 - steady state oscillation frequency
%   cond -Expression for steady-state oscillation condition
%   Ce - LC resonant circuit equivalent capacitance
%   gmR - product of transistor transconductance and total circuit resist.
%
%% 1. Initialization of circuit parameters

syms s
C1 = 100e-12; % F
C2 = 100e-12; % F
L = 3e-6; % H
R_L = 100e3; % Ohm
ro = 50e3; % Ohm
C_epc = 0; % F
R_esr = 0; % Ohm
R = ro*R_L/(ro + R_L); % Ohm
Ce = C1*C2/(C1 + C2);
fs = 11; % font size

%% 2. Root locus plot of an ideal oscillator

gmR = C1/C2

%% 3. Root locus plot of an ideal oscillator
Gos = (L*C_epc*s^2 + R_esr*C_epc*s + 1)/(L*R*(C1*C2 + C1*C_epc + C2*C_epc)*s^3 + (L*(C1 + C_epc) + R*R_esr*(C1*C_epc + C2*C_epc + C1*C2))*s^2 + (R*(C1 + C2) + R_esr*(C1 + C_epc))*s + 1);
[N, D] = numden(Gos);
Af = tf(sym2poly(N), sym2poly(D));

figure(1)
rlocus(Af)
xlim([-3e8 3e8])
ylim([-3e8 3e8])
title("Root Locus of $A_f$ (without parasitics)", 'interpreter','latex');
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
set(axIm,'String','Imaginary Axis', 'interpreter', 'latex');
set(axRe,'String','Real Axis', 'interpreter', 'latex');
text(-0.8e8,0.75e8,'$g_m = 0$', 'interpreter','latex', 'FontSize', fs);
text(-0.8e8,-0.75e8,'$g_m = 0$', 'interpreter','latex', 'FontSize', fs);
text(0.05e8, -0.2e8,'$g_m = 0$', 'interpreter','latex', 'FontSize', fs);
text(-2.9e8, -0.2e8,'$g_m \to \infty$', 'interpreter','latex', 'FontSize', fs);
text(0.5e8, 2.6e8,'$g_m \to \infty$', 'interpreter','latex', 'FontSize', fs);
text(0.5e8, -2.6e8,'$g_m \to \infty$', 'interpreter','latex', 'FontSize', fs);


%% 4. Root locus plot of an ideal oscillator (magnified)
Gos = (L*C_epc*s^2 + R_esr*C_epc*s + 1)/(L*R*(C1*C2 + C1*C_epc + C2*C_epc)*s^3 + (L*(C1 + C_epc) + R*R_esr*(C1*C_epc + C2*C_epc + C1*C2))*s^2 + (R*(C1 + C2) + R_esr*(C1 + C_epc))*s + 1);
[N, D] = numden(Gos);
Af = tf(sym2poly(N), sym2poly(D));

figure(2)

rlocus(Af, linspace(0, 1.1, 1e6))

xlim([-1.55e5 0.06e5])
ylim([-3e8 3e8])
title("Root Locus of $A_f$ (without parasitics)", 'interpreter','latex');
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
set(axIm,'String','Imaginary Axis', 'interpreter', 'latex');
set(axRe,'String','Real Axis', 'interpreter', 'latex');

%% 5. Finding imaginary poles(gain) location for an ideal oscillator
[poles, gains] = rlocus(Af, linspace(0, 1.1, 1e6));
impole_ind = knnsearch(real(poles(1, :))', 0);
impole1 = poles(1, impole_ind)
impole2 = poles(2, impole_ind)
gain = gains(impole_ind)

%% Defining finite values for inductor parasitics

C_epc = 5e-12; % F
R_esr = 2.6; % Ohm

%% 6. Calculating oscillation condition for a non-ideal oscillator

gmR = (C1/C2 + R*R_esr*(1 + C_epc/Ce)*((C1 + C2)/L))


%% 7. Root locus plot of a non-ideal oscillator

Gos = (L*C_epc*s^2 + R_esr*C_epc*s + 1)/(L*R*(C1*C2 + C1*C_epc + C2*C_epc)*s^3 + (L*(C1 + C_epc) + R*R_esr*(C1*C_epc + C2*C_epc + C1*C2))*s^2 + (R*(C1 + C2) + R_esr*(C1 + C_epc))*s + 1);
[N, D] = numden(Gos);
Af = tf(sym2poly(N), sym2poly(D));

figure(3)
rlocus(Af)

xlim([-3e8 3e8])
ylim([-3e8 3e8])
title("Root Locus of $A_f$ (with parasitics)", 'interpreter','latex');
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
set(axIm,'String','Imaginary Axis', 'interpreter', 'latex');
set(axRe,'String','Real Axis', 'interpreter', 'latex');
text(-0.8e8,0.75e8,'$g_m = 0$', 'interpreter','latex', 'FontSize', fs);
text(-0.8e8,-0.75e8,'$g_m = 0$', 'interpreter','latex', 'FontSize', fs);
text(0.05e8, -0.2e8,'$g_m = 0$', 'interpreter','latex', 'FontSize', fs);
text(-2.9e8, -0.2e8,'$g_m \to \infty$', 'interpreter','latex', 'FontSize', fs);
text(-0.9e8, 2.6e8,'$g_m \to \infty$', 'interpreter','latex', 'FontSize', fs);
text(-0.9e8, -2.6e8,'$g_m \to \infty$', 'interpreter','latex', 'FontSize', fs);
% grid on;

%% 8. Root locus plot of a non-ideal oscillator (magnified)

Gos = (L*C_epc*s^2 + R_esr*C_epc*s + 1)/(L*R*(C1*C2 + C1*C_epc + C2*C_epc)*s^3 + (L*(C1 + C_epc) + R*R_esr*(C1*C_epc + C2*C_epc + C1*C2))*s^2 + (R*(C1 + C2) + R_esr*(C1 + C_epc))*s + 1);
[N, D] = numden(Gos);
Af = tf(sym2poly(N), sym2poly(D));

figure(4)
[poles, gains] = rlocus(Af);
rlocus(Af, [linspace(0, 7.3, 1e2) linspace(7.3, 7.4, 1e6) linspace(7.4, 7.7, 1e2) gains(40:55)])


xlim([-5.2e5 0.23e5])
ylim([-3e8 3e8])
title("Root Locus of $T$ (with parasitics)", 'interpreter','latex');
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
set(axIm,'String','Imaginary Axis', 'interpreter', 'latex');
set(axRe,'String','Real Axis', 'interpreter', 'latex');
%% 9. Finding imaginary poles location for a non-ideal oscillator

[poles, gains] = rlocus(Af, linspace(7.3, 7.4, 1e7));
impole_ind = knnsearch(real(poles(1, :))', 0);
impole1 = poles(1, impole_ind)
impole2 = poles(2, impole_ind)
gain = gains(impole_ind)
