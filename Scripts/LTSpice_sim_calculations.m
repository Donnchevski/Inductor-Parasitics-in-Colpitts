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
%   This script calculates the relevant values for the steady-state
%   oscillation condition and oscillation frequency using both of the
%   expressions presented in the paper:
%       *) Expression for steady-state oscillation condition for an ideal
%       inductor without parasitics
%       *) Expression for steady-state oscillation condition for a real
%       inductor with parasitic inductor ESR and inductor EPC
% Two simulation examples are considered. The first example emphasizes the
% influence of inductor ESR on steady-state oscillation condition while the
% second example emphasizes the influence of inductor EPC on steady-state
% oscillation frequency.
%
% Relevant variables:
%   VTO - JFET threshold voltage (SPICE parameter)
%   BETA - JFET transconductance (SPICE parameter)
%   LAMBDA - JFET channel-length modulation (SPICE parameter)
%   C1/C2 - LC resonant circuit capacitances
%   L - LC resonant circuit inductance
%   R - Total circuit resistance
%   R_esr - LC inductor equivalent series parasitic resistance
%   C_epc - LC inductor equivalent parallel parasitic capacitance
%   R_L - Load resistance
%   gm - Transistor transconductance
%   w0 - Steady-state oscillation frequency using ideal expression
%   w0z - Steady-state oscillation frequency using ideal equation
%   V_GS - JFET gate to source voltage (DC analysis)
%   V_DS - JFET drain to source voltage (DC analysis)
%   ro - MOSFET small signal output resistance
%   V_DD - DC voltage supply
%
%% 1. Initialization of circuit parameters (influence of R_esr)

V_DD = 5; % V
VTO = -4; % V
BETA = 0.00315; % A/V^2
LAMBDA = 0.014; % V^(-1)
I_D = 1e-3; % A
R_L = 100e3; % Ohms
C1 = 6e-9; % F
C2= C1; % F
R_esr = 2.6; % Ohms
C_epc = 100e-12; % F
L = 3e-6; % H

%% 2. Calculating relevant values (influence of R_esr)

V_GS = VTO + sqrt(I_D/(BETA*(1 + LAMBDA*V_DS)));
V_DS = V_DD + V_GS;
ro = 1/(LAMBDA*I_D);
gm = 2*BETA*(V_GS - VTO)*(1 + LAMBDA*V_DS);
R = ro*R_L/(ro + R_L);
Ce = C1*C2/(C1 + C2);
w0 = sqrt(1/(L*Ce));
w0z = sqrt(1/(L*(Ce + C_epc))*(1 + R_esr/R*(C1 + C_epc)/(C1 + C2) + gm*R_esr*C_epc/(C1 + C2)));
%% 3. Results (influence of R_esr)

fprintf('\nCalculated steady-state oscillation frequency (ideal equation)')
fprintf('\nFrequency: %.0f Hz', w0/(2*pi))
fprintf('\n\nCalculated values for steady-state oscillation condition (ideal equation)')
fprintf('\nLeft side of equation: %.2f', gm*R)
fprintf('\nRight side of equation: %.2f', C1/C2)
fprintf('\n\nCalculated steady-state oscillation frequency (real equation)')
fprintf('\nFrequency: %.0f Hz', w0z/(2*pi))
fprintf('\n\nCalculated values for steady-state oscillation condition (real equation)')
fprintf('\nLeft side of equation: %.2f', gm*R)
fprintf('\nRight side of equation: %.2f\n', C1/C2 + R*R_esr*(1 + C_epc/Ce)*(C1 + C2)/L)

%% 4. Initialization of circuit parameters (influence of C_epc)

V_DD = 5; % V
VTO = -4; % V
BETA = 0.00315; % A/V^2
LAMBDA = 0.014; % V^(-1)
I_D = 1e-3; % A
R_L = 100e3; % Ohms
C1 = 100e-12; % F
C2= C1; % F
R_esr = 2.6; % Ohms
C_epc = 100e-12; % F
L = 3e-6; % H

%% 5. Calculating relevant values (influence of C_epc)

V_GS = VTO + sqrt(I_D/(BETA*(1 + LAMBDA*V_DS)));
V_DS = V_DD + V_GS;
ro = 1/(LAMBDA*I_D);
gm = 2*BETA*(V_GS - VTO)*(1 + LAMBDA*V_DS);
R = ro*R_L/(ro + R_L);
Ce = C1*C2/(C1 + C2);
w0 = sqrt(1/(L*Ce));
w0z = sqrt(1/(L*(Ce + C_epc))*(1 + R_esr/R*(C1 + C_epc)/(C1 + C2) + gm*R_esr*C_epc/(C1 + C2)));

%% 6. Results (influence of C_epc)

fprintf('\nCalculated steady-state oscillation frequency (ideal equation)')
fprintf('\nFrequency: %.0f Hz', w0/(2*pi))
fprintf('\n\nCalculated values for steady-state oscillation condition (ideal equation)')
fprintf('\nLeft side of equation: %.2f', gm*R)
fprintf('\nRight side of equation: %.2f', C1/C2)
fprintf('\n\nCalculated steady-state oscillation frequency (real equation)')
fprintf('\nFrequency: %.0f Hz', w0z/(2*pi))
fprintf('\n\nCalculated values for steady-state oscillation condition (real equation)')
fprintf('\nLeft side of equation: %.2f', gm*R)
fprintf('\nRight side of equation: %.2f\n', C1/C2 + R*R_esr*(1 + C_epc/Ce)*(C1 + C2)/L)