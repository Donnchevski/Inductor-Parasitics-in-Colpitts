% *************************************************************************
% LICENSE
% *************************************************************************
% Code by Filip Donchevski For Paper, "Effects of Inductor Parasitics on
% Loop Gain in Single-Ended Common-Drain Colpitts Oscillator" by F.
% Donchevski, Z. Kokolanski, and M. Stankovski
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
%   This script derives the steady-state oscillation frequency and
%   oscillation condition for a single-ended common-drain Colpitts
%   oscillator with inductor equivalent series resistance and equivalent
%   parallel capacitance. In order to ensure consistency and error proof
%   derivations, the expressions were derived using three different
%   approaches:
%       *) Finding characteristic equation (CE) through KCL
%       *) Finding CE through loop gain approach
%       *) Finding CE by equating model determinant to zero
%
% Relevant variables:
%   L - LC resonant circuit inductance
%   C1/C2 - LC resonant circuit capacitances
%   Rs - LC inductor equivalent series parasitic resistance
%   Cp - LC inductor equivalent parallel parasitic capacitance
%   R - Total circuit resistance
%   gm - Transistor transconductance
%   w - frequency
%   CE - Expression for characteristic equation
%   CEw - CE where s is substituted with jw
%   re - Real part of CEw
%   im - Imaginary part of CEw
%   w0 - Expression for steady-state oscillation frequency
%   w0approx - Approximated form of w0
%   cond -Expression for steady-state oscillation condition
%   condapprox -Approximated form of cond
%
%% Defining symbolic variables and relevant expressions

syms L C1 R C2 gm w w0 Cp Rs positive
syms s

Xl = L*s + Rs;
Xc = 1/(s*Cp);
X3 = simplify(Xl*Xc/(Xl + Xc));
X1 = 1/(s*C1);
X2 = 1/(s*C2);

%% Derivations of w0approx and condapprox with KCL approach

s1 = simplify(1/(X3 + 1/(s*C1)), 'Steps', 100);
s2 = 1/(X3);
vo = 1/s1*s2;
s3 = (s*C2*R + 1)/R;
s4 = collect(vo*(s1 + s3 + gm) - gm);
[n, d] = numden(s4);
CE = collect(simplify(n, "Steps", 1000))

% substituting jw in place of s
CEw = subs(n, s, 1j*w);
re = real(CEw);
im = imag(CEw);

w0 = simplify(solve(im == 0, w));
w0approx = sqrt(1/((C1*C2/(C1 + C2) + Cp)*L))
cond = simplify(subs(re, w, w0) == 0, "Steps", 1000);
condapprox = simplify(subs(re, w, w0approx) == 0, "Steps", 1000)

%% Derivations of w0approx and condapprox with loop gain approach

Z_L = simplify(X2*(X1 + X3)/(X1 + X2 + X3), "Steps", 1000);
beta = simplify((X3)/(X1 + X3));
Z = Z_L*R/(Z_L + R);
Av = simplify(gm*Z/(1 + gm*Z));
T = simplify(beta*Av, "Steps", 1000);
[n, d] = numden(T);
CE = collect(simplify(d-n, "Steps", 1000))

% substituting jw in place of s
CEw = subs(CE, s, 1j*w);
re = real(CEw);
im = imag(CEw);

w0 = simplify(solve(im == 0, w));
w0approx = sqrt(1/((C1*C2/(C1 + C2) + Cp)*L))
cond = simplify(subs(re, w, w0) == 0, "Steps", 1000);
condapprox = simplify(subs(re, w, w0approx) == 0, "Steps", 1000)

%% Derivations of w0approx and condapprox with determinant approach

A = [R + 1/(s*C2) + gm*R*1/(s*C2), -(gm*R*X3 + 1/(s*C2)*(gm*R + 1));
     -1/(s*C2), X3 + + 1/(s*C1) + 1/(s*C2)];
CE = simplify(collect(det(A)), "Steps", 1000)
[n, d] = numden(CE);

% substituting jw in place of s
CEw = subs(n, s, 1j*w);
re = real(CEw);
im = imag(CEw);
w0 = simplify(solve(im == 0, w));
w0approx = sqrt(1/((C1*C2/(C1 + C2) + Cp)*L))
cond = simplify(subs(re, w, w0) == 0, "Steps", 1000);
condapprox = simplify(subs(re, w, w0approx) == 0, "Steps", 1000)
