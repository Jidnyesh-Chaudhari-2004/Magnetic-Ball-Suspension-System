% Clear workspace and command window
clear; clc;

% Define symbolic variables and unknowns
syms s h2 h1 h0 k1

% Given parameters
zeta = 0.69;
w_n  = 3.0627;
kx   = 2180;
b    = 3694.61;

% Define alpha
alpha = 500 * (zeta * w_n);


% Define the polynomial delta(s)
del_s = expand((s^2 + 2*zeta*w_n*s + w_n^2) * (s + alpha)^2);

% Define gc_1(s) corresponding to:
% s^4 + k1*s^3 - s^2*(kx + b*h2) - s*(k1*kx + b*h1) - b*h0
gc1_s = expand(s^4 + k1*s^3 - s^2*(kx + b*h2) - s*(k1*kx + b*h1) - b*h0);

% Convert delta(s) to a coefficient vector (highest degree first)
p_del = sym2poly(del_s);
% p_del should be:
% [1, 2*(alpha + zeta*w_n), alpha^2 + 4*zeta*w_n*alpha + w_n^2, 2*(zeta*w_n*alpha^2 + w_n^2*alpha), w_n^2*alpha^2]

% In gc1_s, the coefficients are:
% s^4:    1
% s^3:    k1
% s^2:   -(kx + b*h2)
% s^1:   -(k1*kx + b*h1)
% s^0:   -b*h0

% Equate coefficients for corresponding powers of s
eq1 = k1 == p_del(2);            % Coefficient of s^3
eq2 = -(kx + b*h2) == p_del(3);    % Coefficient of s^2
eq3 = -(k1*kx + b*h1) == p_del(4); % Coefficient of s^1
eq4 = -b*h0 == p_del(5);           % Constant term

% Solve for k1, h2, h1, and h0
sol = solve([eq1, eq2, eq3, eq4], [k1, h2, h1, h0]);

% Simplify the solutions
k1_sol = simplify(sol.k1);
h2_sol = simplify(sol.h2);
h1_sol = simplify(sol.h1);
h0_sol = simplify(sol.h0);

% Convert symbolic solutions to numeric values
k1_val = double(vpa(k1_sol, 6));
h2_val = double(vpa(h2_sol, 6));
h1_val = double(vpa(h1_sol, 6));
h0_val = double(vpa(h0_sol, 6));

% Display the computed values
disp('Computed parameter values:');
fprintf('k1 = %f\n', k1_val);
fprintf('h2 = %f\n', h2_val);
fprintf('h1 = %f\n', h1_val);
fprintf('h0 = %f\n', h0_val);

% Define numerator and denominator coefficients of the transfer function
% Numerator: [-b*h2, -b*h1, -b*h0]
% Denominator: [1, k1, -kx, -k1*kx, 0]
num = [-b*h2_val, -b*h1_val, -b*h0_val];
den = [1, k1_val, -kx, -k1_val*kx, 0];

% Create the transfer function model
sys = tf(num, den);

% Plot the Bode plot with gain and phase margins
figure;
margin(sys);
grid on;
title('Bode Plot of the Transfer Function');
