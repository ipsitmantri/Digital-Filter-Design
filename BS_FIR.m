%%
clc; clear all; close all;
%% BS Filter Specs
omega_s1 = 0.3769 * pi;
omega_s2 = 0.5307 * pi;
omega_p1 = 0.3461 * pi;
omega_p2 = 0.5615 * pi;
transition_bw = 0.0308 * pi;
delta = 0.15;
omega_c1 = 0.5 * (omega_p1 + omega_s1);
omega_c2 = 0.5 * (omega_p2 + omega_s2);
%% Kaiser window parameters
A = -20 * log10(delta);
fprintf('A = %f\n', A);

width_min = ceil( 1 + ((A - 8)/(2.285 * transition_bw)));
fprintf('Minimum width = %d\n', width_min);

alpha = -1;
if A < 21
    alpha = 0;
elseif A >=21 && A <= 50
    alpha = 0.5842 * (A - 21) ^ 0.4 + 0.07886 * (A - 21);
elseif A > 50
    alpha = 0.1102 * (A - 8.7);
else
    
end
beta = alpha /width_min;
fprintf('beta = %f\n', beta);
%% Magnitude Response
width = width_min + 13;
w = kaiser(width,beta);
ideal_bsf = ideal_lpf(pi, width) - ideal_lpf(omega_c2, width) + ideal_lpf(omega_c1, width);
fir_bsf = ideal_bsf .* w';
[H,f] = freqz(fir_bsf,1,1024, 260e3);
figure();
plot(f, abs(H));
hold on;
xline(49e3, 'magenta--', 'LineWidth', 1.5);
hold on;
xline(69e3, 'magenta--', 'LineWidth', 1.5);
hold on;
yline(1.15, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.85, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.15, 'red--', 'LineWidth', 1.5); 
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j2 \pi f})');
title('Magnitude Response of Discrete Time FIR Bandstop Filter');
set(gca, 'XTick', [45e3, 49e3, 69e3, 73e3], 'xticklabel', {'f_{p1}', 'f_{s1}', 'f_{s2}', 'f_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
%%
disp(fir_bsf);
%%
fvtool(fir_bsf, 'Analysis', 'Phase');
fvtool(fir_bsf, 'Analysis', 'Impulse');