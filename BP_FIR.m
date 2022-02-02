%%
clc; clear all; close all;
%% BP Filter Specs
omega_s1 = 0.314 * pi;
omega_s2 = 0.4836 * pi;
omega_p1 = 0.3382 * pi;
omega_p2 = 0.4594 * pi;
transition_bw = 0.0242 * pi;
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
%%
width = width_min + 16;
w = kaiser(width,beta);
ideal_bpf = ideal_lpf(omega_c2, width) - ideal_lpf(omega_c1, width);
fir_bpf = ideal_bpf .* w';
[H,f] = freqz(fir_bpf,1,1024, 330e3);
figure();
plot(f, abs(H));
hold on;
xline(55.8e3, 'magenta--', 'LineWidth', 1.5);
hold on;
xline(75.8e3, 'magenta--', 'LineWidth', 1.5);
hold on;
yline(1.15, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.85, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.15, 'red--', 'LineWidth', 1.5); 
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j2 \pi f})');
title('Magnitude Response of Discrete Time FIR Bandpass Filter');
set(gca, 'XTick', [51.8e3, 55.8e3, 75.8e3, 79.8e3], 'xticklabel', {'f_{s1}', 'f_{p1}', 'f_{p2}', 'f_{s2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});

%%
disp(fir_bpf);
%%
fvtool(fir_bpf, 'Analysis', 'Phase');
fvtool(fir_bpf, 'Analysis', 'Impulse');