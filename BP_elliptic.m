% Please execute this code section wise
% so that the values of the variables
% settle and are stable
%

clc; clear ; close all;
%% Elliptic Lowpass filter transfer function and plots
[z, p, H0, B, A] = ellipap2(3, 1.411, 16.48);
syms sl Omega_L;
H_LPF = (1 + 0.6294 * sl^2) / ((1 + 1.6046 * sl) * (1 + 0.2306 * sl + 0.9994 * sl^2));
H_LPF_freq = subs(H_LPF, sl, 1i * Omega_L);
[ns, ds] = numden(H_LPF);
nsl = sym2poly(ns);
dsl = sym2poly(ds);
k = ds(1);
nsl = nsl / k;
ds = ds / k;
disp(nsl);
disp(dsl);
figure();
fplot(abs(H_LPF_freq));
hold on;
fplot(sl - sl - 0.15 + 1, 'k-', 'Markersize', 10);
hold on;
fplot(sl-sl+1, 'k-', 'Markersize', 10);
hold on;
fplot(sl - sl + 0.15, 'r-', 'Markersize', 10);
xline(1, 'magenta-');
hold on;
xline(1.3874,'magenta-');
axis([0 2 0 1.2]);
set(gca, 'XTick', [0, 1, 1.3874], 'xticklabel', {'0', '\Omega_{Lp} = 1','\Omega_{Ls} = 1.3874'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
daspect([1 1 1]);
title('Equivalent Elliptic Lowpass filter mangitude response');
xlabel('\Omega_L');
ylabel('|H(j \Omega_L)|');

%% Analog lowpass to bandpass frequency transformation
syms s Omega;
Omega_0 = 0.7189;
B = 0.2924;
H_BPF = subs(H_LPF, sl, (s^2 + Omega_0^2) / (B * s));
[ns, ds] = numden(H_BPF);
ns = sym2poly(ns);
ds = sym2poly(ds);
k = ds(1);
ns = ns / k;
ds = ds / k;
disp(ns);
disp(ds);
H_BPF_freq = subs(H_BPF, s, 1i * Omega);
figure();
fplot(abs(H_BPF_freq));
set(gca, 'XTick', [0.5375, 0.5875, 0.8799, 0.9498], 'xticklabel', {'\Omega_{s1}', '\Omega_{p1}', '\Omega_{p2}', '\Omega_{s2}'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
hold on;
xline(0.5375, 'magenta-');hold on;
xline(0.5875, 'magenta-');hold on;
xline(0.8799, 'magenta-');hold on;
xline( 0.9498, 'magenta-');
hold on;
fplot(s - s - 0.15 + 1, 'k-', 'Markersize', 10);
hold on;
fplot(s-s+1, 'k-', 'Markersize', 10);
hold on;
fplot(s - s + 0.15, 'r-', 'Markersize', 10);
daspect([1 1 1]);
title('Elliptic Bandpass filter mangitude response');
xlabel('\Omega');
ylabel('|H(j \Omega)|');
axis([0 2 0 1.2]);
%% Analog to z bilinear transformation
syms z;
Hz = subs(H_BPF, s, (z-1)/(z+1));
[Nz, Dz] = numden(Hz);
Nz = sym2poly(Nz);
Dz = sym2poly(Dz);
k = Dz(1);
Nz = Nz / k;
Dz = Dz / k;
disp(Nz);
disp(Dz);
[H,f] = freqz(Nz,Dz,1024*1024, 330e3);
plot(f,abs(H));
axis([0 100e3 0 1.3]);
hold on;
yline(1.00, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.85, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.15, 'red--', 'LineWidth', 1.5);
hold on;
xline(55.8e3, 'magenta--', 'LineWidth', 1.5);
hold on;
xline(75.8e3, 'magenta--', 'LineWidth', 1.5);
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j 2\pi f})|');
title('Magnitude Response of the Discrete Time Bandpass Filter');
set(gca, 'XTick', [51.8e3, 55.8e3, 75.8e3, 79.8e3], 'xticklabel', {'f_{s1}', 'f_{p1}', 'f_{p2}', 'f_{s2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
fvtool(Nz, Dz, 'Analysis', 'Phase');