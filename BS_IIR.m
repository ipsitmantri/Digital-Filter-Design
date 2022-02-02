clc; clear all; close all;
% transformed normalized filter specs
syms x B Omega_0 f(x);
B = 0.6102;
Omega_0 = 0.8568;
f(x) = (B * x) / (-x^2 + Omega_0^2);
disp(vpa(f(0.6725), 5));
disp(vpa(f(1.1014), 5));
disp(vpa(f(0.6044), 5));
disp(vpa(f(1.2146), 5));
%% equilavent lowpass filter specs
D2 = 43.44;
D1 = 0.384;
Omega_ls = 1.4031;
Omega_lp = 1;
epsilon = sqrt(D1);
N = ceil((acosh(sqrt(D2/D1))) / (acosh(Omega_ls / Omega_lp)));
fprintf('N = %d\n', N);
%% finding out the left half poles
k = 0:2 * N - 1;
Ak = (2 * k + 1) .* (pi /(2 * N));
Bk = (1 / N) * asinh(1 / epsilon);
poles = (-1 * sin(Ak) * sinh(Bk)) + 1i * (cos(Ak) * cosh(Bk));
figure();
plot(poles, '.', 'MarkerSize', 20);
daspect([1 1 1]);
hold on;
x = linspace(-pi, pi, 10000);
a = sinh(Bk) .* cos(x);
b = cosh(Bk) .* sin(x);
plot(a,b);
hold on;
plot(0,0,'r*');
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
left_plane_poles = [];
for i=1:size(poles,2)
    hold on;
    plot([0, real(poles(1,i))], [0, imag(poles(1,i))], 'r-');
    if real(poles(1,i)) < 0
        disp(poles(1,i));
        left_plane_poles = [left_plane_poles, poles(1,i)];
    end
end
plot([0, 0],[-1.5, 1.5], 'k-');
plot([-1.5, 1.5],[0, 0], 'k-');
title('Poles of $H_{LPF}(s_L) \, H_{LPF}(-s_L)$ in the $s_L$ plane', 'Interpreter', 'latex');
xlabel('\Sigma_k');
ylabel('\Omega_k');
%% Analog Lowpass transfer function in s_L domain
num = 1;
den = 1;
syms sl;
for i=1:size(left_plane_poles, 2)
    num = num * left_plane_poles(1,i);
    den = den * (sl - left_plane_poles(1,i));
end
H_LPF = num / (den * sqrt(1 + D1));
disp(expand(vpa(den * sqrt(1 + D1), 5)));
disp(num);
syms Omega_L;
H_LPF_freq = subs(H_LPF, sl, 1i * Omega_L);
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
xline(1.4031,'magenta-');
axis([0 2 0 1.2]);
set(gca, 'XTick', [0, 1, 1.4031], 'xticklabel', {'0', '\Omega_{Lp} = 1','\Omega_{Ls} = 1.4031'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
daspect([1 1 1]);
title('Equivalent Chebyshev Lowpass filter mangitude response');
xlabel('\Omega_L');
ylabel('|H(j \Omega_L)|');

figure();
fplot(angle(H_LPF_freq));
xlabel('\Omega_L');
ylabel('\angle H(j \Omega_L)');
title('Equivalent Chebyshev Lowpass filter phase response');
%% Chebyshev lowpass magnitude response
syms s;
H_BSF = subs(H_LPF, sl, (B * s)/(s^2 + Omega_0^2));
[N, D] = numden(H_BSF);
disp(expand(vpa(N, 2)));
disp(expand(vpa(D, 2)));
%% Chebyshev bandstop magnitude response
syms Omega;
N1 = subs(N* 1e-138, s, 1i * Omega);
D1 = subs(D * 1e-138, s, 1i * Omega);
H_BSF_freq_resp = (N1) / (D1);
H_BSF_freq_resp = vpa(H_BSF_freq_resp, 2);

figure();
fplot(abs(H_BSF_freq_resp));
set(gca, 'XTick', [0.6044, 0.6725, 1.1014, 1.2146], 'xticklabel', {'\Omega_{p1}', '\Omega_{s1}', '\Omega_{s2}', '\Omega_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
hold on;
xline(0.6044, 'magenta-');hold on;
xline(0.6725, 'magenta-');hold on;
xline(1.1014, 'magenta-');hold on;
xline(1.2146, 'magenta-');
hold on;
fplot(s - s - 0.15 + 1, 'k-', 'Markersize', 10);
hold on;
fplot(s-s+1, 'k-', 'Markersize', 10);
hold on;
fplot(s -s + 0.15, 'r-', 'Markersize', 10);
daspect([1 1 1]);
title('Chebyshev Bandstop filter mangitude response');
xlabel('\Omega');
ylabel('|H(j \Omega)|');
axis([0 2 0 1.2]);
%% Chebyshev bandstop phase response
figure();
fplot(angle(H_BSF_freq_resp));
daspect([1 1 1]);
title('Chebyshev Bandstop filter phase response');
xlabel('\Omega');
ylabel('\angle H(j \Omega)');
%% Bilinear transformation from analog to digital
syms z;
Hz = subs(H_BSF, s, (z - 1) / (z + 1));
[Nz, Dz] = numden(Hz);
disp(expand(vpa(Nz, 2)));
disp(expand(vpa(Dz, 2)));
Nzz = sym2poly(Nz);
Dzz = sym2poly(Dz);
[H,f] = freqz(Nzz,Dzz,1024*1024, 260e3);
plot(f,abs(H));
axis([0 100e3 0 1.3]);
hold on;
yline(1.15, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.85, 'red--', 'LineWidth', 1.5);
hold on;
yline(0.15, 'red--', 'LineWidth', 1.5);
hold on;
xline(49e3, 'magenta--', 'LineWidth', 1.5);
hold on;
xline(69e3, 'magenta--', 'LineWidth', 1.5);
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j 2\pi f})|');
title('Magnitude Response of the Discrete Time Bandstop Filter');
set(gca, 'XTick', [45e3, 49e3, 69e3, 73e3], 'xticklabel', {'f_{p1}', 'f_{s1}', 'f_{s2}', 'f_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
%% Phase response of DT BSF
figure();
plot(f, angle(H));
xlim([1e4 110e3]);
ylabel('\angle H(e^{j 2\pi f})');
xlabel('f in Hz');
title('Phase Response of the Discrete Time Bandstop Filter');
hold on;
xline(49e3, 'magenta--', 'LineWidth', 1);
hold on;
xline(69e3, 'magenta--', 'LineWidth', 1);
set(gca, 'XTick', [1e4, 45e3, 49e3, 69e3, 73e3, 110e3], 'xticklabel', {1e4 ,'f_{p1}', 'f_{s1}', 'f_{s2}', 'f_{p2}', 110e3});
%%
% fvtool(Nzz, Dzz);
