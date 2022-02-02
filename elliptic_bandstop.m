close all;
clearvars;
clc;

% m is unique for everyone
% find addend such that specs are met by looking at plot
% remove semicolon to see the value printed

% 1. Specs
m = 55;
q_m = floor(0.1*m);
r_m = m - 10*q_m;
BL = 25+1.9*q_m + 4.1*r_m;
BH = BL + 20;
trans_bw = 4;

% 2. Band Edge specifications
fp1 = BL-trans_bw;
fs1 = BL;
fs2 = BH;
fp2 = BH+trans_bw;
f_samp = 260;
wp1_by_pi=fp1*2/f_samp;
ws1_by_pi=fs1*2/f_samp;
ws2_by_pi=fs2*2/f_samp;
wp2_by_pi=fp2*2/f_samp;

% 3. Transformed Band Edge specs using Bilinear Transformation         
Wp1 = tan(fp1/f_samp*pi);
Ws1 = tan(fs1/f_samp*pi);
Ws2 = tan(fs2/f_samp*pi);
Wp2 = tan(fp2/f_samp*pi);
B1=Wp2-Wp1;
W0=sqrt(Wp1*Wp2);

% 4. Frequency transformation 
WLp1=(B1*Wp1)/(W0^2-Wp1^2);
WLs1=(B1*Ws1)/(W0^2-Ws1^2);
WLs2=(B1*Ws2)/(W0^2-Ws2^2);
WLp2=(B1*Wp2)/(W0^2-Wp2^2);

% 5. Lowpass specs
Ws=min(abs(WLs1),abs(WLs2));
Wp=WLp1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% elliptic low pass

Gp = 0.85; 
Gs = 0.15;
ep = sqrt(1/Gp^2 - 1); 
es = sqrt(1/Gs^2 - 1);
k = Wp/Ws ;
k1 = ep/es;
[K,Kp] = ellipk(k);
[K1,K1p] = ellipk(k1);
Nexact = (K1p/K1)/(Kp/K);
N = ceil(Nexact); 
Ap=1.4116;
As=16.4782;
[z,p,H0,B,A] = ellipap2(N,Ap,As);
% k = ellipdeg(N,k1);
% L = floor(N/2); 
% r = mod(N,2); 
% i = (1:L)';
% u = (2*i-1)/N; 
% zeta_i = cde(u,k); 
% za = Wp * j./(k*zeta_i);
% v0 = -j*asne(j/ep, k1)/N;
% pa = Wp * j*cde(u-j*v0, k);
% pa0 = Wp * j*sne(j*v0, k);
% B = [ones(L,1), -2*real(1./za), abs(1./za).^2]; 
% A = [ones(L,1), -2*real(1./pa), abs(1./pa).^2];
% if r==0
%     B = [Gp, 0, 0; B]; A = [1, 0, 0; A];
% else
%     B = [1, 0, 0; B];
%     A = [1, -real(1/pa0), 0; A]; 
% end
% f = linspace(0,10,2001);
% for n=1:length(f)
%     s = j*2*pi*f(n);
%     H(n) = prod((B(:,1) + B(:,2)*s + B(:,3)*s^2)./...
%                (A(:,1) + A(:,2)*s + A(:,3)*s^2));
% end
% figure
% plot(f,abs(H),"r-");  %%%plotting elliptic low pass filter
% xlim([0,10]); 
% ylim([0,1.1]); 
% grid off;
% set(gca, "xtick", 0:1:10); 
% set(gca, "ytick", 0:0.1:1);
% title("Elliptic Lowpass"); 
% xlabel("f"); 
% ylabel("|H(f)|");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% converting back to bandpass
[y1,~] = size(B);
num = B(1,:);
den = A(1,:);
if y1>1
    for i=2:y1
        num=conv(num,B(i,:));
        den=conv(den,A(i,:));
    end
end

syms s z;

analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((B1*s)/(s*s + W0*W0));        %bandpass transformation
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
kd = ds(1);    
kn = ns(1);    
ds = ds/kd;
ns = ns/kn;
k=kn/kd;

% 8. Discrete bpf
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
kd = dz(1);                                             %normalisation factor
k = dz(1);    
dz = dz/k;
nz = nz/k;
fvtool(nz,dz,'Analysis','freq');                        %frequency response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,p,~]=tf2zp(nz,dz);
figure;
plot(real(p),imag(p),'rX');
title("Poles of the transfer function")
xlabel("Re(z)")
ylabel("Im(z)")
axis equal
grid on
t = linspace(0,2*pi,1000);
hold on
plot(cos(t),sin(t),'b-') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")
xlabel("Hz")
ylabel("|H(f)|")
xline(fs1,'--m');
xline(fp1,'--m');
xline(fp2,'--m');
xline(fs2,'--m');
yline(0.85,'--m');
yline(0.15,'--m');
grid