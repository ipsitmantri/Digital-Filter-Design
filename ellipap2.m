function [z,p,H0,B,A] = ellipap2(N,Ap,As)
if nargin==0, help ellipap2; return; end
Gp = 10^(-Ap/20); % passband gain
ep = sqrt(10^(Ap/10) - 1); % ripple factors
es = sqrt(10^(As/10) - 1);
k1 = ep/es;
k = ellipdeg(N,k1); % solve degree equation
L = floor(N/2); r = mod(N,2); % L is the number of second-order sections
i = (1:L)'; ui = (2*i-1)/N; zeta_i = cde(ui,k); % zeros of elliptic rational function
z = j./(k*zeta_i); % filter zeros = poles of elliptic rational function
v0 = -j*asne(j/ep, k1)/N; % solution of sn(jv0NK1, k1) = j/?p
p = j*cde(ui-j*v0, k); % filter poles
p0 = j*sne(j*v0, k); % first-order pole, needed when N is odd
B = [ones(L,1), -2*real(1./z), abs(1./z).^2]; % second-order numerator sections
A = [ones(L,1), -2*real(1./p), abs(1./p).^2]; % second-order denominator sections
if r==0, % prepend first-order sections
B = [Gp, 0, 0; B]; A = [1, 0, 0; A];
else
B = [1, 0, 0; B]; A = [1, -real(1/p0), 0; A];
end
z = cplxpair([z; conj(z)]); % append conjugate zeros
p = cplxpair([p; conj(p)]); % append conjugate poles
if r==1, p = [p; p0]; end % append first-order pole when N is odd
H0 = Gp^(1-r); % dc gain
end