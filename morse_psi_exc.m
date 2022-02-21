function y = morse_psi_exc(r, n, dx)
% Universal constants
u = 1.66e-27;  % kg
h = 6.626e-34;  % J/Hz
hbar = h/(2*pi);  % J/Hz

% Constants for iodine
mu_I2 = 126.90447/2*u;  % kg
re = 3.024e-10;  % m
we_xe = inverse_cm_to_J(0.764);  % J
we = inverse_cm_to_J(125.69);  % J
% De = we^2/(4*we_xe) - we/2 + we_xe/4;
De = we^2/(4*we_xe);
a = sqrt(2/we);
lmb = sqrt(2*mu_I2*De)/(a*hbar);
x = a*r;
xe = a*re;

z = 2*lmb*exp(xe-x);
y = z.^(lmb-n-1/2).*exp(-z/2).*laguerreL(n, 2*lmb-2*n-1,z);

% Normalization
N = sum(y.^2)*dx;
y = y/sqrt(N);

end

