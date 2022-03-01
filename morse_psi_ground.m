function y = morse_psi_ground(r, n, dx)
% Universal constants
u = 1.66e-27;  % kg
h = 6.626e-34;  % J/Hz
hbar = h/(2*pi);  % J/Hz
c = 299792458;  % m/s

SCALING = 8.5e2;

% Constants for iodine
mu_I2 = 126.90447/2*u;  % kg
re = 2.666e-10;  % m
we_xe = 0.614*100;  % m-1
we = 214.50*100;  % m-1
De = h*c*we^2/(4*we_xe);  % J
a = we*2*pi*c*sqrt(mu_I2/(2*De));  % m-1
lmb = sqrt(2*mu_I2*De)/(a*hbar);
x = a*r;
xe = a*re;

z = 2*lmb*exp(xe-x);
y = exp(log(z)*(lmb-n-1/2)-z/2-SCALING).*laguerreL(n, 2*lmb-2*n-1,z);

% Normalization
N = sum(y.^2)*dx;
y = y/sqrt(N);
end

