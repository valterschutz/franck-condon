function y = harmonic_psi_exc(x, n)
% Universal constants
u = 1.66e-27;  % kg
h = 6.626e-34;  % J/Hz
hbar = h/(2*pi);  % J/Hz

% Constants for iodine
mu_I2 = 126.90447/2*u;  % kg
re = 3.024e-10;  % m
we_xe = inverse_cm_to_J(0.764);  % J
we = inverse_cm_to_J(125.69);  % J
De = we^2/(4*we_xe) - we/2 + we_xe/4;
a = sqrt(2/we);

x = x-re;

k=2*De*a^2;
w=sqrt(k/mu_I2);
y = 1/sqrt(2^n*factorial(n))*(mu_I2*w/(pi*hbar))^(1/4)*exp(-mu_I2*w*x.^2/(2*hbar)).*hermiteH(n,sqrt(mu_I2*w/hbar)*x);
end
