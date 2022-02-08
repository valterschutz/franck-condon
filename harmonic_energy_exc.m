function E = harmonic_energy_exc(n)
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

k=2*De*a^2;
w=sqrt(k/mu_I2);

E = hbar*w*(n+1/2);
end
