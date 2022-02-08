function J = inverse_cm_to_J(inverse_cm)
% Universal constants
h = 6.626e-34;  % J/Hz
c = 299792458;  % m/s

k = inverse_cm * 10^2;  % Wavenumber k = 1/lambda
J = h*c*k;
end

