function En = morse_energy_exc(n)
% Universal constants
h = 6.626e-34;  % J/Hz
c = 299792458;  % m/s                       %m

% Constants for iodine
we_xe = inverse_cm_to_J(0.764);  % J
we = inverse_cm_to_J(125.69);  % J

electronic_energy = inverse_cm_to_J(15769.01);

En=(we*(n+1/2)-we_xe*(n+1/2).^2);
end
