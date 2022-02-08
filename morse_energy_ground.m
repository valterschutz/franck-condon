function En = morse_energy_ground(n)
% Universal constants
h = 6.626e-34;  % J/Hz
c = 299792458;  % m/s                       %m

% Constants for iodine
we_xe = inverse_cm_to_J(0.614);  % J
we = inverse_cm_to_J(214.50); % J

En=(we*(n+1/2)-we_xe*(n+1/2).^2);
end

