%% RUN THIS FIRST
% Universal constants
u = 1.66e-27;  % kg
eV = 1.602e-19;  % J
h = 6.626e-34;  % J/Hz
hbar = h/(2*pi);  % J/Hz
kb = 1.380649e-23;  % J/K, Boltzman
c = 299792458;  % m/s

% Constants for iodine
mu_I2 = 126.90447/2*u;  % kg
re_ground = 2.666e-10;  % m

% Run program for several different re_exc. Original is 3.024e-10.
re_exc_og = 3.024e-10;  % m
delta_re_exc = (re_exc_og-re_ground)/10;  % 20% difference
re_exc = linspace(re_exc_og-delta_re_exc, re_exc_og+delta_re_exc, 5);  % List of re_exc values
we_xe_ground = 0.614*100;  % m-1
we_xe_exc = 0.764*100;  % m-1
we_ground = 214.50*100; % m-1
we_exc = 125.69*100;  % m-1
De_ground = h*c*we_ground^2/(4*we_xe_ground);  % J
De_exc = h*c*we_exc^2/(4*we_xe_exc);  % J
a_ground = we_ground*2*pi*c*sqrt(mu_I2/(2*De_ground));  % m-1
a_exc = we_exc*2*pi*c*sqrt(mu_I2/(2*De_exc));  % m-1
lmb_ground = sqrt(2*mu_I2*De_ground)/(a_ground*hbar);
lmb_exc = sqrt(2*mu_I2*De_exc)/(a_exc*hbar);

electronic_energy = inverse_cm_to_J(15769.01);
T = 100;  % Temperature
laser_wavelength = 612e-9;
laser_energy = h*c/laser_wavelength;

% upper_limit_ground = floor(lmb_ground-1/2);
upper_limit_ground = 30;
% upper_limit_exc = floor(lmb_exc-1/2);
upper_limit_exc = 15;

lower_limit_ground = 0;
lower_limit_exc = 0;

N = 2048;
a = 1e-10; b = 9e-10;
dr = (b-a)/N;
r=linspace(a,b,N);  % Seems to fit everything
y1=0*r;
y2=0*r;
%% DO NOT RUN THIS, EXCEPT WHEN CALCULATING NEW INTEGRALS (TAKES TIME)
for m=1:length(re_exc)
    overlap0 = zeros(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1);
    for j=lower_limit_ground:upper_limit_ground
        for k=lower_limit_exc:upper_limit_exc
            y1 = morse_psi_ground(r,j,dr);
            y2 = morse_psi_exc(r,k,dr, re_exc(m));
            % Discrete integral
            val = sum(y1.*y2)*dr;
            fprintf("j=%d, k=%d, val=%f\n",j,k,val)
            overlap0(j+1,k+1)=val;  % One of them should be conjugated but this is not necessary since they are both real
        end
    end
    overlap = overlap0.^2;  % This is what we are interested in
    filename = sprintf("overlap_data_re_%d", m);
    save(filename, 'overlap')
end