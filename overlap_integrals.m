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
re_exc = 3.024e-10;  % m

upper_limit_ground = 30;
upper_limit_exc = 30;

lower_limit_ground = 0;
lower_limit_exc = 0;

N = 2048;
a = 1e-10; b = 9e-10;
dr = (b-a)/N;
r=linspace(a,b,N);  % Seems to fit everything
y1=0*r;
y2=0*r;
%% ONLY RUN TO CALCULATE NEW INTEGRALS (TAKES TIME)
overlap0 = zeros(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1);
for j=lower_limit_ground:upper_limit_ground
    for k=lower_limit_exc:upper_limit_exc
        y1 = morse_psi_ground(r,j,dr);
        y2 = morse_psi_exc(r,k,dr, re_exc);
        % Discrete integral
        val = sum(y1.*y2)*dr;
        fprintf("j=%d, k=%d, val=%f\n",j,k,val)
        overlap0(j+1,k+1)=val;  % One of them should be conjugated but this is not necessary since they are both real
    end
end
overlap = overlap0.^2;  % This is what we are interested in
filename = "data/overlap_data";
save(filename, 'overlap')