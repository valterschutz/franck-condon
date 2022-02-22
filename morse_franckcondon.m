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

% vibration_energy = laser_energy - electronic_energy;  % Energy left for vibration

% abs(vibration_energy - (energy_exc(0:32)-energy_ground(0)))  % This shows
% that the vibrational quantum number must be 0 in the excited state


% upper_limit_ground = floor(lmb_ground-1/2);
upper_limit_ground = 70;
% upper_limit_exc = floor(lmb_exc-1/2);
upper_limit_exc = 5;

lower_limit_ground = 0;
lower_limit_exc = 0;

N = 2048;
a = 1e-10; b = 9e-10;
dr = (b-a)/N;
r=linspace(a,b,N);  % Seems to fit everything
y1=0*r;
y2=0*r;
%%
overlap0 = zeros(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1);
for j=lower_limit_ground:upper_limit_ground
    for k=lower_limit_exc:upper_limit_exc
        fprintf("j=%d, k=%d\n",j,k)
        y1 = morse_psi_ground(x,j,dx);
        y2 = morse_psi_exc(x,k,dx);
        % Discrete integral
        overlap0(j+1,k+1)=sum(y1.*y2)*dx;  % One of them should be conjugated but this is not necessary since they are both real
    end
end
overlap = overlap0.^2;  % This is what we are interested in
% overlap = overlap .* exp(-harmonic_energy_ground(0:upper_limit_exc)'/(kb*(273.15+T)));
%%
clf
load('overlap_data')
load('overlap0_data')
subplot(1,2,1)
% overlap = overlap.*overlap(1,:);
overlap = overlap./max(overlap,[],'all');
imagesc([0 upper_limit_exc], [lower_limit_ground upper_limit_ground], overlap)
% imagesc([0 5], [lower_limit_ground upper_limit_ground], overlap(:,1:5))
title("Probability of transition")
xlabel("Vibrational mode of excited state")
ylabel("Vibrational mode of ground state")
colorbar
%%
energy_difference = electronic_energy*ones(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1-lower_limit_exc);
for j=lower_limit_ground:upper_limit_ground
    for k=lower_limit_exc:upper_limit_exc
        fprintf("j=%d, k=%d\n",j,k)
        energy_difference(j+1,k+1) = energy_difference(j+1,k+1) + morse_energy_exc(k) - morse_energy_ground(j);
    end
end

subplot(1,2,2)
imagesc([0 upper_limit_exc], [lower_limit_ground upper_limit_ground], energy_difference)
% imagesc([0 5], [lower_limit_ground upper_limit_ground], energy_difference(:,1:5))
title("Energy difference E_{exc}(k) - E_{ground}(j)")
xlabel("k")
ylabel("j")
colorbar

%%
clf
overlap = overlap(:,1);
energy_difference = energy_difference(:,1);
flat_energy = reshape(energy_difference, [1 numel(energy_difference)]);
flat_overlap = reshape(overlap, [1 numel(overlap)]);
wavelength = energy_to_m(flat_energy);

% Choose interval and sigma
% interval_start = 0e-6;
% interval_end = 1e-6;
sig = 1e-9;  % sigma for gauss distribution
% 
% wavelength = wavelength(wavelength>interval_start & wavelength<interval_end);
% flat_overlap = flat_overlap(wavelength>interval_start & wavelength<interval_end);


x = linspace(500e-9,1300e-9,10000);
y = 0*x;

for j=1:length(wavelength)
    y = y + flat_overlap(j)*gauss(x,wavelength(j),sig);
end


% y(x<610e-9) = 0;
y = y/max(y);

plot(x,y)
% axis([100e-9 60000e-9 0 1])
%% For testing
clf
subplot(2,1,1)
y = morse_psi_ground(x,70,dx);
plot(x,y)
% plot(x,morse_psi_ground(x,10,dx)), hold off
% plot(x,morse_psi_ground(x,5))
title("Ground state")
subplot(2,1,2)
plot(r,morse_psi_exc(x,5,dr)), hold on
% plot(x,morse_psi_exc(x,1,dx))
% plot(x,morse_psi_exc(x,2,dx)), hold off
title("Excited state")