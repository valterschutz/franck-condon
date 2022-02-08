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
we_xe_ground = inverse_cm_to_J(0.614);  % J
we_xe_exc = inverse_cm_to_J(0.764);  % J
we_ground = inverse_cm_to_J(214.50); % J
we_exc = inverse_cm_to_J(125.69);  % J
De_ground = we_ground^2/(4*we_xe_ground) - we_ground/2 + we_xe_ground/4;
De_exc = we_exc^2/(4*we_xe_exc) - we_exc/2 + we_xe_exc/4;
a_ground = sqrt(2/we_ground);
a_exc = sqrt(2/we_exc);
lmb_ground = sqrt(2*mu_I2*De_ground)/(a_ground*hbar);
lmb_exc = sqrt(2*mu_I2*De_exc)/(a_exc*hbar);

electronic_energy = inverse_cm_to_J(15769.01);
T = 100;  % Temperature
sig = 1e-9;  % sigma for gauss distribution
laser_wavelength = 612e-9;
laser_energy = h*c/laser_wavelength;

vibration_energy = laser_energy - electronic_energy;  % Energy left for vibration

% abs(vibration_energy - (energy_exc(0:32)-energy_ground(0)))  % This shows
% that the vibrational quantum number must be 0 in the excited state


% upper_limit_ground = floor(lmb_ground-1/2);
% upper_limit_exc = floor(lmb_exc-1/2);

lower_limit_ground = 0;
lower_limit_exc = 0;
upper_limit_ground = 50;
upper_limit_exc = 0;



N = 500;
a = 2e-10; b = 3.5e-10;
dx = (b-a)/N;
x=linspace(a,b,N);  % Seems to fit everything

%%
overlap = zeros(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1);
for j=lower_limit_ground:upper_limit_ground
    for k=lower_limit_exc:upper_limit_exc
        fprintf("j=%d, k=%d\n",j,k)
        y1 = harmonic_psi_ground(x,j);
        y2 = harmonic_psi_exc(x,k);
        % Discrete integral
        overlap(j+1,k+1)=sum(conj(y1).*y2)*dx;
    end
end
overlap = overlap.^2;
% overlap = overlap .* exp(-harmonic_energy_ground(0:upper_limit_exc)'/(kb*(273.15+T)));
%%
clf
subplot(1,2,1)
imagesc([0 upper_limit_exc], [lower_limit_ground upper_limit_ground], overlap)
title("Probability of transition")
xlabel("Vibrational mode of excited state")
ylabel("Vibrational mode of ground state")
colorbar
%%
energy_difference = electronic_energy*ones(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1);
for j=lower_limit_ground:upper_limit_ground
    for k=0:upper_limit_exc
        fprintf("j=%d, k=%d\n",j,k)
        energy_difference(j+1,k+1) = energy_difference(j+1,k+1) + (harmonic_energy_exc(k) - harmonic_energy_ground(j));
    end
end

subplot(1,2,2)
imagesc([0 upper_limit_exc], [lower_limit_ground upper_limit_ground], energy_difference)
title("Energy difference E_{exc}(k) - E_{ground}(j)")
xlabel("k")
ylabel("j")
colorbar

%%
clf
flat_energy = reshape(energy_difference, [1 numel(energy_difference)]);

flat_overlap = reshape(overlap, [1 numel(overlap)]);
% flat_overlap = flat_overlap / max(flat_overlap);
wavelength = energy_to_m(flat_energy);


x = linspace(min(wavelength),max(wavelength),1000);
y = 0*x;

for j=1:length(wavelength)
    y = y + flat_overlap(j)*gauss(x,wavelength(j),sig);
end

y = y/max(y);
y(x<610e-9) = 0;

plot(x,y)
% axis([600e-9 900e-9 0 1])
%% For testing
clf
subplot(2,1,1)
plot(x,harmonic_psi_ground(x,4)), hold on
plot(x,harmonic_psi_ground(x,5))
% plot(x,harmonic_psi_ground(x,6))
title("Ground state")
subplot(2,1,2)
plot(x,harmonic_psi_exc(x,10)), hold on
% plot(x,harmonic_psi_exc(x,1))
title("Excited state")
% plot(x,harmonic_psi_exc(x,10))