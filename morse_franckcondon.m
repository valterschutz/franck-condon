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

laser_wavelength = 612e-9;
laser_energy = h*c/laser_wavelength;
vibration_energy = laser_energy - electronic_energy;  % Energy left for vibration

% vibration_energy - (morse_energy_exc(0:32)-morse_energy_ground(0))  % Exactly enough energy for fifth vibrational level in excited state 


% upper_limit_ground_theory = floor(lmb_ground-1/2);
upper_limit_ground = 70;
% upper_limit_exc_theory = floor(lmb_exc-1/2);
upper_limit_exc = 15;

lower_limit_ground = 0;
lower_limit_exc = 0;

N = 2048;
a = 1e-10; b = 9e-10;
dr = (b-a)/N;
r=linspace(a,b,N);  % Seems to fit everything
y1=0*r;
y2=0*r;
%% RUN THIS
clf
load('overlap_data/overlap_data_final')
subplot(1,2,1)
overlap = overlap./max(overlap,[],'all');
imagesc([lower_limit_exc upper_limit_exc], [lower_limit_ground upper_limit_ground], overlap)
title("Probability of transition",'fontsize',16)
xlabel("v'",'fontsize',14)
ylabel("v''",'fontsize',14)
colormap(gca,cool)
colorbar

energy_difference = electronic_energy*ones(upper_limit_ground+1-lower_limit_ground, upper_limit_exc+1-lower_limit_exc);
for j=lower_limit_ground:upper_limit_ground
    for k=lower_limit_exc:upper_limit_exc
        energy_difference(j+1,k+1) = energy_difference(j+1,k+1) + morse_energy_exc(k) - morse_energy_ground(j);
    end
end

subplot(1,2,2)
imagesc([0 upper_limit_exc], [lower_limit_ground upper_limit_ground], energy_difference)
title("\Delta E = E_{electronic} + E_{exc}(v') - E_{ground}(v'') [eV]",'fontsize',16)
xlabel("v'",'fontsize',14)
ylabel("v''",'fontsize',14)
colormap(gca,hot)

energy_labels = (0.5:0.5:2);
energy_ticks = energy_labels*eV;
colorbar('Ticks',energy_ticks, 'TickLabels', energy_labels)

% Plot energies that are experimentally measured (output1.csv)
% exp_wavelengths = [5.81746e-07, 5.89153e-07, 5.96184e-07, 6.03996e-07, 6.11903e-07, 6.19653e-07, 6.27996e-07, 6.36465e-07, 6.44871e-07,6.53434e-07, 6.62496e-07, 6.71465e-07, 6.80371e-07, 6.8984e-07, 6.92653e-07, 6.99715e-07, 7.02965e-07, 7.06528e-07, 7.09715e-07, 7.24153e-07, 7.29778e-07, 7.4009e-07, 7.50809e-07, 7.61746e-07, 7.72903e-07, 7.77684e-07, 7.84246e-07];
% exp_intensities = [0.014068, 0.0153546, 0.0162136, 0.035622, 1, 0.162708, 0.21963, 0.166461, 0.0365643, 0.242126, 0.0957773, 0.0896708, 0.0163895, 0.0285911, 0.0152358, 0.0759446, 0.0147762, 0.0343232, 0.0375311, 0.0186084, 0.0293969, 0.0361266, 0.0129227, 0.0109756, 0.0226675, 0.0178466, 0.0186132];
exp_wavelengths = [5.81746e-07, 5.89153e-07, 5.96184e-07, 6.03996e-07, 6.19653e-07, 6.27996e-07, 6.36465e-07, 6.44871e-07,6.53434e-07, 6.62496e-07, 6.71465e-07, 6.80371e-07, 6.8984e-07, 6.92653e-07, 6.99715e-07, 7.02965e-07, 7.06528e-07, 7.09715e-07, 7.24153e-07, 7.29778e-07, 7.4009e-07, 7.50809e-07, 7.61746e-07, 7.72903e-07, 7.77684e-07, 7.84246e-07];
exp_intensities = [0.014068, 0.0153546, 0.0162136, 0.035622, 0.162708, 0.21963, 0.166461, 0.0365643, 0.242126, 0.0957773, 0.0896708, 0.0163895, 0.0285911, 0.0152358, 0.0759446, 0.0147762, 0.0343232, 0.0375311, 0.0186084, 0.0293969, 0.0361266, 0.0129227, 0.0109756, 0.0226675, 0.0178466, 0.0186132];
exp_intensities = exp_intensities / max(exp_intensities);
exp_energies = h*c./exp_wavelengths;

% subplot(1,3,3)
% eps = 1e-20;
% for j=1:length(exp_energies)
%     A = abs(energy_difference - exp_energies(j));
%     A(A > min(A,[],'all')+eps) = 0;
%     A = A./max(A,[],'all') * log(1+exp_intensities(j));
%     imagesc([0 upper_limit_exc], [lower_limit_ground upper_limit_ground], A, 'AlphaData', A), hold on
% end
% hold off
% title("Measured intensity")
% xlabel("Predicted k")
% ylabel("Predicted j")


%% AND FINALLY THIS
clf

% Choose interval and sigma
interval_start = min(exp_wavelengths);
interval_end = max(exp_wavelengths);
sig = 2e-9;  % sigma for gauss distribution
x = linspace(interval_start,interval_end,10000);
y = 0*x;

% First plot experimental data

for j=1:length(exp_wavelengths)
    y = y + exp_intensities(j)*gauss(x,exp_wavelengths(j),sig);
end

y = y/max(y);
subplot(2,1,1)
% bar(exp_wavelengths, exp_intensities)
plot(x,y)
title("Experimental data")
xlabel("Wavelength [m]")
ylabel("Intensity")
axis([interval_start interval_end 0 1])

% Now plot theoretical data
% Change this to look at radiation from different vibrational levels
overlap_plot = overlap(:,5); energy_difference_plot = energy_difference(:,5);
flat_energy = reshape(energy_difference_plot, [1 numel(energy_difference_plot)]);
flat_overlap = reshape(overlap_plot, [1 numel(overlap_plot)]);
wavelength = energy_to_m(flat_energy);




y = 0*x;

for j=1:length(wavelength)
    y = y + flat_overlap(j)*gauss(x,wavelength(j),sig);
end

y = y/max(y);

subplot(2,1,2)
plot(x,y)
xlabel("Wavelength [m]")
ylabel("Intensity")
title("Expected intensity from Franck-Condon theory")
axis([interval_start interval_end 0 1])

calibration = h*c/inverse_cm_to_J(15797.976);  % Should be strongest line
% xline(calibration,'r','linewidth',2)
% xline(h*c/(electronic_energy+morse_energy_exc(6)-morse_energy_ground(3)),'linewidth',2)
hold off
%% For testing
clf
re_exc_og = 3.024e-10;  % m
subplot(2,1,1)
plot(r,morse_psi_ground(r,1,dr),':','linewidth',2), hold on
plot(r,morse_psi_ground(r,3,dr),'linewidth',1), hold off
legend("v''=0","v''=3")
title("Ground state", 'fontsize',16)
ylabel('\psi','fontsize',16, 'fontweight','bold')
grid on
xlim([1e-10 4e-10])
xticks((1:0.5:4)*1e-10)
xticklabels([])
yticklabels([])
subplot(2,1,2)
plot(r,morse_psi_exc(r,4,dr,re_exc_og),':','linewidth',2), hold on
plot(r,morse_psi_exc(r,6,dr,re_exc_og),'linewidth',1),hold off
legend("v'=4","v'=6")
% plot(x,morse_psi_exc(x,2,dx)), hold off
title("Excited state",'fontsize',16)
grid on
xlim([1e-10 4e-10])
ylabel('\psi','fontsize',16, 'fontweight','bold')
xlabel('Internuclear distance [Å]')
xticks((1:0.5:4)*1e-10)
xticklabels(1:0.5:4)
yticklabels([])