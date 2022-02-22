clf
data = readmatrix('data/output1.csv');
data(:,1) = data(:,1) * 1e-9;
data(:,2) = data(:,2) - min(data(:,2));
data(:,2) = data(:,2) / max(data(:,2));
plot(data(:,1),data(:,2)), hold on
title("Experimental data")
xlabel("Wavelength [m]")
ylabel("Intensity")