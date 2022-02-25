clf
data = readmatrix('data/output1.csv');
data(:,1) = data(:,1) * 1e-9;
data(:,2) = data(:,2) - min(data(:,2));
data(:,2) = data(:,2) / max(data(:,2));
subplot(1,2,1)
plot(data(:,1),data(:,2))
title("Experimental data")
xlabel("Wavelength [m]")
ylabel("Intensity")

strong_lines_wavelength = [5119.29, 5161.20, 5245.71, 5338.22, 5345.15, 5435.83, 5464.62, 5625.69, 5690.91, 5710.53, 5950.25,6074.98, 6127.49, 6619.66, 6812.57, 7402.06, 7468.99, 8043.74, 8240.05, 8393.30, 8857.50, 9022.40, 9058.33, 9113.91, 9426.71, 9427.15 9653.06, 9731.73, 10466.54]*1e-10;
strong_lines_intensity = [130,150,150,500,250,150,100,500,100,200,250,100,100,70,200,70,70, 130,50,130,40,70,200,150,50,40,40,70,70];
strong_lines_intensity = strong_lines_intensity / max(strong_lines_intensity);
subplot(1,2,2)
bar(strong_lines_wavelength, strong_lines_intensity)
title("Experimental data from online source")