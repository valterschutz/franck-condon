data = readmatrix('data/output4.csv');
data(:,1) = data(:,1) * 1e-9;
data(:,2) = data(:,2) - min(data(:,2));
data(:,2) = data(:,2) / max(data(:,2));
plot(data(:,1),data(:,2))