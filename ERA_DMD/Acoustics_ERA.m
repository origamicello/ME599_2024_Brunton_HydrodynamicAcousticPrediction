clear all; close all; clc

%% Read data from simulation output
T = readtable('p_lock_offset.txt', 'HeaderLines', 1);
% Convert table to numerical array
yFull = table2array(T)
yFull = yFull(:, 2:6);
YY = permute(yFull,[2 3 1]);

%% Plot impulse for each output location
total_time = 0.3902; % total time in seconds
num_points = 1951; % total number of data points
time_step = total_time / (num_points - 1); % calculate time step

% Define time vector
t = 0:time_step:total_time;
t = t'; % transpose into column vector

% Create a figure with a larger size
figure('Position', [500, 100, 1000, 700]); % Adjust position and size as needed

% Define subplot titles
subplot_titles = {'Inlet', 'Outlet', 'Top Wall', 'Bottom Wall', 'Outlet'};

% Plot impulse response for each output
for i = 1:size(yFull, 2)
    subplot(size(yFull, 2), 1, i);
    plot(t, yFull(:, i), 'LineWidth', 1.2);
    title(['Impulse Response for ', subplot_titles{i}]);
    xlabel('Time (seconds)');
    ylabel('SPL (dB)');
    grid on;
end

%% Compute additional terms & ERA
mco = floor((length(yFull)-1)/2);
numInputs = 1;
numOutputs = 5;
r = 100; %choose rank to truncate at

%ERA
[Ar,Br,Cr,Dr,HSVs] = ERA(YY,mco,mco,numInputs,numOutputs,r);
sysERA = ss(Ar,Br,Cr,Dr,-1);

%% Plot Hankel Singular Values (HSVs)
figure;
semilogy(1:numel(HSVs), HSVs, 'r-', 'LineWidth', 1.4);
hold on;

% Plot a scatter point for the location of 'r' on the Hankel plot
scatter(r, HSVs(r), 100, 'black', 'o',LineWidth=2); % Plot a blue filled circle at (r, HSVs(r))

legend('', 'rank r');
title('Hankel Singular Values');
xlabel('Index');
ylabel('Magnitude');
grid on;

%% Plot impulse responses for sysERA approximation
% Define number of points for impulse response
num_points = 1951;

% Compute impulse response using sysERA
[y2, t2] = impulse(sysERA, num_points);

% Define the number of outputs
num_outputs = size(yFull, 2);

% Create a figure with a larger size
figure('Position', [700, 100, 1000, 700]); % Adjust position and size as needed

% Iterate over each output and plot its impulse response
for i = 1:num_outputs
    % Plot the impulse response for the current output
    subplot(num_outputs, 1, i);
    plot(yFull(:,i), 'LineWidth', 1.5);
    hold on;
    plot(y2(:,i), '--r', 'LineWidth', 1);
    grid on;
    title(subplot_titles{i});
    legend('Model', 'ERA');
end

% Set the title for the entire plot
sgtitle(['ERA using rank ', num2str(r), '/', num2str(mco)], 'FontName', 'Times New Roman');














