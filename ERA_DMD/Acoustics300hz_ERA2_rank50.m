clear all; close all; clc

%% Read data from simulation output
T = readtable('300hz_single_Quad_Sinc_timeImpulse_1_edited_1.xlsx', 'HeaderLines', 1);
% Convert table to numerical array
X = table2array(T);
total_time = X(10000,1);
yFull = X(:, 2:5);
YY = permute(yFull,[2 3 1]);

%% Plot impulse for each output location
num_points = 10000; % total number of data points
time_step = total_time / (num_points - 1); % calculate time step

% Define time vector
t = 0:time_step:total_time;
t = t'; % transpose into column vector

% Create a figure with a larger size
figure('Position', [500, 100, 1000, 700]); % Adjust position and size as needed

% Define subplot titles
subplot_titles = {'Origin', 'Inlet', 'Top Wall', 'Outlet'};

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
numOutputs = 4;
r = 25; %choose rank to truncate at

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
num_points = 10000;

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


%% Plot impulse responses for sysERA approximation
% Define number of points for impulse response
num_points = 1951;

% Compute impulse response using sysERA
[y2, t2] = impulse(sysERA, num_points);

% Define the number of outputs
num_outputs = size(yFull, 2);

% Ensure yFull and y2 have the same number of rows
min_length = min(size(yFull, 1), size(y2, 1));
yFull_trunc = yFull(1:min_length, :);
y2_trunc = y2(1:min_length, :);

% Titles for the subplots
subplot_titles = {'Inlet', 'Outlet'};

% Create a figure with a larger size
figure('Position', [100, 50, 1200, 800]); % Adjust position and size as needed

% Define positions for the subplots to make them larger
subplot_positions = [
    0.1, 0.55, 0.85, 0.4; % [left, bottom, width, height] for first subplot
    0.1, 0.05, 0.85, 0.4  % [left, bottom, width, height] for second subplot
];

% Iterate over each output and plot its impulse response
for i = 1:2
    % Calculate R² value
    y_actual = yFull_trunc(:,i);
    y_predicted = y2_trunc(:,i);
    ss_res = sum((y_actual - y_predicted).^2);
    ss_tot = sum((y_actual - mean(y_actual)).^2);
    r_squared = 1 - (ss_res / ss_tot);
    
    % Plot the impulse response for the current output
    subplot('Position', subplot_positions(i, :));
    plot(y_actual, 'LineWidth', 2);
    hold on;
    plot(y_predicted,'--','Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
    grid on;
    title(subplot_titles{i}, 'FontSize', 12, 'FontWeight', 'bold');
    legend(['Model (R² = ', num2str(r_squared, '%.2f'), ')'], 'ERA', 'FontSize', 12);
    ylabel('SPL (dB)', 'FontSize', 12);
    xlabel('Iteration', 'FontSize', 12); % Add X label to each subplot
end













