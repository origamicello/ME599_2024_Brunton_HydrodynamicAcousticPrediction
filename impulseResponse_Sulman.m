clear all
close all

% Read data from CSV file
T = readtable('p_lock_offset.txt', 'HeaderLines', 1);
% Convert table to numerical array
yFull = table2array(T);

%% 
total_time = 0.3902; % total time in seconds
num_points = 1951; % total number of data points
time_step = total_time / (num_points - 1); % calculate time step

% Define time vector
t = 0:time_step:total_time;

% Plot impulse response for each output
for i = 2:size(yFull, 2) % When you use size(yFull, 2), it returns the size of yFull along the second dimension, 
                         % which is typically the number of columns. 
    subplot(size(yFull, 2), 1, i);
    plot(t, yFull(:, i));
    title(['Impulse Response for Output ', num2str(i)]);
    xlabel('Time (seconds)');
    ylabel('Response');
    grid on;
end

%% Compare ERA with impulse response
mco = floor((length(yFull)-1)/2);
numInputs = yFull(:,1);
numOutputs = yFull(:,2:6);
r = 5;
[Ar,Br,Cr,Dr,HSVs] = ERA(yFull,mco,mco,numInputs,numOutputs,r);
sysERA = ss(Ar,Br,Cr,Dr,-1);
