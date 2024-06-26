clear all; close all; clc

%% Read data from CSV file
T = readtable('p_lock_offset.txt');

% Convert table to numerical array
X = table2array(T);
X = X(:, 2:end-1);

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
for i = 1:size(X, 2)
    subplot(size(X, 2), 1, i);
    plot(t, X(:, i), 'LineWidth', 1.2);
    title(['Impulse Response for ', subplot_titles{i}]);
    xlabel('Time (seconds)');
    ylabel('SPL (dB)');
    grid on;
end

%% Perform SVD
[U,S,V] = svd(X,'econ');

% Now you can proceed with the rest of your computations...
X2 = table2array(T); % time shifted X'
X2 = X2(:,3:end);

%%  Compute DMD (Phi are eigenvectors)
r = 4;  % truncate at 4 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%%  Plot DMD spectrum
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%% Visualize DMD dominant modes using Atilde and time series

% X_k+1 = Atilde(i)*X
X_k1 = Atilde(1,1)*t'
X_k2 = Atilde(2,1)*t'

figure
plot(t,X_k1,t,X_k2)

