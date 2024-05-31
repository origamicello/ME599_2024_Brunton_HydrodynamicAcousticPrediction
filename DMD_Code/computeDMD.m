% DMD is best fit linear regression
clear all
close all
clc

% Read data from CSV file
T = readtable('p_lock_offset.txt', 'HeaderLines', 1);
T = T(:, 2:6);

% Convert table to numerical array
X = table2array(T);
X1=X(:,:);
X2=X(:, 2:end-1);

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
for i = 1:size(X1, 2)
    subplot(size(X1, 2), 1, i);
    plot(t, X1(:, i), 'LineWidth', 1.2);
    title(['Impulse Response for ', subplot_titles{i}]);
    xlabel('Time (seconds)');
    ylabel('SPL (dB)');
    grid on;
end
%% Perform SVD
[U,S,V] = svd(X2,'econ'); %computes dominaant coherent structures, my POD modes
semilogy(diag(S), 'linewidth', 2) % good idea to plot my singular values, how they look like, you will see values in pairs (sine and cosine pair) which are the harmonics mode pair for my acosutics signal
xlabel('modes')
%what i can do with my SVD to decide how many modes to keep how many to throw away
%% Create DMD data matrices
X3 = table2array(T); % Convert table to numerical array
X3 = X3(:,3:end);

%%  Compute DMD (Phi are eigenvectors)
r = 3;  % truncate at 3 modes
U = U(:,1:r); %reduced rank
S = S(1:r,1:r); %reduced rank
V = V(:,1:r); %reduced rank
Atilde = U'*X3*V*inv(S); % (1) gives spatial-temporal coherent modes and (2) Eigenvalues which are my time dynamics shows how my modes evolve in time
[W,eigs] = eig(Atilde); %eigenvectors W and eigenvalues eigs of my reduced dynamic operator Atilda
Phi = X3*V*inv(S)*W; % (eigenfunction) gives me big eigenvector of my orignal A matrix(not reduced)

%%  Plot DMD spectrum
figure;
lambdaStr = 'λ';
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
xlabel(['Real ' lambdaStr]) % properly concatenating the string
ylabel(['Imag ' lambdaStr]) % properly concatenating the string
hold on
grid off
scatter(real(diag(eigs)),imag(diag(eigs)),'or', 'linewidth', 2)
axis([-1.1 1.1 -1.1 1.1]);


%% DMD Spectra
lambda = diag(eigs);
omega = log(lambda)/time_step;
plot(omega, '+', 'LineWidth', 4);
%% Compute DMD Solution (1) (tells you how these eigenvalues evolve in time)
xi = linspace(-1, 1, size(X, 1)); % Adjust as necessary
x1 = X(:, 1);
b = Phi\x1;
time_dynamics = zeros(r,length(t));
for iter = 1:length(t),
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end;
X_dmd = Phi*time_dynamics; %eigenvalues for the big matrix

subplot(2,2,4); 
surfl(real(X_dmd')); 
shading interp; colormap(gray); view(-20,60);
set(gca, 'YTick', numel(t)/4 * (0:4)), 
set(gca, 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca, 'XTick', linspace(1,numel(xi),3)), 
set(gca, 'Xticklabel',{'1', '976', '1951'});

set(gcf, 'Color', 'w', 'Position', [500 500 400 300]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 3], 'PaperPositionMode', 'manual');

%% Compute DMD Solution (2) (plot each row of X_dmd as a line:)
% This method plots each row of the X_dmd matrix as a line, showing how the solution evolves for different values of xi.

xi = linspace(0, total_time, size(X, 1)); % Adjust as necessary
x1 = X(:, 1);
b = Phi\x1;
time_dynamics = zeros(r, length(t));
for iter = 1:length(t)
    time_dynamics(:, iter) = (b .* exp(omega * t(iter)));
end
X_dmd = Phi * time_dynamics;

% Plot each row of X_dmd as a separate line
plot(xi, real(X_dmd), 'LineWidth', 1.2);
xlabel('Time (sec)');
ylabel('Real(X_{dmd})');
title('DMD Modes');
set(gca, 'XTick', linspace(-1, 1, 3));
set(gca, 'YTick', linspace(0, size(X_dmd, 1), 5));
set(gca, 'XTick', linspace(0, total_time, 5));
grid off;
hold off;

% Create a legend for each mode
legend_labels = arrayfun(@(k) sprintf('Mode %d', k), 1:r, 'UniformOutput', false);
legend(legend_labels, 'Location', 'Best');

set(gcf, 'Color', 'w', 'Position', [500 500 400 300]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 3], 'PaperPositionMode', 'manual');
%% FFT of the Strong Modes
% Number of columns to compute FFT for
num_columns = 3; % Number of DMD modes solution = 1951, but only three modes have the largest power

% Initialize figure
figure;
hold on;

% Sampling frequency
Fs = 1 / time_step;

% Number of samples
N = length(t);

% Frequency vector
freqs = Fs * (0:(N/2)) / N;

% Loop through each column and compute FFT
for col_idx = 1:num_columns
    % Compute FFT of the real part of the current DMD mode
    xhat = fft(real(X_dmd(:, col_idx)));
    
    % Compute power spectrum
    xpower = abs(xhat(1:N/2+1)) * 2 / N;
    
    % Plot power spectrum
    plot(freqs, xpower, 'LineWidth', 1.2, 'DisplayName', sprintf('Mode %d', col_idx));
end

xlabel('Frequency (Hz)');
ylabel('Power');
title('FFT of DMD Modes');
xlim([0, 1000]);
grid on;
hold off;

% Add legend
legend('show', 'Location', 'Best');
set(gcf, 'Color', 'w', 'Position', [500 500 400 300]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 3], 'PaperPositionMode', 'manual');

%% Spectogram of Strongest Modes
% Given data
total_time = 0.3902; % total time in seconds
num_points = 1951; % total number of data points

% Calculate time step
dt = total_time / (num_points - 1);

% Create time vector
t = (0:dt:total_time)';

% Sampling frequency
Fs = 1 / dt;

% Number of samples
N = length(t);

% Frequency vector
freqs = Fs * (0:(N/2)) / N;

% Assume y is the data vector X
y = real(X_dmd(:, 1:3));

% Fourier Transform
FT = fft(y);

% Define sampling rate and other parameters for wavelet transform
fs = 10000; % sampling rate
w = 3.5; % width of the mother wavelet in the frequency domain

% Create frequency vector
freq = linspace(1, fs/2, 3000);

% Compute widths for the Morlet wavelet
widths = (w * fs) ./ (2 * freq * pi);

% Perform Continuous Wavelet Transform (CWT)
cwtm = cwt(y, widths, 'morl', 'SamplingFrequency', fs);

% Compute the magnitude of the coefficients
magnitude = abs(cwtm);

% Plot the spectrogram in 2D
figure('Position', [100, 100, 800, 400]);
imagesc(t, freq, magnitude);
axis xy;
h = colorbar; % Get the handle for the colorbar
colormap jet;
ylabel(h, 'Magnitude of wavelet coeffecients'); % Label the colorbarcolormap jet;
caxis([min(magnitude(:)), max(magnitude(:))]);
title('Wavelet Transform Spectrogram of all 3 modes - Lock');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0, 500]);
%% FFT of time impulse
% Number of columns to compute FFT for
num_columns = 5;

% Initialize figure
figure;
hold on;

% Sampling frequency
Fs = 1/time_step;

% Number of samples
N = length(t);

% Frequency vector
freqs = Fs * (0:(N/2)) / N;

% Loop through each column and compute FFT
for col_idx = 1:num_columns
    % Compute FFT of the current column
    xhat = fft(X(:, col_idx));
    
    % Compute power spectrum
    xpower = abs(xhat(1:N/2+1)) * 2 / N;
    
    % Plot power spectrum
    plot(freqs, xpower, 'LineWidth', 1.2, 'DisplayName', sprintf('Output %d', col_idx));
end

% Set x-axis limits
xlim([0 600]);

% Set other axis properties
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Outputs');
grid on;
hold off;

% Add legend
legend('show', 'Location', 'Best');
set(gcf, 'Color', 'w', 'Position', [500 500 400 300]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 3], 'PaperPositionMode', 'manual');

%% Compute Spectogram B/w time and frequency
% Given data
total_time = 0.3902; % total time in seconds
num_points = 1951; % total number of data points

% Calculate time step
dt = total_time / (num_points - 1);

% Create time vector
t = (0:dt:total_time)';

% Assume y is the data vector X
y = X(:,1);

% Fourier Transform
FT = fft(y);

% Define sampling rate and other parameters for wavelet transform
fs = 5000; % sampling rate
w = 3.5; % width of the mother wavelet in the frequency domain

% Create frequency vector
freq = linspace(1, fs/2, 200);

% Compute widths for the Morlet wavelet
widths = (w * fs) ./ (2 * freq * pi);

% Perform Continuous Wavelet Transform (CWT)
cwtm = cwt(y, widths, 'morl', 'SamplingFrequency', fs);

% Compute the magnitude of the coefficients
magnitude = abs(cwtm);

% Scale to dB and clip values
p0 = 1.0e-6;
log_scaling_SPL = (20.0 * log10(sqrt(magnitude))) / p0;

% Plot the spectrogram in 2D
figure('Position', [100, 100, 800, 400]);
imagesc(t, freq, log_scaling_SPL);
axis xy;
colorbar;
colormap jet;
caxis([min(log_scaling_SPL(:)), max(log_scaling_SPL(:))]);
title('Wavelet Transform Spectrogram, Exact - Source');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0, 1000]);
% grid on; % Uncomment if you want to show the grid lines

