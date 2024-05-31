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

[U,S,V] = svd(X2,'econ'); %computes dominaant coherent structures, my POD modes
%what i can do with my SVD to decide how many modes to keep how many to throw away

X3 = table2array(T); % Convert table to numerical array
X3 = X3(:,3:end);

r = 3;  % truncate at 3 modes
U = U(:,1:r); %reduced rank
S = S(1:r,1:r); %reduced rank
V = V(:,1:r); %reduced rank
Atilde = U'*X3*V*inv(S); % (1) gives spatial-temporal coherent modes and (2) Eigenvalues which are my time dynamics shows how my modes evolve in time
[W,eigs] = eig(Atilde); %eigenvectors W and eigenvalues eigs of my reduced dynamic operator Atilda
Phi = X3*V*inv(S)*W; % (eigenmodes) gives me big eigenvector of my orignal A matrix(not reduced)
lambda = diag(eigs);
omega = log(lambda)/time_step; %eigenvalues

%% Reconstructs my signal using Phi
n =X(:,4)';
t = 0:time_step:total_time;
a = exp(omega(1,1)*t);
c = mean(Phi(:,1));
b = real((c*a)'\n'); % effectively finds the coefficients b that best reconstruct the signal n using the Dynamic Mode Decomposition (DMD) modes Φ
d = n;

% figure;
% hold on;
% grid on
% plot(t,d,'green')
% xlabel('time');
% ylabel('state');
% title('Reconstruction of Lock response usingn DMD vs COMSOL data')
% scatter(t,X(:,4), 5, 'red', 'filled')
% legend('green line - DMD', 'scatter - COMSOL data');

%% Reconstructs my signal using Phi
n =X(:,4)';
t = 0:time_step:total_time;
a = exp(omega(1,1)*t);
c = mean(Phi(:,1));
b = real((c*a)'\n'); % effectively finds the coefficients b that best reconstruct the signal n using the Dynamic Mode Decomposition (DMD) modes Φ

xdmd = b*c*a;

figure;
hold on;
grid on
plot(t,xdmd,'green')
xlabel('time');
ylabel('state');
title('Reconstruction of Lock response usingn DMD vs COMSOL data')
scatter(t,X(:,4), 5, 'red', 'filled')
legend('green line - DMD', 'scatter - COMSOL data');

%% Augmented DMD

Xaug = [n(1:end-2);
    n(2:end-1)];
Xaug2 = [n(2:end-1);
    n(3:end)];

[U,S,V] = svd(Xaug,'econ');
Atilde = U'*Xaug2*V*inv(S);
[W,Lambda] = eig(Atilde);
Omega = diag(log(diag(Lambda)))/time_step
Phi = Xaug2*V*inv(S)*W;

b = Phi\Xaug(:,1);
for k=1:length(t)
    xaugdmd(:,k) = Phi*exp(Omega*t(k))*b;
end
plot(t,real(xaugdmd(1,:)),'b--','LineWidth',1.5)
%ylim([-1 1]), grid on
legend('Data','DMD','Augmented DMD')
xlabel('Time'), ylabel('State')


%set(gcf,'Position',[100 100 500 200])
%set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', 'standingwave');