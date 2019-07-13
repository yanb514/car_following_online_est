% Unscented Kalman filter for parameter estimation.
% Estimate the spacing, velocity and Bando OV parameters of the
% car-following system

clear all
% close all
clc
tic

% plot data
data = csvread('data/VehA_dip_min.csv');
data = data(300:1200,:);
data(:,1) = data(:,1) - data(1,1); % normalize time stamps
% save('data.csv','data')

% figure;
% plot(data(:,1),data(:,2)); % follower speed m/s
% hold on
% plot(data(:,1),data(:,3)); % leader speed m/s
% % plot(data(:,1),data(:,4)); % spacing m
% legend('Follower','Leader');
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% set(gca,'fontsize',20)

dt = 0.1;
T = data(end,1);
timestep = dt:dt:T;
num_tstep = length(timestep);

%% Matrices initialize
w = [2 1.5 0.1 0.1 0.1 0.1]; % process noise stdev
Q = diag(w.^2); % covariance of process noise
v = [0.5 0.2]; % measurement noise stdev
% v = [3 2];
R = diag(v.^2); % covariance of measurement noise
H = [1 0 0 0 0 0; 0 1 0 0 0 0]; % measurement equation
L = eye(6); % covariance of process error
M = eye(2); % covariance of measurement error

theta_0 = [4.3049 26.8558 31.0387 11.9599];
x = [data(1,4); data(1,2); theta_0']; % initial state
P = diag([3 2 2 6 6 4].^2); % initial covariance
% P = diag([0.5 0.3 2 6 6 4].^2); % initial covariance
% Initialize arrays for later plotting
x_prior = x;
x_posterior = x;

u = data(:,3); % input is the lead vehicle's speed
y_obs = [data(:,4)';data(:,2)']; % s and v are measured

n = numel(x);
m = numel(measure(x));
alpha=1e-3;                                 %default, tunable
ki=0;                                       %default, tunable
beta=2;                                     %default, tunable
lambda=alpha^2*(n+ki)-n;                    %scaling factor
c=n+lambda;
Wm=[lambda/c 0.5/c+zeros(1,2*n)];           %weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);
c=sqrt(c);

%% filter
figure;
for k = 2:num_tstep+1
    % Prior update.
    % Generate sigma points
    A = c*chol(P)';
    Y = x(:,ones(1,numel(x))); % a temporary variable
    s_x = [x Y+A Y-A]; 
    
    % Unscented transform
    s_x_mean = zeros(6,1);
    for i = 1:size(s_x,2) % for all the sigma points
        s_x(:,i) = process(s_x(:,i),u(k-1),dt); % assume additive noise, so sigma points have no noise
        s_x_mean=s_x_mean+Wm(i)*s_x(:,i); 
    end
    e_x = s_x-s_x_mean(:,ones(1,13));
    P = e_x*diag(Wc)*e_x' + Q; 
    x_prior = [x_prior s_x_mean];

    % Measurement update.
    % Generate sigma points
    A = c*chol(P)';
    Y = x(:,ones(1,numel(x))); % a temporary variable
    s_x = [x Y+A Y-A]; 
    % Unscented transform
    s_z_mean = zeros(2,1);
    for i = 1:size(s_x,2) % for all the sigma points
        s_z(:,i) = measure(s_x(:,i)); % assume additive noise, so sigma points have no noise
        s_z_mean = s_z_mean + Wm(i) * s_z(:,i); 
    end
    e_z = s_z-s_z_mean(:,ones(1,13));
    P_zz = e_z*diag(Wc)*e_z' + R; 
    P_xz = e_x*diag(Wc)*e_z';
    
    K = P_xz * inv(P_zz);
    x = x + K*(y_obs(:,k)-s_z_mean);
    P = P - K*P_zz*K';
    x_posterior = [x_posterior x];
    
%     subplot(211)
%     scatter(t,x_posterior(1,end),'r.')
%     hold on
%     scatter(t,data(t,4),'k.')
%     scatter(t,x_prior(1,end),'b.')
%     xlim([0 900])
%     legend('posterior','measured','prior');
%     title('spacing')
% 
%     subplot(212)
%     scatter(t,x_posterior(2,end),'r.')
%     hold on
%     scatter(t,data(t,2),'k.')
%     scatter(t,x_prior(2,end),'b.')
%     xlim([0 900])
%     legend('posterior','measured','prior');
%     title('velocity')
%     
%     drawnow;
end

%% Plot results
% t = 0 : dtPlot : tf;
% 
% figure;
% plot(t, xArray(3,:) - xhatArray(3,:));
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Seconds'); ylabel('w_n^2 Estimation Error');
% 
% figure;
% plot(t, P3Array);
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Seconds'); ylabel('Variance of w_n^2 Estimation Error');
% 
% disp(['Final estimation error = ', num2str(xArray(3,end)-xhatArray(3,end))]);

figure;
subplot(211)
plot(data(:,1)',x_posterior(1,:))
hold on
plot(data(:,1)',data(:,4)')
plot(data(:,1)',x_prior(1,:))
hold off
legend('posterior','measured','prior');
title('spacing')

subplot(212)
plot(data(:,1)',x_posterior(2,:))
hold on
plot(data(:,1)',data(:,2)')
plot(data(:,1)',x_prior(2,:))
hold off
legend('posterior','measured','prior');
title('velocity')

figure;
subplot(311)
plot(data(:,1)',x_posterior(4,:))
title('a')

subplot(312)
plot(data(:,1)',x_posterior(5,:))
title('b')

subplot(313)
plot(data(:,1)',x_posterior(6,:))
title('hm')
