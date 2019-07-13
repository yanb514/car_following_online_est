% Extended Kalman filter for parameter estimation.
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

figure;
plot(data(:,1),data(:,2)); % follower speed m/s
hold on
plot(data(:,1),data(:,3)); % leader speed m/s
% plot(data(:,1),data(:,4)); % spacing m
legend('Follower','Leader');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
set(gca,'fontsize',20)

dt = 0.1;
T = data(end,1);
timestep = dt:dt:T;
num_tstep = length(timestep);

%% Matrices initialize
w = [2 1.5 0.1 0.1 0.1 0.1]; % process noise stdev
Q = diag(w.^2); % covariance of process noise
v = [0.5 0.2]; % measurement noise stdev
R = diag(v.^2); % covariance of measurement noise
H = [1 0 0 0 0 0; 0 1 0 0 0 0]; % measurement equation
L = eye(6); % covariance of process error
M = eye(2); % covariance of measurement error

theta_0 = [4.3049 26.8558 31.0387 11.9599];
x = [data(1,4); data(1,2); theta_0']; % initial state
% x_0 = 1 * x; % initial state estimate
P = diag([3 2 2 6 6 4].^2); % initial error

% Initialize arrays for later plotting
x_prior = x;
x_posterior = x;

u = data(:,3); % input is the lead vehicle's speed
y_obs = [data(:,4)';data(:,2)']; % s and v

%%
figure;
for t = 2:num_tstep+1
    % Prior update.
    A = A_calc(x);
    x = process(x,u(t-1),dt) + w'*randn;
    P = A*P*A' + L*Q*L';
    x_prior = [x_prior x];
    
    % Measurement update.
%     if mod(t,10)==0
        K = P*H'* inv(H*P*H' + M*R*M');
        x = x + K*(y_obs(:,t)+v'*randn-measure(x));
        P = (eye(size(P,1))-K*H)*P;
        x_posterior = [x_posterior x];
%     end
    
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
xlabel('time (sec)')
legend('posterior','measured','prior');
title('spacing (m)')

subplot(212)
plot(data(:,1)',x_posterior(2,:))
hold on
plot(data(:,1)',data(:,2)')
plot(data(:,1)',x_prior(2,:))
hold off
xlabel('time (sec)')
legend('posterior','measured','prior');
title('velocity')

figure;
subplot(411)
plot(data(:,1)',x_posterior(3,:))
title('alpha')

subplot(412)
plot(data(:,1)',x_posterior(4,:))
title('a')

subplot(413)
plot(data(:,1)',x_posterior(5,:))
title('b')

subplot(414)
plot(data(:,1)',x_posterior(6,:))
title('hm')

%%
% x_test(:,1) = [data(1,4); data(1,2); theta_0'];
% for t = 1: num_tstep
%     x_test(:,t+1) = process(x_posterior(:,t),data(t,3),0.1);
% end
