clear all
% close all
clc
tic

%% plot data
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

%% parameters setup
%Model parameters
model = params_OV;
model.Nv = 2; % number of vehicles
% Optimal Velocity Model parameters
% model.alpha = 1.1; % model sensitivity
% model.a = 12; % 10 affect the converging point (height)
% model.hm = 3; % positive height
% model.b = 4;
% https://link.springer.com/content/pdf/10.1 631%2Fjzus.A0900370.pdf pg524
model.m = 2; % measurement dimension
% Parameters in Particle Filter
param = params_PF;
param.Np = 1000; % number of particles
param.init_stdev = [3;2;2;6;6;4]; %10 5
param.model_stdev = [2;1.5;0.1;0.1;0.1;0.1]; % process error std dev, position and velocity. w_k
param.meas_stdev = [0.5;0.2]; % measurement error std dev. eta_k

param.dt = 0.1; % sample frequency 10 Hz
param.T = size(data,1)*param.dt; % total simulation time in sec
% param.eps = [0.1;0.3]; % uncertainties of initial x
param.init_spacing = 35; % m
param.init_vel = 24; % in m/s 
R = diag(param.meas_stdev.^2); % for weight calculation

timestep = 0:param.dt:param.T;
num_tstep = length(timestep);

% prior distributions
% a_dist = makedist('Uniform','lower',2,'upper',8);
% b_dist = makedist('Uniform','lower',12,'upper',25);
% hm_dist = makedist('Uniform','lower',9,'upper',18);
a_dist = makedist('Uniform','lower',5,'upper',50);
b_dist = makedist('Uniform','lower',5,'upper',50);
hm_dist = makedist('Uniform','lower',-10,'upper',20);
alpha_dist = makedist('Uniform','lower',0,'upper',15);
% a_stdev = 0.1;
% b_stdev = 0.1;
% hm_stdev = 0.1;
% alpha_stdev = 0.0;

directory = pwd;
foldername = 'parameter_both';
mkdir(foldername);
directory = fullfile(directory,foldername);

%% Initialize particles and weights
% x is the column state vector. x = [s;v;theta]
x = zeros(6,param.Np,num_tstep); % estimated states of the follower vehicle
wt = ones(1, param.Np)/param.Np; % particle weights
theta_0 = [4.3049 26.8558 31.0387 11.9599];
x_init = [data(1,4); data(1,2);theta_0']; % states of the second vehicle
for i = 1:param.Np 
    x(:,i,1) = x_init + normrnd(0, param.init_stdev);
end 
% x(x<0) = 0;
x(x(1:2,:,:)<0) = 0;
y = measure(x(:,:,1));

%% Measurements
y_obs = [data(:,4)';data(:,2)']; % s and v

%% initialize arrays
s_est = zeros(1,num_tstep - 1); p_max = s_est; p_min = s_est; % stores position estimation
v_est = s_est; v_max = s_est; v_min = s_est; % stores velocity estimation
s_est(1) = mean(x(1,:,1)); p_max(1) = max(x(1,:,1)); p_min(1) = min(x(1,:,1));
v_est(1) = mean(x(2,:,1)); v_max(1) = max(x(2,:,1)); v_min(1) = min(x(2,:,1));
x_n = x(:,:,1);
y_n = y;
M = zeros(2,1);
N_eff = zeros(1,num_tstep-1);

%% Perform time propagation to obtain priori particles
% x_axis = 5:0.1:25;

% a_prev = mean(a_dist);
% b_prev = mean(b_dist);
% hm_prev = mean(hm_dist);
% alpha_prev = mean(alpha_dist);

measured_state = 'both';

figure;
for k = 1:num_tstep-1
    x_axis = timestep(k).* ones(1,param.Np);
%     a_rand = random(a_dist,[param.Np 1]) + normrnd(0,a_stdev,[param.Np 1]);
%     b_rand = random(b_dist,[param.Np 1]) + normrnd(0,b_stdev,[param.Np 1]);
%     hm_rand = random(hm_dist,[param.Np 1]) + normrnd(0,hm_stdev,[param.Np 1]);
%     alpha_rand = random(alpha_dist,[param.Np 1]) + normrnd(0,alpha_stdev,[param.Np 1]);

%     a_rand = a_dist.mean + normrnd(0,a_stdev,[param.Np 1]);
%     b_rand = b_dist.mean + normrnd(0,b_stdev,[param.Np 1]);
%     hm_rand = hm_dist.mean + normrnd(0,hm_stdev,[param.Np 1]);
%     alpha_rand = alpha_dist.mean + normrnd(0,alpha_stdev,[param.Np 1]);

    for i = 1:param.Np
        x_n(:,i) = process(x(:,i,k),data(k,3),param.dt) + normrnd(0,param.model_stdev);
        y_n(:,i) = measure(x_n(:,i)) + normrnd(0,param.meas_stdev);
        
        switch measured_state
            case 'both'
                wt(i) = (2 * pi)^(-model.m/2) * (sqrt(sum(sum(abs(R).^2))))^(-1/2) * ...
                     exp(-1/2 * (y_obs(:,k) - y_n(:,i))'* R^(-1) * (y_obs(:,k) - y_n(:,i)));
            case 's'
                wt(i) = (2 * pi)^(-1/2) * (param.meas_stdev(1))^(-1/2) * ...
                    exp(-1/2 * (y_obs(1,k) - y_n(1,i))^2* (param.meas_stdev(1))^(-1));
            case 'v'
                wt(i) = (2 * pi)^(-1/2) * (param.meas_stdev(2))^(-1/2) * ...
                 exp(-1/2 * (y_obs(2,k) - y_n(2,i))^2* (param.meas_stdev(2))^(-1)); 
        end
    end
    % normalize weight
    wt = wt./sum(wt); % make sure all the weights sum up to one
    % effective particle size
    N_eff(k) = (sum(wt.^2))^(-1);
    
%     resampling
    for i = 1:param.Np
        index = find(rand <= cumsum(wt),1);
        x_n(:,i) = x_n(:,index);
%         a_rand(i) = a_rand(index);
%         b_rand(i) = b_rand(index);
%         hm_rand(i) = hm_rand(index);
    end
    
%     b_dist = fitdist(b_rand,'Normal'); %b_dist = truncate(b_dist,0,inf);
%     a_dist = fitdist(a_rand,'Normal'); %a_dist = truncate(a_dist,0,inf);
%     hm_dist = fitdist(hm_rand,'Normal'); %hm_dist = truncate(hm_dist,0,inf);
%     alpha_dist = fitdist(alpha_rand,'Normal'); %alpha_dist = truncate(alpha_dist,0,inf);
    
    % ************** plot parameters **************
%     subplot(411)
%     plot_prob(k,a_prev,a_dist,'r'); hold on
%     title('a');
%     set(gca,'Fontsize',15)
%     subplot(412)
%     plot_prob(k,b_prev,b_dist,'g'); hold on
%     title('b');
%     set(gca,'Fontsize',15)
%     subplot(413)
%     plot_prob(k,hm_prev,hm_dist,'b'); hold on
%     title('hm');
%     set(gca,'Fontsize',15)
%     subplot(414)
%     plot_prob(k,alpha_prev,alpha_dist,'b'); hold on
%     title('alpha');
%     set(gca,'Fontsize',15)
    
%     drawnow
% %     filename = sprintf('param_%03d',k);
% %     path = fullfile(directory,filename);
% %     saveas(gca,path,'png')
    % ************** plot parameters **************
    
    switch measured_state
        case 'both'
            s_est(k) = mean(x_n(1,:)); p_max(k) = max(x_n(1,:)); p_min(k) = min(x_n(1,:));
            v_est(k) = mean(x_n(2,:)); v_max(k) = max(x_n(2,:)); v_min(k) = min(x_n(2,:));
            alpha_est(k) = mean(x_n(3,:)); alpha_max(k) = max(x_n(3,:)); alpha_min(k) = min(x_n(3,:));
            a_est(k) = mean(x_n(4,:)); a_max(k) = max(x_n(4,:));a_min(k) = min(x_n(4,:));
            b_est(k) = mean(x_n(5,:)); b_max(k) = max(x_n(5,:)); b_min(k) = min(x_n(5,:));
            hm_est(k) = mean(x_n(6,:)); hm_max(k) = max(x_n(6,:)); hm_min(k) = min(x_n(6,:));

            %         case 's'
%             s_est(k) = mean(x_n(1,:)); p_max(k) = max(x_n(1,:)); p_min(k) = min(x_n(1,:));
%             v_est(k) = (s_est(k)-s_est(k-1))/param.dt;
%         case 'v'
%             v_est(k) = mean(x_n(2,:)); v_max(k) = max(x_n(2,:)); v_min(k) = min(x_n(2,:));
%             s_est(k) = s_est(k-1)+v_est(k-1)*param.dt;
    end
    
    x(:,:,k+1) = x_n;
%     a_prev = mean(a_dist);
%     b_prev = mean(b_dist);
%     hm_prev = mean(hm_dist);
%     alpha_prev = mean(alpha_dist);
    
    % ************** plot particles **************
%     if mod(k,10) == 0
%         subplot(211)
%         scatter(x_axis,x(1,:,k),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.01);
%         xlim([0 120]);
%         hold on
%         subplot(212)
%         scatter(x_axis,x(2,:,k),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.01);
%         xlim([0 120]);
%         ylim([20 26]);
%         hold on
%         drawnow
%     end
    % ************** plot particles **************

end

%% plot effective particle size
figure
plot(3:num_tstep-1,N_eff(3:end),'LineWidth',2,'color',[0.8,0.61,0]);
xlabel('Simulation time step');
ylabel('Effective Particle Size');
title('Effective Particle Size');
set(gca,'fontsize',22)

%% plot state 1 (spacing)
figure
plot(timestep(1:end-1), data(:,4)', '-b','lineWidth',2);hold on
plot(timestep(1:end-1), s_est, '-.r','lineWidth',2);
plot(timestep(1:end-1), p_max, '-.k');
plot(timestep(1:end-1), p_min, '-.k');
% for k = 1:num_tstep-1 % plot all the particles
%     x_axis_1 = timestep(k).* ones(1,param.Np);
%     scatter(x_axis_1,x1(:,k),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.04)
% end

set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('Spacing (m)');
legend('Measured state','Particle filter mean','max','min');

%% plot state 2 (velocity)
figure
plot(timestep(1:end-1), data(:,2)', '-b','lineWidth',3);hold on
plot(timestep(1:end-1), v_est, '-.r','lineWidth',3);
plot(timestep(1:end-1), v_max, '-.k');
plot(timestep(1:end-1), v_min, '-.k');

% for k = 1:num_tstep-6000 % plot all the particles
%     x_axis_2 = timestep(k).* ones(1,param.Np);
%     scatter(x_axis_2,x(2,:,k),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.1)
% end

set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Measured state','Particle filter mean','max','min');
toc

%% plot parameters
figure
subplot(411)
plot(timestep(1:end-1), alpha_est, '-.r','lineWidth',3); hold on;
plot(timestep(1:end-1), alpha_max, '-.k');
plot(timestep(1:end-1), alpha_min, '-.k');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('alpha');
legend('Particle filter mean','max','min');

subplot(412)
plot(timestep(1:end-1), a_est, '-.r','lineWidth',3); hold on;
plot(timestep(1:end-1), a_max, '-.k');
plot(timestep(1:end-1), a_min, '-.k');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('a');
legend('Particle filter mean','max','min');

subplot(413)
plot(timestep(1:end-1), b_est, '-.r','lineWidth',3); hold on;
plot(timestep(1:end-1), b_max, '-.k');
plot(timestep(1:end-1), b_min, '-.k');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('b');
legend('Particle filter mean','max','min');

subplot(414)
plot(timestep(1:end-1), hm_est, '-.r','lineWidth',3); hold on;
plot(timestep(1:end-1), hm_max, '-.k');
plot(timestep(1:end-1), hm_min, '-.k');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('hm');
legend('Particle filter mean','max','min');
%% compute RMSE
RMSE_p = sqrt(mean((x_true(1,1:end-1)-s_est(1:end-1)).^2))
RMSE_v = sqrt(mean((x_true(2,1:end-1)-v_est(1:end-1)).^2))
% normalize
% RMSE_p = sqrt(mean((x_true(1,1:end-1)-p_est(1:end-1)).^2))/mean(x_true(1,1:end-1))
% RMSE_v = sqrt(mean((x_true(2,1:end-1)-v_est(1:end-1)).^2))/mean(x_true(2,1:end-1))

