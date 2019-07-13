clear all
close all
clc
%% parameters setup
%Model parameters
model = params_OV;
model.Nv = 2; % number of vehicles. assume only two cars for now
% Optimal Velocity Model parameters
model.alpha = 1; % model sensitivity
% model.a = 12; % 10 affect the converging point (height)
% model.hm = 18; % positive height
% model.b = 12; %4
% https://link.springer.com/content/pdf/10.1 631%2Fjzus.A0900370.pdf pg524
model.alpha_true = 1; % model sensitivity
model.a_true = 8.97; % ±0.12 affect the converging point
model.hm_true = 12.78; % ±0.17
model.b_true = 20.01; % ±0.36
model.m = 2; % 2 states
% Parameters in Particle Filter
param = params_PF;
param.Np = 1000; % number of particles
param.init_stdev = [3;2]; %10 5
param.model_stdev = [0.5;0.5]; % [3;2]process error std dev, position and velocity. w_k
param.meas_stdev = [0.5;0.2]; %[1;0.5] measurement error std dev. eta_k
param.lead_stdev = 0;
param.T = 30; % total simulation time in sec
param.dt = 0.1; % time interval
param.eps = [0.1;0.3]; % uncertainties of initial x
param.init_spacing = 15; % m
param.init_vel = 10; % in m/s 10
R = eye(length(param.meas_stdev)).*(param.meas_stdev * param.meas_stdev'); % for weight calculation
timestep = 0:param.dt:param.T;
num_tstep = length(timestep);
% prior distributions
a_dist = makedist('Uniform','lower',2,'upper',8);
b_dist = makedist('Uniform','lower',12,'upper',25);
hm_dist = makedist('Uniform','lower',9,'upper',18);

directory = pwd;
foldername = 'parameter_both';
mkdir(foldername);
directory = fullfile(directory,foldername);
%% Initialize leading vehicle starting positions and velocities
vel_vals = zeros(model.Nv, num_tstep-1);
pos_vals = zeros(model.Nv, num_tstep-1);

for i = model.Nv:-1:1
    spacing = param.init_spacing + param.eps(1) * (rand*2-1);
    vel_vals(i,1) = param.init_vel; % first row first vehicle, second row second vehicle
    pos_vals(end+1-i,1) = spacing * i;
end

% duplicate for a biased leading vehicle profile that to be passed to the
% estimator
vel_bias = vel_vals;
pos_bias = pos_vals;

%% leading vehicle velocity and position profiles
for i = 0 : num_tstep - 2
    vel_vals(1, i + 1) = vel_fcn(i * param.dt, param.init_vel);
    vel_bias(1, i + 1) = vel_fcn_bias(i * param.dt);
end
pos_vals(1,:) = Euler(pos_vals(1,:), vel_vals(1,:), param.dt);
pos_bias(1,:) = Euler(pos_bias(1,:), vel_bias(1,:), param.dt);

%% Initialize particles and weights
x = zeros(2,param.Np); % estimated states of the follower vehicle
wt = ones(1, param.Np)/param.Np; % particle weights
x_init = [pos_vals(2,1); vel_vals(2,1)]; % states of the second vehicle
for i = 1:param.Np 
    x(:,i) = x_init + normrnd(0, param.init_stdev);
end 
x(x<0) = 0;
y = measure(x);

%% True states
x_true = zeros(2,num_tstep-1); % each column is a vector of state variables, 2 variables
x_true(:,1) = x_init;
for k = 2:num_tstep-1
    x_true(:,k) = process_true(x_true(:,k-1),pos_vals(1,k-1),model,param);
end
y_true = measure_true(x_true);

%% initialize arrays
p_est = zeros(1,num_tstep - 1); p_max = p_est; p_min = p_est; % stores position estimation
v_est = p_est; v_max = p_est; v_min = p_est; % stores velocity estimation
p_est(1) = mean(x(1,:)); p_max(1) = max(x(1,:)); p_min(1) = min(x(1,:));
v_est(1) = mean(x(2,:)); v_max(1) = max(x(2,:)); v_min(1) = min(x(2,:));
x_n = x;
y_n = y;
M = zeros(2,1);
N_eff = zeros(1,num_tstep-1);

%% Perform time propagation to obtain priori particles
x_axis = 5:0.1:25;
a_prev = mean(a_dist);
b_prev = mean(b_dist);
hm_prev = mean(hm_dist);

measured_state = 'both';
for k = 2:num_tstep-1
    a_rand = random(a_dist,[param.Np 1]);
    b_rand = random(b_dist,[param.Np 1]);
    hm_rand = random(hm_dist,[param.Np 1]);
    for i = 1:param.Np
        x_n(:,i) = process(x(:,i),pos_bias(1,k-1),model,param,b_rand(i),a_rand(i),hm_rand(i)) + normrnd(0,param.model_stdev);
        y_n(:,i) = measure(x_n(:,i)) + normrnd(0,param.meas_stdev);
        switch measured_state
            case 'both'
                wt(i) = (2 * pi)^(-model.m/2) * (sqrt(sum(sum(abs(R).^2))))^(-1/2) * ...
                     exp(-1/2 * (y_true(:,k) - y_n(:,i))'* R^(-1) * (y_true(:,k) - y_n(:,i)));
            case 'p'
                wt(i) = (2 * pi)^(-1/2) * (param.meas_stdev(1))^(-1/2) * ...
                    exp(-1/2 * (y_true(1,k) - y_n(1,i))^2* (param.meas_stdev(1))^(-1));
            case 'v'
                wt(i) = (2 * pi)^(-1/2) * (param.meas_stdev(2))^(-1/2) * ...
                 exp(-1/2 * (y_true(2,k) - y_n(2,i))^2* (param.meas_stdev(2))^(-1)); 
        end
    end
    % normalize weight
    wt = wt./sum(wt); % make sure all the weights sum up to one
    % effective particle size
    % HOW TO AVOID THE CURSE OF DIMENSIONALITY: SCALABILITY OF PARTICLE FILTERS WITH AND WITHOUT IMPORTANCE WEIGHTS
    N_eff(k) = (sum(wt.^2))^(-1);
    
%     resampling
    for i = 1:param.Np
        x_n(:,i) = x_n(:,(find(rand <= cumsum(wt),1)));
        a_rand(i) = a_rand(find(rand <= cumsum(wt),1));
        b_rand(i) = b_rand(find(rand <= cumsum(wt),1));
        hm_rand(i) = hm_rand(find(rand <= cumsum(wt),1));
    end
    
    b_dist = fitdist(b_rand,'Normal'); b_dist = truncate(b_dist,0,inf);
    a_dist = fitdist(a_rand,'Normal'); a_dist = truncate(a_dist,0,inf);
    hm_dist = fitdist(hm_rand,'Normal'); hm_dist = truncate(hm_dist,0,inf);
    
    % ************** plot parameters **************
    subplot(311)
    plot_prob(k,a_prev,a_dist,'r'); hold on
    title(sprintf('a: true value = %.02f',model.a_true));
    set(gca,'Fontsize',15)
    subplot(312)
    plot_prob(k,b_prev,b_dist,'g'); hold on
    title(sprintf('b: true value = %.02f',model.b_true));
    set(gca,'Fontsize',15)
    subplot(313)
    plot_prob(k,hm_prev,hm_dist,'b'); hold on
    title(sprintf('h_m: true value = %.02f',model.hm_true));
    set(gca,'Fontsize',15)
    
    drawnow
    filename = sprintf('param_%03d',k);
    path = fullfile(directory,filename);
    saveas(gca,path,'png')
    % ************** plot parameters **************
    switch measured_state
        case 'both'
            p_est(k) = mean(x_n(1,:)); p_max(k) = max(x_n(1,:)); p_min(k) = min(x_n(1,:));
            v_est(k) = mean(x_n(2,:)); v_max(k) = max(x_n(2,:)); v_min(k) = min(x_n(2,:));
        case 'p'
            p_est(k) = mean(x_n(1,:)); p_max(k) = max(x_n(1,:)); p_min(k) = min(x_n(1,:));
            v_est(k) = (p_est(k)-p_est(k-1))/param.dt;
        case 'v'
            v_est(k) = mean(x_n(2,:)); v_max(k) = max(x_n(2,:)); v_min(k) = min(x_n(2,:));
            p_est(k) = p_est(k-1)+v_est(k-1)*param.dt;
    end
    
    x = x_n;
    a_prev = mean(a_dist);
    b_prev = mean(b_dist);
    hm_prev = mean(hm_dist);

end

%% plot effective particle size
figure
plot(3:num_tstep-1,N_eff(3:end),'LineWidth',2,'color',[0.8,0.61,0]);
xlabel('Simulation time step');
ylabel('Effective Particle Size');
% title('Effective Particle Size');
set(gca,'fontsize',22)

%% plot state 1 (position)
figure
plot(timestep(1:end-1), x_true(1,:), '-b','lineWidth',2);hold on
plot(timestep(1:end-1), p_est, '-.r','lineWidth',2);
plot(timestep(1:end-1), p_max, '-.k');
plot(timestep(1:end-1), p_min, '-.k');
% for k = 1:num_tstep-1 % plot all the particles
%     x_axis_1 = timestep(k).* ones(1,param.Np);
%     scatter(x_axis_1,x1(:,k),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.04)
% end

set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('Position (m)');
legend('True state','Particle filter mean','max','min');

%% plot state 2 (velocity)
figure
plot(timestep(1:end-1), x_true(2,:), '-b','lineWidth',3);hold on
plot(timestep(1:end-1), v_est, '-.r','lineWidth',3);
plot(timestep(1:end-1), v_max, '-.k');
plot(timestep(1:end-1), v_min, '-.k');

% for k = 1:num_tstep-1 % plot all the particles
%     x_axis_2 = timestep(k).* ones(1,param.Np);
%     scatter(x_axis_2,x2(:,k),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.1)
% end

set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('True state','Particle filter mean','max','min');

%% compute RMSE
RMSE_p = sqrt(mean((x_true(1,1:end-1)-p_est(1:end-1)).^2))
RMSE_v = sqrt(mean((x_true(2,1:end-1)-v_est(1:end-1)).^2))
% normalize
% RMSE_p = sqrt(mean((x_true(1,1:end-1)-p_est(1:end-1)).^2))/mean(x_true(1,1:end-1))
% RMSE_v = sqrt(mean((x_true(2,1:end-1)-v_est(1:end-1)).^2))/mean(x_true(2,1:end-1))

