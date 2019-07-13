gap = 0:70;
ov_model = zeros(1,length(gap));
ov_true = zeros(1,length(gap));

t = 300;
for i=1:length(gap)
%     ov_model(i) = ov_calc(a_dist.mean,b_dist.mean,hm_dist.mean,gap(i));
    ov_model(i) = ov_calc(xArray(4,t),xArray(5,t),xArray(6,t),gap(i));
end

figure
plot(gap,ov_model,'r');
xlabel('spacing (m)');
ylabel('Optimal Velocity (m/s)');
% hold on;
% plot(gap,ov_true,'b');
% legend('model','true');
hold on
scatter(data(:,4),data(:,2),'.');
% scatter(s_est,v_est,'r.')
scatter(xArray(1,:),xArray(2,:),'r.')

title(sprintf('V(s) alpha=%.2f,a=%.2f b=%.2f hm=%.2f, v_m = %.2f',alpha_dist.mean,a_dist.mean,b_dist.mean,hm_dist.mean,a_dist.mean*(1+tanh(hm_dist.mean/b_dist.mean))))
set(gca,'fontsize',20,'fontname','Times')

% lambda = lambda_calc(model.alpha,a_dist.mean,b_dist.mean,hm_dist.mean,param);
% if lambda < 0
%     disp('string stable')
% else
%     disp('string unstable')
% end

%% plot measured vs. simulated
% OV = ov_calc(a_dist.mean,b_dist.mean,hm_dist.mean,data(:,4)');
% accr = alpha_dist.mean * (OV-data(:,2)');
% v_sim = Euler(data(1,2),accr,param.dt);
% s_dot = data(:,3)'-data(:,2)';
% s_sim = Euler(data(1,4),s_dot,param.dt);

% x_sim(:,1) = [data(1,4);data(1,2)];
% v_vals(1) = data(1,2);
% s_vals(1) = data(1,4);

% theta = [alpha_dist.mean,a_dist.mean,b_dist.mean,hm_dist.mean];
theta = [xArray(3,t),xArray(4,t),xArray(5,t),xArray(6,t)];
x_sim(:,1) = [data(1,4);data(1,2);theta'];
% theta = [10.84,33.7,25,-4.88];
for t=1:num_tstep
%     v = v_vals(t-1);
%     s = s_vals(t-1);
%     v_l = data(t-1,3);
%     a = alpha_dist.mean*(ov_calc(a_dist.mean,b_dist.mean,hm_dist.mean,s) - v);
%     v_vals(t) = v + a*param.dt;
%     s_vals(t) = s + (v_l-v)*param.dt;
    x_sim(:,t+1) = process(x_sim(:,t),data(t,3),0.1);
end
figure;
subplot(211)
plot(data(:,1)',x_sim(1,:))
hold on
plot(data(:,1)',data(:,4)')
hold off
legend('simulated','measured');
title('spacing')

subplot(212)
plot(data(:,1)',x_sim(2,:))
hold on
plot(data(:,1)',data(:,2)')
hold off
legend('simulated','measured');
title('velocity')

%%
subplot(211)
plot(data(:,1)',xArray(1,:))
hold on
plot(data(:,1)',data(:,4)')
hold off
legend('simulated','measured');
title('spacing')

subplot(212)
plot(data(:,1)',xArray(2,:))
hold on
plot(data(:,1)',data(:,2)')
hold off
legend('simulated','measured');
title('velocity')
% plot(data(:,1)',data(:,4)')


%% to ballpark the model stddev
error_arr = zeros(2,10000);
for i = 1:10000
    p_input = 17.73+rand*(599-17.73);
    v_input = 6.014+rand*(19.46-6.014);
    time = randi([1 600]);
    error_arr(:,i) = process([p_input;v_input],pos_bias(1,time),model,param)-process_true([p_input;v_input],pos_vals(1,time),model,param);
end

%%
meas_output = zeros(2,600);
for i = 1:600
    meas_output(:,i) = measure_true(x_true(:,i));
end
RMSE_p = sqrt(mean((x_true(1,1:end-1)-meas_output(1,1:end-1)).^2))
RMSE_v = sqrt(mean((x_true(2,1:end-1)-meas_output(2,1:end-1)).^2))
    
%%
p1 = (s_x-x)*(s_x-x)'./(2*size(x,1));

p2 = 0;
for i = 1:12
    p2 = p2 + ans(:,i)*ans(:,i)';
end


p2 = p2/12;