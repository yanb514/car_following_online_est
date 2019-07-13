% function [x_output] = process(x,s,alpha,b_rand,a_rand,hm_rand,param)
% % next move after passing through the process equation.
% % input is a vector of state variables, x is a column vector
% % model.b = random(b_dist);
% % dynamic equation
%     x1_dot = x(2);
%     if x1_dot < 0
%         disp('negative velocity');
%     end
%     x2_dot = (alpha) * (ov_calc(a_rand,b_rand,hm_rand,s)-x(2));
%     x_1 = x(1) + param.dt * x1_dot; % first state: position
%     x_2 = x(2) + param.dt * x2_dot; % second state: velocity
%     x_output = [x_1; x_2];
% end

function [x_next] = process(x,input,dt)
% next move after passing through the process equation.
% input is a vector of state variables, x is a column vector
% model.b = random(b_dist);
% dynamic equation
    alpha = x(3);
    a = x(4);
    b = x(5);
    hm = x(6);
    
    s = x(1);
    v = x(2);
    v_l = input;
    s_dot = v_l - v;
    v_dot = alpha * (ov_calc(a,b,hm,s) - v);
    
    if v < 0
        disp('negative velocity');
    end

    s_next = s + dt * s_dot; % first state: position
    v_next = v + dt * v_dot; % second state: velocity
    x_next = [s_next; v_next; x(3:6)];
end