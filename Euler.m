% Euler method (first order integration)
% function [output] = Euler(x,v,sim_t,params,v_lead)
%     x_next = double.empty(0,length(x));
%     v_next = double.empty(0,length(v));
%     derivatives = calcf(x,v,sim_t,params,v_lead);
%     dx = derivatives(1,:);
%     dv = derivatives(2,:);
%     for i = 1:params.N
%         x_next(i) = x(i) + dx(i) * params.dt;
%         v_next(i) = v(i) + dv(i) * params.dt;
%     end
%     output = [x_next;v_next];
function [array] = Euler(initial,grad_array,time_step)
array = zeros(size(grad_array));
array(1) = initial;
    for i = 1:length(array)-1   
        array(i+1) = array(i) + grad_array(i) * time_step;
    end
end
