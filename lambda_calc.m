% reference: http://www1.maths.leeds.ac.uk/~jaward/publications/JAW-UTSG-2010.pdf
%function [output] = lambda_calc(a_val,b_val,hm_val,params)
function [output] = lambda_calc(alpha_model,a,b,hm,params)
    clear a s alpha v dv h_m b 
    syms a s h_m b v dv
    ov_fcn = alpha_model*(a*(tanh((s-h_m)/b)+tanh(h_m/b))-v);
    fs = diff(ov_fcn,s);
    fdv = diff(ov_fcn,dv);
    fv = diff(ov_fcn,v);
    lambda(h_m,b,a) = fs/(fv^3)*(fv^2/2-fdv*fv-fs);
    lambda_2 = subs(lambda,{a,s,h_m,b},{a,params.init_spacing,hm,b});
    lambda_2 = double(lambda_2);
    output = lambda_2;
end