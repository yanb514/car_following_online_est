function A = A_calc(x)
    alpha = x(3);
    a = x(4);
    b = x(5);
    hm = x(6);
    s = x(1);
    derivative = a/b*(1-(tanh((s-hm)/b))^2);
    
    A = [0 -1 0 0 0 0;
        alpha*derivative -alpha 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
end

