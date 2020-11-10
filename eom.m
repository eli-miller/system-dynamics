function phi = eom(z, zdot, t)
   
    
global m k b F

    x1 = z(1);
    x2 = z(2);
    f1 = z(3);
    f2 = z(4);

    x1dot = zdot(1);
    x2dot = zdot(2);
    f1dot = zdot(3);
    f2dot = zdot(4);

    phi = zeros(4<1);

    M = [1, 0, 0, 0;
        0, 1, 0, 0;
        0, 0, 0, 0;
        0, 0, 0, m];
    
    V = [-f1;
        -f2;
        k*x1 - b*(f2-f1);
        b*(f2-f1) - F];

    phi = M*zdot + V;
    
end
    