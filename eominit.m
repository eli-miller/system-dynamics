function zd0 = eominit(z)
    global m k b F

    x1 = z(1);
    x2 = z(2);
    f1 = z(3);
    f2 = z(4);
    
    zd0 = zeros(size(z));
    zd0(1:2) = z(3:4);
    
    
    M = [0, 0; 0, m];
    V =  [k*x1 - b*(f2-f1);
        b*(f2-f1) - F];
    zd0(3,4) = M\V;
    
end