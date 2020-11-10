function cord = pendout(yout)
%   function for easy definition of x and y coordinates from ride
    theta1 = yout(:,1);
    theta2 = yout(:,2);
    global l1 l2
    x1 = l1*sin(theta1);
    y1 = -l1*cos(theta1);
    x2 = x1 + l2*sin(theta2);
    y2 = y1 - l2*cos(theta2);

    cord = [x1, y1, x2, y2];
end