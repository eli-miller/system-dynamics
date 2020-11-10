function phi = pend(z, zdot, t)
%EOM of double pendelum
global m1 m2 l1 l2 g
theta1 = z(1);
theta2 = z(2);
f1 = z(3);
f2 = z(4);

thetadot1 = zdot(1);
thetadot2 = zdot(2);
f1dot = zdot(3);
f2dot = zdot(4);

phi = zeros(4,1);
%Inertia Matrix
M = [1, 0, 0, 0;
     0, 1, 0, 0;
     0, 0, l1^2*(m1+m2), m1*l1*l2*cos(theta1-theta2);
     0, 0, l1*cos(theta1-theta2), l2];
 
V = [-f1;
     -f2;
     m2*l1*l2*f2.^2*sin(theta1-theta2) + g*l1*(m1+m2)*sin(theta1);
     -l1*f1.^2*sin(theta1-theta2) + g*sin(theta2)];
 
 phi = M*zdot + V;

end

