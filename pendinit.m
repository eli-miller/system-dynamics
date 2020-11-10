function zinit = pendinit(z0)
%Function to initialize the system for zdot0 
%Solves using linearized euqation that arizes at t0
global m1 m2 l1 l2 g


theta1 = z0(1);
theta2 = z0(2);
f1 = z0(3);
f2 = z0(4);

zdot = zeros(size(z0));
zdot(1) = z0(3);
zdot(2) = z0(4);

M =[
    l1^2*(m1+m2), m1*l1*l2*cos(theta1-theta2);
    l1*cos(theta1-theta2), l2];
 
V = [
     m2*l1*l2*f2.^2*sin(theta1-theta2) + g*l1*(m1+m2)*sin(theta1);
     -l1*f1.^2*sin(theta1-theta2) + g*sin(theta2)];

zdot(3:4) = M\V;

zinit = zdot;
end