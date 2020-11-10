clear all; close all; clc;

global m k b F

m = 0.3;
k = 45;
b = 0.75;
F = 0.5;

z0 = [0, 0, 0, 0]';
zdot0 = [0, 0, 0, F/m]';
zdot0 = zeros(size(z0));
tspan = linspace(0, 5,501)';

[t,yout,INFO] = ride('eom', '', tspan, z0, zdot0);


legvec = ["x_{1}", "x_{2}", "xdot_{1}", "xdot_{1}"];

fig = figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(t, yout(:,i))
    legend(legvec(i))
end
saveas(fig,'Problem5.png')

