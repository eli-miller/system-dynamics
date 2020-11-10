clear all; close all; clc;

%Variables to be used in all associated functions

global m1 m2 l1 l2 g

m1 = 2;
m2 = 3;
l1 = 1;  
l2 = 2;
g = 9.81;

%plotting and animation options
bnd = l1+l2;
titlevec = ["Part a","Part b","Part c","Part d"];
time = .2;
N = 1001;
tstop = 15;

tspan = linspace(0,tstop,N);

zcond = [pi/8, pi/8, 0, 0;
        0, pi/2, 0, 0;
        pi/2, 0, 0, 0;
        3*pi/4, pi/2, 0, 0];
 
zdotcond = [zcond(:,3:end),zeros(4,2)];
for iter = 1:4
    %for each case, compute the flow derivitives
    zd0 = pendinit(zcond(iter,:))';
    part = iter;
    z0 = zcond(part,:)';
    
    %Solve DAE
    
    [tout,yout,INFO] = ride('pend', '', tspan, z0, zd0);
    cord = pendout(yout);
    
    ax = figure(1)
    
    %add a subplot for every case, and label it as such
    subplot(2,2,iter)
    plot(cord(:,1),cord(:,2),'b',cord(:,3),cord(:,4),'r')
    axis([-bnd, bnd, -bnd, bnd])
    axis equal
    title(titlevec(iter))
 end
%Save file as png
% saveas(ax,'Problem6.png')
%%
%Some animation plotting practice!
%Not part of homework submission!



figure(2)
    for i = 1:length(tout)
        plot([0,cord(i,1)],[0,cord(i,2)],'k',[cord(i,1),cord(i,3)],[cord(i,2),cord(i,4)],'k',cord(i,1),cord(i,2),'ro',cord(1:i,3),cord(1:i,4),'r.-','LineWidth',2)
        axis([-bnd, bnd, -bnd, bnd])
        pause(time/N)
    end