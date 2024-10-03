% Some parameters
N = 100;
M = 100;
L = 1;
H = 1;
c = [L/2, H/2];
r = 0.1;

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Channel flow');
xlim([0,L]);
ylim([0,H]);
axis equal