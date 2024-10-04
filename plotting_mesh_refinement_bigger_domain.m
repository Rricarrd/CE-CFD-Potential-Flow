
% Some parameters
N = 150;
M = N
L = 3;
H = 3;
c = [L/2, H/2];
r = 0.1;

% Output data
data = readtable('output_mesh_refining_bigger_domain\150_output.csv');
X = table2array(data(:,1));
Y = table2array(data(:,2));
U = table2array(data(:,3));
V = table2array(data(:,4));
S = table2array(data(:,5));
rho = table2array(data(:,6));
Cp = table2array(data(:,7));
solid = table2array(data(:,8));

% Analytic output data
data_a = readtable('output_mesh_refining_bigger_domain/150_analytic_output.csv');
X_a = table2array(data_a(:,1));
Y_a = table2array(data_a(:,2));
S_a = table2array(data_a(:,5));
rho_a = table2array(data_a(:,6));
Cp_a = table2array(data_a(:,7));
solid_a = table2array(data_a(:,8));


% CONTOUR PLOT
figure(2)
%Contour plotting
[x_grid,y_grid] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
p_grid = griddata(X, Y, S ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contour(x_grid,y_grid,p_grid,'LevelList',linspace(min(S),max(S)))

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Streamline value (m^2/s)';

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Numerical stream lines ($\psi$)','Interpreter','latex');
grid on
axis equal
colormap cool


% ANALYTIC CONTOUR PLOT
figure(3)
%Contour plotting
[x_grid_a,y_grid_a] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
p_grid_a = griddata(X_a, Y_a, S_a ,x_grid_a,y_grid_a); %interpolates surface from  mesh and streamline values (cubic interpolation)
contour(x_grid_a,y_grid_a,p_grid_a,'LevelList',linspace(min(S),max(S)))

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure Coefficient (Cp)';

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Analytic stream lines ($\psi$)','Interpreter','latex');
grid on
axis equal
colormap cool



%% ERROR PLOT
figure(4);
% Error streamlines
data_e = readtable('output_mesh_refining_bigger_domain/150_error_output.csv');
[x_grid_e,y_grid_e] = meshgrid(linspace(0,L,M-2),linspace(0,H,N-2)); 
S_e = abs(table2array(data_e(:,5)));
S_e_mat = reshape(S_e,N-2,M-2);
pcolor(x_grid_e,y_grid_e,S_e_mat/max(S)*100)
xlabel('Y-axis [m]');
ylabel('Y-axis [m]');
title('Relative error on the domain');
colormap cool


rel_error_average = mean2(S_e_mat/max(S)*100)

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Relative error (%)';
colormap cool

% Alotting along y
figure(5);
plot(Y(1:148),S_e(1:148)/max(S)*100)
%Plot parameters
xlabel('Y-axis [m]');
ylabel('Relative error (%)');
title('Relative error along Y axis (X = 0.1)');
colormap cool

% %% RESULTS CONVERGENCE
% L = 0;
% H = 0;
% R = 0;
% N = 0;
% M = 0;
% Cl = 0;
% Cd = 0;
% Circ = 0;
% 
% for n=[15,25,50,75,100,150]
%    name = 'output_mesh_refining_bigger_domain\';
%    elems = sprintf('%i_results.csv',n);
%    data = table2array(readtable(append(name,elems)));
%    L = [L, str2double(data{1})];
%    H = [H, str2double(data{2})];
%    R = [R, str2double(data{3})];
%    N = [N, str2double(data{4})];
%    M = [M, str2double(data{5})];
%    Cl = [Cl, str2double(data{6})];
%    Cd = [Cd, str2double(data{7})];
%    Circ = [Circ, str2double(data{8})];
%     
% end
%    data = table2array(readtable('output_mesh_refining_bigger_domain\300_results.csv'));
%    L = [L, data(1)];
%    H = [H, data(2)];
%    R = [R, data(3)];
%    N = [N, data(4)];
%    M = [M, data(5)];
%    Cl = [Cl, data(6)];
%    Cd = [Cd, data(7)];
%    Circ = [Circ, data(8)];
% 
% plot(N, Cl)
