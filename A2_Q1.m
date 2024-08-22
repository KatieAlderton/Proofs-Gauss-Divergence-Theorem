%{
Name: Katie Alderton
Date: 19/04/2024
Description: Visualisation 

Input(s):
- Number of spatial divisions (n_div)
- Number of faces (n_face)

Output(s): 
- A figure showing a graph of the unit sphere, its normal vectors and the vector field
it is located within

%}
%
clear 
close all
clc
%
% Set up 3D grid
n_div=20; % number of equal spatial divisions
x=linspace(-1,1,n_div);
y=linspace(-1,1,n_div);
z=linspace(-1,1,n_div);
[X,Y,Z]=meshgrid(x,y,z); % Sets up a 3D grid based on x,y and z
%
% Define the vector field F as denoted in the Assignment brief
F_X=X.*(Z.^2);
F_Y=X.*Y;
F_Z=Y.*Z;
%
% Generate coordinates on the spherical surface
n_face=30; % Number of faces on the sphere (this will be 30x30 as per the question)
[s_x,s_y,s_z]=sphere(n_face); % Generates coordinates on the spherical surface
[n_x,n_y,n_z]=surfnorm(s_x,s_y,s_z); % Generates surface normals
%
% Visualisation of the unit sphere
figure(1)
hold on
quiver3(X,Y,Z,F_X,F_Y,F_Z,'Color','b') % Plots the vector field F
surf(s_x,s_y,s_z) % Plots the spherical Gaussian surface 
quiver3(s_x,s_y,s_z,n_x,n_y,n_z,'Color','r') % Plots the normals to the spherical surface
view(-37.5,30) % Sets a 3D view point.
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('Visualisation of a the unit sphere in vector field F')

