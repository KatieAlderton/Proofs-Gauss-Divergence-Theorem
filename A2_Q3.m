%{
Name: Katie Alderton
Date: 19/04/2024
Description: Numerical Proof

Input(s):
-LHS:
   -n_div
-RHS:
   -n_faces
   - r=1 (radius [m])
Output(s): 
- A display of the approximation of the volume and surface integrals
%}
%
clear
clc
close all
%
% Set up 3D grid
n_div=100; % Number of divisions
x=linspace(-1,1,n_div);
y=linspace(-1,1,n_div);
z=linspace(-1,1,n_div);
[X,Y,Z]=meshgrid(x,y,z); 
%
% Determines vector field, F, components at each point in 3D grid
F_X=X.*(Z.^2);
F_Y=X.*Y;
F_Z=Y.*Z;
F=[F_X F_Y F_Z]; % Defines the fector field F in terms of the 3D grid
%
% Determines numerical divergence of F at each point in 3D grid
div_F=divergence(X,Y,Z,F_X,F_Y,F_Z);
% Determines the volume of a volume element,delta_V,corresponding to the 3D
% grid
Dx=x(2)-x(1); % [m] x length of a volume element
Dy=y(2)-y(1); % [m] y length of a volume element
Dz=z(2)-z(1); % [m] z length of a volume element
delta_V=Dx*Dy*Dz; % [m^3] volume of a volume element
%
% Unwrap relevant matrices prior to entering for loop
X=X(:);
Y=Y(:);
Z=Z(:);
div_F=div_F(:); % unwraps divergence matrice 
%
% Initialise numerical approximation of volume integral
lhsvolume=0;
%
% Usage of a for loop to sum the values of (divF*delta_V) for all points
% within the unit sphere (x^2+y^2+z^2=1). note r (radius)<1
for i=1:length(X)
    pos_curr=[X(i) Y(i) Z(i)]; % extracts the coordinates of the current position
    r_curr=norm(pos_curr); % Calculates the distance of the current position from the origin (0,0,0)
    if r_curr<1 % that is: if the distance is within the unit sphere (r<1)
        lhsvolume=lhsvolume+(div_F(i)*delta_V); % adds volume element to running numerical approximation
    else
    end
end


%% RHS Surface Integral
% A numerical approximation of the surface integral is as follows:
%
% use sphere function to generate points that define surface of a unit
% sphere (30x30 faces) and also the surfnormal function to generate normal
% vectors to the sphere
r=1; % [m] radius
n_faces=30; % number of faces
[sx,sy,sz]=sphere(n_faces); % generates coordinates on a spherical surface
[n_x,n_y,n_z]=surfnorm(sx,sy,sz); % generates the surface normals
%
% Determines the components of the vector field F on the surface of the
% unit sphere
S_Fx=sx.*sz.^2; % x component of F on spherical surface
S_Fy=sx.*sy; % y component of F on spherical surface
S_Fz=sy.*sz; % z component of F on spherical surface
%
% Surface area elements: Calculates sphere area and assumes each face on the unit sphere has equal surface area
S=4*pi*r^2; % surface area of sphere
delta_S=S/(n_faces^2); % assumes each surface area element has the same area (n_face^2 faces)
%
% Unwrap elements prior to loop
S_Fx=S_Fx(:);
S_Fy=S_Fy(:);
S_Fz=S_Fz(:);
n_x=n_x(:);
n_y=n_y(:);
n_z=n_z(:);
%
% Use for loop to compute the sum of values of (fdotnhat)*(delta_S)
rhssurface=0; % initialises running sum for the surface integral
for i=1:length(S_Fx)
    F_curr=[S_Fx(i) S_Fy(i) S_Fz(i)]; %current field vector position
    n_curr=[n_x(i) n_y(i) n_z(i)]; %current surface normal vector
    F_dot_nhat=dot(F_curr, n_curr); % calculates F dot nhat
    F_dot_nhat_deltaS=F_dot_nhat*delta_S; % calculates (F dot nhat)*(delta S)
    rhssurface=rhssurface+F_dot_nhat_deltaS; %adds to running rhs sum
end
%
%% Display results for numerical proof
volumedisplay='The volumetric sum of div(F) is %.4f.\n'; % display sentence to be used in fprinft for LHS
surfacedisplay='The surface sum of F dot ncap is %.4f.\n'; % display sentence to be used in fprintf for RHS
% Uses fprintf to display results
fprintf(volumedisplay,lhsvolume)
fprintf(surfacedisplay,rhssurface)
%
