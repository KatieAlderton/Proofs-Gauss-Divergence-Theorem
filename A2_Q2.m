%{
Name: Katie Alderton
Date: 19/04/2024
Description: Analytical Proof

Input(s):
- Symbolic variables (x, y, z, r, theta, phi)
- Vector field F, F(x,y,z)=[xz^2 xy yz]


Output(s): 
- Volume integral to compute flux (using spherical coordinates and
divergence)  
- Surface integral to compute flux
- A displayed text to check calculated volume and surface integrals are
equal, hence proving the divergence theorem via analytical methods.
%}
%
clear 
close all
clc
%
% Define symbolic variables using syms 
syms x y z r theta phi  
%
% Declare each variable as real
assume(x,'real');
assume(y,'real');
assume(z,'real');
assume(theta,'real');
assume(phi,'real');
%
% Define the vector field F(x,y,z)
Fx=x*(z)^2;
Fy=x*y;
Fz=y*z;
F=[Fx Fy Fz]; % F defines the vector field F(x,y,z)
%
%% LHS volume integral
% Define Cartesian coordinates in terms of spherical coordinates
x_sph=r*cos(theta)*sin(phi);
y_sph=r*sin(theta)*sin(phi);
z_sph=r*cos(phi);
%
% Convert these coordinates to vector form so they can be used in the
% jacobian function
cart_to_sph=[x_sph y_sph z_sph];
sph=[r phi theta];
%
%  Compute the jacobian matrix and determinant
J_matrix=jacobian(cart_to_sph,sph);
detJ=det(J_matrix);
%
% Determine the divergence of F in cartesian coordinates
div_F_cart=divergence(F,[x y z]); 
%
% Convert div_F to spherical coordinates via appropriate subsitutions
div_F_sph=subs(div_F_cart,[x,y,z],[cart_to_sph]); 
%
% Determine the volume integral (integrand=1*detJ) using of div f using spherical
% coordinates
% NOTE: the spherical coordinates are: 0<=r<=1, 0<=phi<=pi, 0<=theta<=2pi
integral1=int(div_F_sph*detJ,r,0,1);
integral2=int(integral1,phi,0,pi);
integral3=int(integral2,theta,0,2*pi);
% Result from LHS of Gauss's divergence theorem:
volumeintegral=integral3;
%
%% RHS surface integral 
% Define the vector normal to a spherical surface
n_vector=[x y z];
%
% Normalise n_vector into a unit vector (n_norm)
n_hat=n_vector/norm(n_vector);
%
% Determine the dot of vector field F(x,y,z) and the normal in unit vector
% form
F_dot_nhat_cart=dot(F,n_hat);
%
% Use substitutiomns to convert F_dot_norm_cart to spherical coordinates 
F_dot_nhat_sph=subs(F_dot_nhat_cart,[x,y,z],[cart_to_sph]); 
%
% Set up integrand:
%
% Define the integrand for the surface integral by multiplying F.nhat with
% Jacobian determinant 
integrand=F_dot_nhat_sph*detJ; 
% sub r=1 into integrand expression
integrand_r1=subs(integrand,r,1);
% Simplify this integral using the simplify function prior to proceeding
% with integration
integrand_r1_simple=simplify(integrand_r1);
%
% Determine the surface integral of F_dot_nhat_sph using spherical
% coordinates, note as we substitute r=1 above we only need to worry about
% the double integral with respect to phi and theta respectively. 
%
% for simplicity I have split up the double integral as follows
surf_integral1=int(integrand_r1_simple, phi, 0, pi);
surf_integral2=int(surf_integral1,theta, 0, 2*pi);
% Which gives the final computed surface integral:
surfaceintegral=surf_integral2;
%
%% The Proof
% Set up sentences for fprintf 
% if equal:
sentence_equal='The 2 integrals are equal, with a value of %.4f: Divergence Theorem shown.\n';
sentence_not='The volume integral is %.4f while the surface integral is %.4f. There is a mistake somewhere.\n';
%
% use an if else statement to validate that the LHS and RHS of the
% divergence theorem return equal results
%
if volumeintegral==surfaceintegral
    fprintf(sentence_equal,surfaceintegral)
else
    fprintf(sentence_not,volumeintegral,surfaceintegral)
end


