% File: sb_eig_01.m
% Date: September 27, 2008
% Author: Jason Moore
clear all
close all
clc
global   CF CR G H IBXX IBXZ IBYY IBZZ IEXX IEXZ IEYY IEZZ IFXX IFYY IHXX IHXZ IHYY IHZZ IRXX IRYY LAMBDAF LAMBDAR M_B M_E M_F M_H M_R RF RR W XB XE XH ZB ZE ZH;
global   T_DELTAF T_DELTAR T_PHI T_THETAR;
%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
CF                              =  0.08;                   % M                   Constant
CR                              =  0.08;                   % M                   Constant
G                               =  9.81;                   % KG                  Constant
H                               =  0.8;                    % M                   Constant
IBXX                            =  9.2;                    % KG*M^2              Constant
IBXZ                            =  2.4;                    % KG*M^2              Constant
IBYY                            =  11;                     % KG*M^2              Constant
IBZZ                            =  2.8;                    % KG*M^2              Constant
IEXX                            =  0.0053;                 % KG*M^2              Constant
IEXZ                            =  0.0;                    % KG*M^2              Constant
IEYY                            =  0.1;                    % KG*M^2              Constant
IEZZ                            =  0.1;                    % KG*M^2              Constant
IFXX                            =  0.1405;                 % KG*M^2              Constant
IFYY                            =  0.28;                   % KG*M^2              Constant
IHXX                            =  0.05892;                % KG*M^2              Constant
IHXZ                            = -0.00756;                % KG*M^2              Constant
IHYY                            =  0.06;                   % KG*M^2              Constant
IHZZ                            =  0.00708;                % KG*M^2              Constant
IRXX                            =  0.0603;                 % KG*M^2              Constant
IRYY                            =  0.12;                   % KG*M^2              Constant
LAMBDAF                         =  0.3141592653589793;     % RAD                 Constant
LAMBDAR                         =  0.3141592653589793;     % RAD                 Constant
M_B                             =  85;                     % KG                  Constant
M_E                             =  3;                      % KG                  Constant
M_F                             =  3;                      % KG                  Constant
M_H                             =  4;                      % KG                  Constant
M_R                             =  2;                      % KG                  Constant
RF                              =  0.35;                   % M                   Constant
RR                              =  0.3;                    % M                   Constant
W                               =  1.02;                   % M                   Constant
XB                              =  0.3;                    % M                   Constant
XE                              =  0.59;                   % M                   Constant
XH                              =  0.9;                    % M                   Constant
ZB                              = -0.9;                    % M                   Constant
ZE                              = -0.8;                    % M                   Constant
ZH                              = -0.7;                    % M                   Constant
%-------------------Construct the Stability Matrix and Calculate Eigenvalues for
%-------------------Various Velocities
delta = 1e-11; % perturbance value
vmax = 10; % max foward velocity of the rear wheel to be calculated
n = 1000; % number of iterations
for i=1:n
    v(i) = (i-1)/n*vmax; % ith velocity
    NU5(i) = -v(i)/RR; % ith angular velocity of the rear wheel 
    % "bfr_evalprimes" computes the derivatives of the state variables, the
    % equations of motion were generated in Autolev
    % compute nominal solution for the ith velocity
    nominal(:,i) = sb_ep_01([zeros(10,1);NU5(i);0;0]);
    % build the stability matrix by numerically calculating the partial
    % derivatives of each differential equation with respect to each state
    % variable
    for j=1:13;
        perturb1 = [zeros(10,1);NU5(i);0;0]; %initialize function input
        perturb2 = [zeros(10,1);NU5(i);0;0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        % solve differential equations for perturbance
        prime1 = sb_ep_01(perturb1);
        prime2 = sb_ep_01(perturb2);
        m(:,j) = (prime1-prime2)./2./delta;  % compute partial derivative
    end
    % reduce stability matrix
    N = [4 7 8 10 12 13]';
    stab=zeros(length(N));
    for k1=1:length(N)
        for k2 = 1:length(N)
            stab(k1,k2)=m(N(k1),N(k2));
        end
    end
    A(:,:,i)=stab;
    % calculate the eigenvalues for the reduced stability matrix
    % calculate the eigenvalues for stability matrix
    [V,D]=eig(stab);
    eigval(1:length(diag(D)),i)=diag(D);
    eigvec(:,:,i)=V;
end
%-------------------Plot the Eigenvalues
figure(1)
hold on
for i=1:length(eigval(:,1))
    plot(v,real(eigval( i,1:length(v))),'.k')
%     plot(v,imag(eigval( i,1:length(v))),'.r')
end
plot(v,zeros(length(v),1),'k')
hold off

% Create color gradient vector
color = colormap(cool(length(v)));
% Plot eigenvalues as a function of speed
figure(2)
hold on
plot([min(min(real(eigval)));max(max(real(eigval)))],[0;0],'k',[0;0],[min(min(imag(eigval)));max(max(imag(eigval)))],'k')
for i=1:length(v)
    plot(real(eigval(:,i)),imag(eigval(:,i)),'.','markeredgecolor',color(i,:),'markersize',5)
end
title('Eigenvalue Loci as a Function of Speed')
xlabel('Re(\lambda) [1/s]')
ylabel('Imag(\lambda) [1/s]')
axis image
box on
hold off
colormap cool
caxis([0 vmax])
colorbar('YTickLabel',...
    {'0 m/s','1','2','3',...
     '4','5','6','7','8','9','10 m/s'})