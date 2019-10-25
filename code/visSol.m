% sachin natesh
% time snaps of solution
clear all; close all;clc;
addpath(genpath('./'));
set(groot, 'defaultLineLineWidth', 1.5);
set(groot,'defaultLineMarkerSize',10);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',20);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',20);

level = 5;
n = 2^level+1;
ops = load('./ellip/operators200P1Blocks1_level5.mat');
M = load('./parabolic/MassMatrix_level5.mat');
M = M.M;
[X,Y] = meshgrid(linspace(0,1,n));
A = ops.ACell{15}; f = ops.f{15};
u0 = randn(size(A,1),1);
t0 = 0;
tf = 0.01;
nsteps = 5000;
dt = (tf-t0)/nsteps;
%t0 = 0; tf = .01; dt = 0.00001;
[U,t] = backwardEuler(u0,A,M,f,t0,tf,dt,true);
%%
% up = max(max(U));
% lo = min(min(U));
% for j = 1:length(t)
%     u = reshape(U(:,j),n,n);
%     surf(X,Y,u);
%     zlbl = strcat('$u(x,t=',num2str(t(j)),')$');
%     xlabel('$x$'); ylabel('$y$'); zlabel(zlbl);
%     axis([0,1,0,1,lo,up])
%     drawnow; pause(0.001);
% end


% for j = 1:4
%     u = reshape(U(:,(j-1)*25+1),n,n);
%     subplot(2,2,j);
%     surfl(X,Y,u);
%     zlbl = strcat('$u(x,t=',num2str(t((j-1)*25+1)),')$');
%     xlabel('$x$'); ylabel('$y$'); zlabel(zlbl,'Rotation',0);
% end
