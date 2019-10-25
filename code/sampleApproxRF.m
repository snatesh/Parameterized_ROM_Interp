% plots for KL-expansion of Gaussian random field
clear all; close all; clc;
addpath(genpath('./ellip'));
set(groot, 'defaultLineLineWidth', 1.5);
set(groot,'defaultLineMarkerSize',10);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',20);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',20);

%% KL of Gaussain RF
rng('default');
rho = 1;
sigma = 2;
% covariance kernel
sepexpk = @(x,y) prod(exp(-abs(x-y)),2);
sqexpk = @(x,y) sigma^2*exp(-norm(x-y)^2/(2*rho^2));
n = 30;
[S,X,Y] = kernMat2D(sqexpk, n, 0,1);

% eigendecomposition
[U, V] = eig(S);
[V, I] = sort(diag(V), 'descend');
U = U(:, I);

% plot eigenfunctions
figure(1);
subplot(2,2,1);
surfl(X,Y,reshape(sqrt(V(1))*U(:,1),n,n)); grid on;
xlabel('x'); ylabel('y'); zlabel('$\psi_1$');
subplot(2,2,2);
surfl(X,Y,reshape(sqrt(V(1))*U(:,10),n,n)); grid on;
xlabel('x'); ylabel('y'); zlabel('$\psi_{10}$');
subplot(2,2,3);
surfl(X,Y,reshape(sqrt(V(1))*U(:,30),n,n)); grid on;
xlabel('x'); ylabel('y'); zlabel('$\psi_{30}$');
subplot(2,2,4);
surfl(X,Y,reshape(sqrt(V(1))*U(:,40),n,n)); grid on;
xlabel('x'); ylabel('y'); zlabel('$\psi_{40}$');

%%
% finde candidate eigenvectors and eigen functions 
% for rank truncation of square exp kernel to rank=2
Is = find(sqrt(V)./sqrt(V(1))>1e-3);
ni = length(Is);
if ni >= 2
    Is=Is(1);
    rIs =  randn(length(Is),1)
    subplot(1,2,1);
    G = U(:,Is)*diag(sqrt(V(Is)))*rIs;    
    surf(X,Y,reshape(G,n,n));
    xlabel('$x$');ylabel('$y$'); zlabel('$G(x;\lambda)$');
    title('KL Expansion of Gaussian Random Field Truncated After 1 Term');
    % compute log normal
    K = exp(G);
    subplot(1,2,2);
    surf(X,Y,reshape(K,n,n));
    xlabel('$x$'); ylabel('$y$'); zlabel('$k(x;\lambda)$');
    title('Corresponding Log-Normal Field');
end

function [K,X,Y] = kernMat2D(covk,n,a,b)
    x = linspace(a,b,n);
    [X,Y] = meshgrid(x);
    K = zeros(n^2,n^2);
    for i = 1:n^2
        x1 = [X(i),Y(i)];
        for j = 1:n^2
            x2 = [X(j),Y(j)];
            K(i,j) = covk(x1,x2);
        end
    end
end
