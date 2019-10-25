% precompute operators for block with 1
% and log-normal coefficient field
clear all; close all; clc;

% first define square exp covariance kernel
% then take KL of gaussian RF for several pairs (xi1,xi2) ~ N(0,1)
% last, sample log-normal field and compute operator A and rhs f
% currently working at level 5
level = 5;
% 65^2 points in [0,1]^2 for square exp covariance kernel
n = 65;
rho = 1;
sigma = 2;
% covariance kernel
sqexpk = @(x,y) sigma^2*exp(-norm(x-y)^2/(2*rho^2));
[S,X,Y] = kernMat2D(sqexpk, n, 0,1);

% eigendecomposition
[U, V] = eig(S);
[V, I] = sort(diag(V), 'descend');
U = U(:, I);
% find candidate eigenvectors and eigen functions 
% for rank truncation of square exp kernel to rank=2
Is = find(sqrt(V)./sqrt(V(1))>1e-3);
ni = length(Is);
rng('default');
if ni >= 1
    xis = zeros(200,1);
    tEvalAf = zeros(200,1);
    Iss = 0;
    parfor j = 1:200 
        tic 
        Iss=Is(1);
        % sample standard nomral to get p = (xi1)
        xi =  randn(length(Iss),1); 
        xis(j) = xi;
        G = U(:,Iss)*diag(sqrt(V(Iss)))*xi;    
        % compute log normal
        K = exp(G);
        K = reshape(K,n,n);
        % compute interpolant
        interpK = @(x,y) interp2(X,Y,K,x,y,'cubic');
        [ACell{j}, f{j}] = genOperators2DLogNormal(@(x,y) interpK(x,y),level);
        tEvalAf(j) = toc;
        disp(num2str(j));        
    end
    eval(['save operators200P1Blocks1_level', num2str(level),'.mat ACell f']);
    save('tEvalAf200.mat','tEvalAf');
    save('xis200.mat','xis');
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
