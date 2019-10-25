% cross-validation of ROM for non-affine/log-normal coeff field

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

rng('default');
% load all precomputed operators (200 of them)
ops = load('./ellip/operators200P1Blocks1_level5.mat');
Nops = size(ops.ACell,2);
% load mass matrix
M = load('./parabolic/MassMatrix_level5.mat');
M = M.M;
% random initial data for parabolic pde
u0 = randn(size(M,1),1);
% load operating points
xis = load(['./ellip/xis200.mat']);
xis = xis.xis;
% size of full order model
Nf = size(ops.ACell{1}, 1);
% number of operating points xi, i.e. # of RBs to genereate
Nr = 50;
% of these 50, find 5 that are near each other for interpolation.
Ninterp = 5;
% number of tests
Ntests = 1;%10;
% number of validations
Nvalidations = 1;%20;
% number of basis vectors after truncation
maxb = 10; 
Nphi_list = maxb;%1:maxb;
speedups = zeros(maxb,1);
errQevals = zeros(maxb,1);
t0 = 0; tf = 1; nsteps = 100; dt = tf/nsteps;
% cross validate over 20 permuations of training and test sets
for nv = 1:Nvalidations
    % random indices from which we determine test and train data
    q = randperm(Nops);
    speedup = zeros(maxb,1);
    errQeval = zeros(maxb,1);    
    for Nphi = Nphi_list
        % initialize container for RBs at each operating point
        Phis = zeros(Nf, Nphi, Nr);
        %t0 = 0; tf = 1; dt = 0.05;
        % use first Nr random indices as training set 
        % but it becomes Ninterp random indices after we find nearest RBs
        j = 1;
        for i = q(1:Nr)
            % solve parabolic pde for ith operating point
            A = ops.ACell{i}; f = ops.f{i};
            i
            [U,~] = backwardEuler(u0,A,M,f,t0,tf,dt,false);     
            % generate RB by truncating svd of [u_1 ... u_nt]
            [Phi, ~, ~] = svd(U, 0);
            Phis(:,:,j) = Phi(:, 1:Nphi);  
            j = j + 1;
        end
        [PhisNeighbors,I] = findNearestRBs(Phis,Ninterp);
        % clear useless variables
        Phis(:) = []; U(:) = []; Phi(:) = []; A(:) = []; f(:) = [];
        % add operating points to training set
        Lambdas = xis(q(I))';
        % compute projections of training RBs onto tangent space of first Phi 
        % so we can interpolate to get test RBs
        [Gammas,Phi_1] = computeLogMaps(PhisNeighbors,Nf,Nphi,Ninterp);
        % define new operating points (test points)
        newOpPnts = q(Nr+10:Nr+10+Ntests);
        nNewOpPnts = 0;
        k = 1;
        for newOpPnt = newOpPnts
            LambdaNr = xis(newOpPnt);
            % if new op point out of range, continue
            if (LambdaNr < min(Lambdas) || LambdaNr > max(Lambdas))
                continue
            end
            % otherwise, increment number of new op points and interp
            nNewOpPnts = nNewOpPnts + 1;
            % interpolate Gammas to find that for test point
            GammaNr = interpolateGammas(Gammas, Lambdas, xis, Nf, Nphi, Ninterp, newOpPnt);
            % now use exponential map to project back to Grassmanian
            Phi_nr = computeExpMap(GammaNr,Phi_1);
            % project FOM operators and initial cond onto RB
            A = ops.ACell{newOpPnt}; f = ops.f{newOpPnt};
            % including projection ops as part of time for ROM eval
            tic
            Mrb = Phi_nr' * M * Phi_nr;
            Arb = Phi_nr' * A * Phi_nr;
            frb = Phi_nr' * f;
            u0rom = Phi_nr'*u0;
            [Urom,~] = backwardEuler(u0rom,Arb,Mrb,frb,t0,tf,dt,false);
            % evaluate QoI for each time
            Qrom = Urom' * frb;
            timeRom = toc;
            tic
            [Ufom,~] = backwardEuler(u0,A,M,f,t0,tf,dt,false);    
            Qfom = Ufom' * f;
            timeFom = toc;
            % total time and error over test points
            speedup(Nphi) = speedup(Nphi) + timeFom/timeRom;
            errQeval(Nphi) = errQeval(Nphi) + norm(Qrom-Qfom);
            if k == 1
                for j = 1:(nsteps+1)
                    j
                    subplot(2,1,1);
                    pcolor(reshape(Ufom(:,j),n,n));  shading interp;
                    ylabel('Y Index');
                    ttl = strcat("Time = ",num2str(t(j-1)));
                    title(ttl);
                    subplot(2,1,2);
                    pcolor(reshape(Urom(:,j),n,n)); shading interp;
                    ylabel('Y Index'); xlabel('X Index');
                    set(gcf, 'Position',  [0,0, 1920, 1080])
                    drawnow; pause(0.001);
                    %ttl = strcat("Time=",num2str(t(j-1)),'.png');
                    %saveas(gcf,ttl);  
                end
            end
        end
        % average time and error over test points
        if (nNewOpPnts > 0)
            speedup(Nphi) = speedup(Nphi)/nNewOpPnts;
            errQeval(Nphi) = errQeval(Nphi)/nNewOpPnts;
        end
    end
    % total time and error over validations
    speedups = speedups + speedup;
    errQevals = errQevals + errQeval;
    disp(num2str(nv));
end

% average time and error over validations
speedups = speedups/Nvalidations;
errQevals = errQevals/Nvalidations;

figure(1);
semilogy(Nphi_list,errQevals,'-ro','DisplayName',"$||Q_{fom}(0:T) - Q_{rom}(0:T))||$");
legend show;
xlabel('Reduced Basis Dimension');
title('Cross-validated Average 2-norm Error');
figure(2);
plot(Nphi_list,speedups,'-ro','DisplayName','speedup');
legend show;
xlabel('Reduced Basis Dimension');
title('Cross-validated Average Speedup');



% compute distance between points on Grassmanian
function d = grassmanDist(X,Y)
    [~,S,~] = svd(X'*Y,0);
    Theta = diag(acos(diag(S)));
    d = norm(Theta);
end

% find nearest neighbors of Phi_1
function [PhisNearest,I] = findNearestRBs(Phis,nNearest)
    sz = size(Phis);
    PhisNearest = zeros(sz(1),sz(2),nNearest);
    for i = 1:sz(3)
        dists(i) = grassmanDist(Phis(:,:,1),Phis(:,:,i));
    end   
    [~,I] = sort(dists);
    PhisNearest = Phis(:,:,I(1:nNearest));
    I = I(1:nNearest);
end

% project RBs onto tangent space of Phi_1
function [Gammas,Phi_1] = computeLogMaps(Phis,Nf,Nphi,Nr)
    % reference point on Grassmanian
    Phi_1 = Phis(:,:,1);
    % container of projections of Phi_i onto tangent space of Phi_i0
    Gammas = zeros(Nf, Nphi, Nr);
    I = speye(size(Phi_1, 1));
    Phi1_Phi1T = Phi_1 * Phi_1';
    % compute logarithm maps, i.e. projections onto tangent space
    for i = 1:Nr
        Phi_i = Phis(:,:,i);
        % slight variation of eq. 20 in Amsallem and Farhat to 
        % solve lin sys instead of compute inv
        preLogMap = (Phi_i' * Phi_1) \ (Phi_i' - Phi_i' * Phi1_Phi1T);
        [U,S,V] = svd(preLogMap',0);
        Gammas(:,:,i) = U * diag(atan(diag(S))) * V';
    end 
end

% interpolate in tangent space to find projection of RB at test point
% on tangent space to Phi_1
function GammaNr = interpolateGammas(Gammas, Lambdas, xis, Nf, Nphi, Nr, newOpPnt)
    % define test point
    LambdaNr = xis(newOpPnt);
    % initialize logmap for this test point
    GammaNr = zeros(Nf,Nphi);
    % interpolate entries of Gammas to populate GammaNr
    for i = 1:Nf
        for j = 1:Nphi 
            data = zeros(Nr,1);
            for k = 1:Nr
                data(k) = Gammas(i,j,k);
            end
            % extrap included just incase
            GammaNr(i,j) = interp1(Lambdas',data,LambdaNr,'linear','extrap');
        end
    end
end

% exponential map to project interpolated test Gamma (on tangent space to Phi_1)
% back to Grassmanian
function Phi_nr = computeExpMap(GammaNr,Phi_1)
    [Unr,Snr,Vnr] = svd(GammaNr,0);
    cosSnr = diag(cos(diag(Snr)));
    sinSnr = diag(sin(diag(Snr)));
    Phi_nr = Phi_1 * Vnr * cosSnr  + Unr * sinSnr ;
end


% forward euler, unstable except for prohibitive time steps (not used)
function fwdEuler(u0,A,M,f,t0,tf,dt,X,Y)
    t = t0:dt:tf;
    nt = length(t);
    % precompute step operator
    stepOp = M\(M-dt*A);
    % step to final time
    uprev = u0;
    for j = 1:(nt-1)
        surf(X,Y,reshape(uprev,33,33));drawnow;
        u = stepOp*uprev + dt*(M\f);
        uprev = u;
    end
end
