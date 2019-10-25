% backward euler to step parabolic system forward in time
% Input: 
%   - u0 initial condition
%   - A  stiffness matrix
%   - M  mass matrix
%   - f  load vector
%   - t0 start time
%   - tf end time
%   - dt time step
% Output:
%   - U  solution for each time
%   - t  times

function [U,t] = backwardEuler(u0,A,M,f,t0,tf,dt,plt)
    t = t0:dt:tf;
    nt = length(t);
    % precompute step operator
    MpA = M+dt*A;
    stepOp = MpA\M;
    % step to final time
    n = sqrt(size(u0,1));
    U = zeros(n^2,nt);    
    %uprev = u0;
    U(:,1) = u0;  lo = min(u0); up = max(u0);  
    for j = 2:nt
        if plt
            surfl(reshape(U(:,j-1),n,n));  shading interp;
            ttl = strcat("Time = ",num2str(t(j-1)));
            title(ttl);
            set(gcf, 'Position',  [0,0, 1920, 1080])
            ttl = strcat("Time=",num2str(t(j-1)),'.png');
            %saveas(gcf,ttl);
            drawnow;  pause(0.000001); 
            
        end    
        u = stepOp*U(:,j-1) + dt* (MpA \ f);
        U(:,j) = u;
    end
    if plt
        surf(reshape(u,n,n));  shading interp;
        ttl = strcat("Time = ",num2str(t(j)));
        title(ttl);
        set(gcf, 'Position',  [0,0, 1920, 1080])
        ttl = strcat("Time=",num2str(t(j)),'.png');
        %saveas(gcf,ttl);
        drawnow; 
    end     
end
