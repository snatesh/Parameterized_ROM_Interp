function [ A, f, gridPts ] = genOperators2DLogNormal( a, level )
% Dirichlet top
% Homogeneous Neumann side
% Neumann 1 at bottom
% In
%   a           ...     coefficient function
%   level       ...     2^level number of grid points in each dimension
% Out
%   A           ...     discretized Laplace operator
%   f           ...     right-hand side
%   gridPts     ...     grid points

N = 2^level + 1;
h = 1/2^level;

dphili = @(x, i)1/h*dphi((x - i*h)/h);
phi = @(x)max(1 - abs(x), 0);
phili = @(x, i)phi((x - i*h)/h);
dphiProd = @(x, y, i, j, l, k)(dphili(x, i).*dphili(x, l).*phili(y, j).*phili(y, k) + phili(x, i).*phili(x, l).*dphili(y, j).*dphili(y, k)); %dphilili(x, y, i, j).*dphilili(x, y, l, k);
philili = @(x, y, i, j)phili(x, i).*phili(y, j);

%% compute grid
[X, Y] = meshgrid(linspace(0, 1, N), linspace(0, 1, N));
gridPts = [X(:), Y(:)];

%% setup discretized laplace operator
A = sparse(N^2, N^2);
f = sparse(N^2, 1);
for j=0:N-2
    for i=0:N-1
        gIndx = getIndx([i, j], N);
        ptsList = [
            i - 1, j;
            i + 1, j;
            i, j - 1;
            i, j + 1;
            i - 1, j - 1;
            i + 1, j - 1;
            i - 1, j + 1;
            i + 1, j + 1
            ];
        for ptsIter=1:size(ptsList, 1)
            pts = ptsList(ptsIter, :);
            if(any(pts < 0) || any(pts > N-1))
                continue
            end
            ggIndx = getIndx(pts, N);
            A(gIndx, ggIndx) = integral2(@(x, y)a(x, y).*dphiProd(x, y, i, j, pts(1), pts(2)), max(0, (i - 1)*h), min(1, (i + 1)*h), max(0, (j - 1)*h), min(1, (j + 1)*h), 'Method', 'iterated');
        end
        A(gIndx, gIndx) = integral2(@(x, y)a(x, y).*dphiProd(x, y, i, j, i, j), max(0, (i - 1)*h), min(1, (i + 1)*h), max(0, (j - 1)*h), min(1, (j + 1)*h), 'Method', 'iterated');
        % impose Neumann 1 on base
        if(j == 0)
            % added coefficient field to load functional
            f(gIndx) = integral(@(x)a(x,0).*phili(x, i), max(0, (i - 1)*h), min(1, (i + 1)*h));
        end
    end
end
% Top: Dirichlet 0
j = N-1;
for i=0:N-1
    gIndx = getIndx([i, j], N);
    A(gIndx, gIndx) = 1;
    f(gIndx) = 0;
end
end

function gIndx = getIndx(pts, N)

gIndx = pts(1)*N + pts(2) + 1;

end

function f = dphi(x)
[n, m] = size(x);
x = x(:);
f = zeros(length(x), 1);
f(x < 0) = 1;
f(x >= 0) = -1;
f(x < -1 | x > 1) = 0;
f = reshape(f, n, m);
end
