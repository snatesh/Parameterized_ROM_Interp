function M = genMassMatrix2D( level )
% In
%   level       ...     2^level number of grid points in each dimension
% Out
%   A           ...     Mass matrix 

N = 2^level + 1;
h = 1/2^level;

phi = @(x)max(1 - abs(x), 0);
phili = @(x, i)phi((x - i*h)/h);
philili = @(x, y, i, j) phili(x, i).*phili(y, j);
phiProd = @(x, y, i, j, l, k) philili(x,y,i,j).*philili(x,y,l,k);
%% compute grid
[X, Y] = meshgrid(linspace(0, 1, N), linspace(0, 1, N));
gridPts = [X(:), Y(:)];

%% setup mass matrix
M = sparse(N^2, N^2);
for j=0:N-1
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
            M(gIndx, ggIndx) = integral2(@(x, y)phiProd(x, y, i, j, pts(1), pts(2)), max(0, (i - 1)*h), min(1, (i + 1)*h), max(0, (j - 1)*h), min(1, (j + 1)*h), 'Method', 'iterated');
        end
        M(gIndx, gIndx) = integral2(@(x, y)phiProd(x, y, i, j, i, j), max(0, (i - 1)*h), min(1, (i + 1)*h), max(0, (j - 1)*h), min(1, (j + 1)*h), 'Method', 'iterated');
    end
end
end


function gIndx = getIndx(pts, N)

gIndx = pts(1)*N + pts(2) + 1;

end

