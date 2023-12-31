function f = chebfun3double(f, op, dom, pref, isEqui)
%CHEBFUN3DOUBLE   CHEBFUN3 constructor for discrete tensor of values
%   This algorithm is equivalent to chebfun3classic for double inputs.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The input is a discrete tensor of values.
pseudoLevel = pref.cheb3Prefs.chebfun3eps;
tech = pref.tech();

f = chebfun3();
if ( ~isEqui && numel(op) == 1 )
    f = constructor(f, @(x,y,z) op + 0*x, dom);
    return;
end

% N.B. We cannot detect if MESHGRID was used to generate values unless we
% know also the (x,y,z) points used to generate those values. If we knew
% beforehand that ALL users WILL generate their tensor of values ONLY from
% meshgrid pts, all we need is to say vals = permute(vals,[2 1 3]); to
% generate a tensor corresponding to ''meshgrid'', in which case a copy of
% ``tensorGrid`` should also be used accordingly. op = permute(op,[2 1 3]);
if ( ~isEqui )
    m = size(op, 1);
    n = size(op, 2);
    p = size(op, 3);
    out = tech.tensorGrid([m, n, p], dom);
    xx = out{1};
    yy = out{2};
    zz = out{3};
else
    % Equispaced points from ndgrid, not meshgrid!
    x = linspace(dom(1), dom(2), size(op, 1));
    y = linspace(dom(3), dom(4), size(op, 2));
    z = linspace(dom(5), dom(6), size(op, 3));
    [xx, yy, zz] = ndgrid(x, y, z);
end

% Calculate a tolerance and find numerical rank to this tolerance: The
% tolerance assumes the samples are generated by NDGRID from a function.
% It depends on the size of the sample tensor, hscale of domain, vscale of
% the samples, condition number of the function, and the accuracy target in
% chebfun3 preferences.
[relTol, absTol] = getTol3D(xx, yy, zz, op, max(size(op)), dom, pseudoLevel);
pref.chebfuneps = relTol;

% Perform 3D ACA with complete pivoting:
factor = 0;
[colsValues, rowsValues, pivotVals2D, ~, tubesValues, pivotVals3D, ~, ~, ...
    ~] = completeACA3D(op, absTol, factor, dom, pref);

sepRank = numel(pivotVals3D); % first separation rank
diagValues2D = cell(sepRank, 1);
for k=1:sepRank
    diagValues2D{k} = diag(1./pivotVals2D{k});
end

% BTD ---> Tucker compression:
[core, colsValues, rowsValues] = btd2tucker(colsValues, rowsValues, ...
    diagValues2D, pivotVals3D, absTol);

% Construct a CHEBFUN3 object and call a sampleTest.
if ( ~isEqui )
    f.cols = chebfun(colsValues, dom(1:2), pref);
    f.rows = chebfun(rowsValues, dom(3:4), pref);
    f.tubes = chebfun(tubesValues, dom(5:6), pref);
else
    f.cols = chebfun(colsValues, dom(1:2), 'equi', pref);
    f.rows = chebfun(rowsValues, dom(3:4), 'equi', pref);
    f.tubes = chebfun(tubesValues, dom(5:6), 'equi', pref);
end
f.core = core;
f.domain = dom;
return

end

%%
function [relTol, absTol] = getTol3D(xx, yy, zz, vals, grid, dom, ...
    pseudoLevel)

relTol = 2*grid^(4/5) * pseudoLevel; % this should be vscale and hscale invariant
vscale = max(abs(vals(:)));
[m,n,p] = size(vals);
% Remove some edge values so that df_dx, df_dy and df_dz have the same size. 
% xx changes in the first mode:
df_dx = diff(vals(:, 1:n-1, 1:p-1), 1, 1) ./ diff(xx(:, 1:n-1, 1:p-1), 1, 1);
% yy changes row-wise (2nd mode):
df_dy = diff(vals(1:m-1, :, 1:p-1), 1, 2) ./ diff(yy(1:m-1, :, 1:p-1), 1, 2);
% zz changes tube-wise (3rd mode):
df_dz = diff(vals(1:m-1, 1:n-1, :), 1, 3) ./ diff(zz(1:m-1, 1:n-1, :), 1, 3);
gradNorms = [max(abs(df_dx(:))), max(abs(df_dy(:))), max(abs(df_dz(:)))];
% A vector of gradient information over the domain.
if ( isempty(gradNorms) )
    % This happens if the input in not a trivariate function in which case
    % we basically disable using gradient information:
    gradNorms = 1;
end
domDiff = [diff(dom(1:2)) diff(dom(3:4)) diff(dom(5:6))];
absTol = max(max(gradNorms.*domDiff), vscale) * relTol;
% absTol should depend on the vscale of the function while it also uses
% derivative information to prevent issues like the one mentioned in
% https://github.com/chebfun/chebfun/issues/1491.

end

%%
function [colsBtd, rowsBtd, pivotValues2D, pivotIndices2D, fibers, ...
    pivotValues3D, pivotIndices3D, ifail3D, ifail2D] = completeACA3D(A, ...
    tol, factor, dom, pref)
%   Non-adaptive (fixed-size) MACA, i.e., a 3D analogue of Gaussian 
%   elimination with complete pivoting.
%
%   INPUTS:     A:        A given tensor of function values at 3D chebpts.
%
%               tol:      A given tolerance on the magnitude of the pivots.
%
%               factor:   The ratio between the width of A and the number 
%                         of iterations (rank) allowed in Gaussian elimination.
%
%  OUTPUTS:     colsBtd: A cell-array containing skeleton columns of slices 
%                        in block term decomposition.
%
%               rowsBtd: A cell-array containing skeleton rows of slices in
%                        block term decomposition.
%
%               pivotValues2D: A cell array containing the values of pivots
%                              in 2D ACAs.
%
%               pivotIndices2D: A cell array containing indices i.e.,
%                               locations of of 2D pivot points.
%
%                       fibers: A matrix of size n1 x ITER. Each of its columns
%                             contains the values of the updated (residue) 
%                             tensor at the pivot fiber.
%
%               pivotValues3D: A row vector containing the values of the 
%                             pivot entries during the iterations. 
%
%               pivotIndices3D: A matrix of size rank x 3, where rank = iter. 
%                         Each of its rows contain the index of one 3D pivotValues.
%
%                      ifail: We fail if iter >= (width/factor).

% Developer Note: The output of this code should satisfy the following 
% slice decomposition:
%       AA \approx temp, 
% where AA is a copy of A from input, and temp is computed as follows:
%   temp = zeros(size(A2)); 
%   for i = 1:3, 
%    temp = temp + chebfun3.outerProd(slices(:,:,i),fibers(:,i)./pivotValues3D(i));
%   end
%   norm(AA(:) - temp(:))
%
% An analogous BTD decomposotion should also hold.

pseudoLevel = pref.cheb3Prefs.chebfun3eps;
tech = pref.tech();

% Set up output variables.
[n1, n2, n3] = size(A);
width = min(n3, n1*n2);        % Use to tell us how many pivots we can take
                               % See Developer note in the following.
pivotValues3D = zeros(1);      % Store an unknown number of Pivot values
pivotIndices3D = zeros(1, 3);  % Store (col, row, tube) = entries of pivot location
ifail3D = 1;                   % Assume we fail in 3D ACA
ifail2D = 1;                   % Assume we also fail in the 2D ACAs
globalTol = [];
sliceDim = [1 2];              % See Developer note in the following.

% Main algorithm
iter = 0;                      % Count number of interpolated rows/slices.
[infNorm, ind] = max(abs(reshape(A, numel(A), 1))); % Complete pivoting
[col, row, tube] = ind2sub(size(A), ind);

scl = infNorm;
% If the function is the zero function.
if ( scl == 0 )
    pivotValues3D = 0;
    fibers = 0;
    colsBtd{1} = 0;
    rowsBtd{1} = 0;
    pivotValues2D = 0;
    pivotIndices2D = [0 0];
    ifail3D = 0;
    ifail2D = 0;
else
    fibers(:,1) = zeros(size(A, 3), 1);
    colsBtd{1} = zeros(size(A, sliceDim(1)), 1);
    rowsBtd{1} = zeros(size(A, sliceDim(2)), 1);
    pivotValues2D{1} = 0;
    pivotIndices2D{1} = [0 0];
    slices(:,:,1) = zeros(size(A,sliceDim(1)), size(A, sliceDim(2)), 1);
end
dom2D = dom(1:4);

while ( ( infNorm > tol ) && ( iter < width / factor) ...
        && ( iter < width ) )
    fibers(:, iter+1) = A(col, row, :);  % Extract skeleton tubes. Each
    % column in "fibers" is N3 x 1.
    
    slices(:,:,iter+1) = A(:,:,tube);    % Extract skeleton slices to be 
    % decomposed further. Each slice in "slices" is N1 x N2.
    
    % Developer Note: As the above lines show, we are always separating the
    % last variable from the first two. The point is that the function 
    % handle has already been permuted outside this subroutine and 
    % therefore the tensor A here contains values of the permuted function.
    % In this sense, here we are in essense separating the variable chosen 
    % by the dimension clustering step, and not necessarily the last 
    % variable z.
    
    PivVal3D = A(col, row, tube);        % = f(X(col), Y(row), Z(tube)) in
                                         % the 1st iteration.
    
    % Use the first slice to compute globalTol for 2D ACAs applied to all 
    % slices.
    if ( iter == 0 )
        out = tech.tensorGrid([n1, n2], dom2D);
        xx2D = out{1};
        yy2D = out{2};
        globalTol = GetTol2D(xx2D, yy2D, slices(:,:,1), dom2D, pseudoLevel);
    end
    
    % Apply 2D ACA to each slice to form columns and rows in block term
    % decomposition:
    [colsBtd{iter+1}, pivotValues2D{iter+1}, rowsBtd{iter+1}, ...
        pivotIndices2D{iter+1}, ifail2DIter] = ...
        chebfun2ACA(slices(:, :, iter+1), globalTol, factor);
    
    % Developer Note: Since we use globalTol for slices after 1st 
    % iteration, it might be that these 2D ACA's don't fail, while with a 
    % localTol they would fail. So, it is auaully the 1st slice which shows
    % whether or not we got the 2D rank right. So, just use that one:
    if ( iter == 0 )
        ifail2D = ifail2DIter;
    end
    
    % Update the tensor, i.e., compute the residual tensor:
    A = A - chebfun3.outerProd(colsBtd{iter+1} * ...
        (diag(1./(pivotValues2D{iter+1}))) * (rowsBtd{iter+1}.'), ...
        fibers(:,iter+1) ./ PivVal3D);
    % Equivalently, we have: 
    % A = A - chebfun3.outerProd(slices(:,:,iter+1),fibers(:,iter+1)./PivVal3D);
    
    % Keep track of progress in 3D ACA:
    iter = iter + 1;              % One more fiber and slice are removed from A
    pivotValues3D(iter) = PivVal3D;           % Store pivot value in 3D ACA
    pivotIndices3D(iter, :)=[col row tube];    % Store pivot location in 3D ACA
    
    % Find next 3D pivot value and its location:
    [infNorm, ind] = max(abs(A(:)));
    [col, row, tube] = ind2sub(size(A), ind);
end

if ( ( iter > 0 ) && ( all(pivotValues2D{iter} == 0) ) )
    % If the last 2D pivot was zero, remove it
    colsBtd=colsBtd(1:iter-1);
    rowsBtd=rowsBtd(1:iter-1);
    pivotValues2D = pivotValues2D(1:iter-1);
    pivotIndices2D = pivotIndices2D(1:iter-1);
    fibers = fibers(:, 1:iter-1);
    pivotValues3D = pivotValues3D(1:iter-1);
    pivotIndices3D=pivotIndices3D(1:iter-1,:);
    infNorm = 0; % If the last 2D pivot was zero, infNorm will be NaN and 
    % it makes the next statement result in ifail3D \neq 0; So, put it 0 to
    % get ifail3D = 0;
end

if ( infNorm <= tol )
    ifail3D = 0;                               % We didn't fail in 3D ACA.
    if ( iter == 0 )
        ifail2D = 0;    
    end
end
if ( iter > width/factor )
    ifail3D = 1;                               % We did fail in 3D ACA.
end

end
%%

%%
function absTol = GetTol2D(xx, yy, vals, dom, pseudoLevel)
% GETTOL2D   Calculate a tolerance for the Chebfun2 constructor.
%
%  This is the 2D analogue of the tolerance employed in the chebtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (4/5) exponent. 

[m, n] = size(vals); 
grid = max(m, n);
relTol = grid^(4/5) * pseudoLevel; % this should be vscale and hscale invariant

% Remove some edge values so that df_dx and df_dy have the same size. 
% xx is generated by ndgrid, i.e., xx changes in the first mode:
dfdx = diff(vals(:, 1:n-1), 1, 1) ./ diff(xx(:, 1:n-1), 1, 1);
% yy is generated by ndgrid, i.e., yy changes row-wise (2nd mode):
dfdy = diff(vals(1:m-1, :), 1, 2) ./ diff(yy(1:m-1, :), 1, 2);
gradNorms = [max(abs(dfdx(:))), max(abs(dfdy(:)))];
% A vector of gradient information over the domain.
if ( isempty(gradNorms) )
    % This happens if the input in not a trivariate function.
    gradNorms = 1;
end

vscale = max(abs(vals(:)));
domDiff = [diff(dom(1:2)) diff(dom(3:4))];
absTol = max(max(gradNorms.*domDiff), vscale) * relTol;
end

%%
function [col, pivotVals, row, pivotLoc, ifail2D] = chebfun2ACA(op, ...
    tol, factor)
% Perform GE with complete pivoting:

if ( factor ~= 0 )
    % FACTOR in the 3D steps is either 0 (in case of constructionFromDoubles) 
    % or 2sqrt(2) otherwise. For 2D steps however, we are happy with 
    % FACTOR = 0 or 2 as in Chebfun2. This IF conditional, makes it
    % possible to rewrite FACTOR in 2D steps, but at the same time keeping 
    % it zero fro constructionFromDoubles.
    factor = 2;
end
[pivotVals, pivotLoc, row, col, ifail2D] = completeACA2D(op, tol, factor); 
end


%%
function [pivotValue, pivotElement, rows, cols, ifail2D] = ...
    completeACA2D(A, tol, factor) 
% 2D ACA with complete pivoting which is the continuous analogue of 
% Gaussian elimination with complete pivoting.
% We attempt to adaptively find the numerical rank of function in the 2D
% level. This is _almost_ the same as the one in chebfun2/constructor.


% Set up output variables.
[nx, ny] = size(A);
width = min(nx, ny);        % Use to tell us how many pivots we can take.
pivotValue = zeros(1);      % Store an unknown number of Pivot values.
pivotElement = zeros(1, 2); % Store (j,k) entries of pivot location.
ifail2D = 1;                  % Assume we fail.

% Main algorithm
zRows = 0;                  % count number of zero cols/rows.
[infNorm, ind] = max(abs(reshape(A, numel(A), 1)));
[row, col] = chebfun3.myind2sub(size(A) , ind);

% Bias toward diagonal for square matrices (see reasoning below):
if ( ( nx == ny ) && ( max(abs(diag(A))) - infNorm) > -tol )
    [infNorm, ind] = max(abs(diag(A)));
    row = ind;
    col = ind;
end

scl = infNorm;
% If the function is the zero function
if ( scl == 0 )
    pivotValue = 0;
    cols = 0;
    rows = 0;
    ifail2D = 0;
else
    cols(:,1) = zeros(size(A, 1), 1);
    rows(1,:) = zeros(1, size(A, 2));
end

while ( ( infNorm > tol ) && ( zRows < width / factor) ...
        && ( zRows < min(nx, ny) ) )

    cols(:, zRows+1) = A(:, col);             % Extract skeleton columns
    rows(zRows+1, :) = A(row, :);             % Extract skeleton rows
    PivVal = A(row, col);
    A = A - cols(:, zRows+1)*(rows(zRows+1,:)./PivVal); % One step of GE
    
    % Keep track of progress.
    zRows = zRows + 1;                       % One more row is interpolated
    pivotValue(zRows) = PivVal;              % Store value of 2D pivot
    pivotElement(zRows, :)=[row col];        % Store index of 2D pivot
    
    % Find value and index of next 2D pivot
    [infNorm, ind] = max(abs(A(:))); % Slightly faster
    [row, col] = chebfun3.myind2sub(size(A), ind);
    
    % Have a bias towards the diagonal of A, so that it can be used as a test
    % for nonnegative definite functions. (Complete GE and Cholesky are the
    % same as nonnegative definite functions have an absolute maximum on the
    % diagonal, except there is the possibility of a tie with an off-diagonal
    % absolute maximum. Bias toward diagonal maxima to prevent this.)
    if ( ( nx == ny ) && ( max(abs(diag(A))) - infNorm) > -tol )
        [infNorm, ind] = max(abs(diag(A)));
        row = ind;
        col = ind;
    end
end

if ( infNorm <= tol )
    ifail2D = 0;                               % We didn't fail in 2D ACA
end
if ( zRows >= (width/factor) )
    ifail2D = 1;                               % We did fail in 2D ACA
end

rows = rows.';                               % To unify all the columns, 
                                             % rows and tubes, store 
                                             % skeleton rows also as column 
                                             % vectors.
end

%%
function [core, colsTucker, rowsTucker] = btd2tucker(colsValues, ...
    rowsValues, diagValues2D, pivotVals3D, absTol)
allCols = []; 
allRows = []; 
allDiags = [];
nn = numel(pivotVals3D);
sizeIndex = zeros(nn, 1);

for kkk = 1:nn
    sizeIndex(kkk) = size(rowsValues{kkk}, 2);
    allCols = [allCols, colsValues{kkk}];
    allRows = [allRows, rowsValues{kkk}];
    allDiags = [allDiags; diag(diagValues2D{kkk})];
end
sizeIndex = cumsum([0; sizeIndex]);
        
% Compress allCols and allRows:
if ( size(allCols, 2) > 1 )
    [Su, ~, Vu, colsTucker] = completeACA2D(allCols, absTol, 0);
    % factor = 0, because we want the ACA to be applied even if op is not 
    % low-rank. In contrast to Chebfun2, we now have 
    % allCols = colsTucker * diag(1./Su) * Vu'.
    % Moreover, we use the same absTol to chop columns, the same tolerance
    % used to chop fibers so that ranks are less inconsistent for symmetric
    % functions.
    Vu = Vu*diag(1./Su);
else
    colsTucker = allCols; 
    Vu = 1;
end

if ( size(allRows, 2) > 1 )
    [Sv, ~, Vv, rowsTucker] = completeACA2D(allRows, absTol, 0);
    Vv = Vv*diag(1./Sv);
else
    rowsTucker = allRows; 
    Vv = 1;
end

% Form the core tensor:
core = zeros(size(colsTucker, 2), size(rowsTucker, 2), nn);
for kkk = 1:nn
    core(:,:,kkk) = Vu(sizeIndex(kkk)+1:sizeIndex(kkk+1), :).' * ...
        diag(allDiags(sizeIndex(kkk)+1:sizeIndex(kkk+1))) * ...
        Vv(sizeIndex(kkk)+1:sizeIndex(kkk+1), :)./pivotVals3D(kkk);
end

if ( (max(abs(allCols(:))) == 0) || (max(abs(allRows(:))) == 0) ) 
    % zero input
    core = 0;
end

end
