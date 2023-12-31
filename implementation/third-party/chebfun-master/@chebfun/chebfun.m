classdef chebfun
%CHEBFUN   CHEBFUN class for representing functions on [a,b].
%
%   Class for approximating functions defined on finite, semi-infinite, or
%   doubly-infinite intervals [a,b]. Functions may be smooth, piecewise smooth,
%   weakly singular, or blow up on the interval.
%
% CHEBFUN(F) constructs a CHEBFUN object representing the function F on the
% interval [-1,1]. F may be a string, e.g., 'sin(x)', a function handle, e.g.,
% @(x) x.^2 + 2*x + 1, or a vector of numbers. In the first two instances, F
% should be "vectorized" in the sense that it may be evaluated at a column
% vector of points x(:) in [-1,1] and return an output of size NxM where N =
% length(x(:)). If this is not possible then the flag CHEBFUN(F, 'vectorize')
% should be passed. CHEBFUN(F, 'vectorcheck', 'off') disables the automatic
% checking for vector input. Additionally, F may be a CHEBFUN, in which case
% CHEBFUN(F) is equivalent to CHEBFUN(@(X) FEVAL(F, X)). CHEBFUN() returns an
% empty CHEBFUN object.
%
% CHEBFUN(F, [A, B]) specifies an interval [A,B] on which the CHEBFUN is
% defined, where A and/or B may be infinite. CHEBFUN(F, ENDS), where ENDS is a
% 1x(K+1) vector of unique ascending values, specifies a piecewise smooth
% CHEBFUN defined on the interval [ENDS(1), ENDS(K+1)] with additional interior
% breaks at ENDS(2), ..., ENDS(K). Specifying these breaks can be particularly
% useful if F is known to have discontinuities. For example,
%   CHEBFUN(@(x) abs(x), [-1, 0, 1]).
% If a domain is passed to the constructor, it should always be the 2nd input.
%
% CHEBFUN(F, N) constructs a CHEBFUN object obtained by interpolating F on an N
% point Chebyshev grid of the second kind in [-1,1]. Note that this is
% different from CHEBFUN(F, 'trunc', N), which is described below.
% 
% CHEBFUN(A) or CHEBFUN(A, 'chebkind', 2), where A is an Nx1 matrix, constructs
% a CHEBFUN object which interpolates the data in A on an N-point Chebyshev grid
% of the second kind (see >> help chebpts). CHEBFUN(A, 'chebkind', 1) and
% CHEBFUN(A, 'equi') are similar, but here the data is assumed to come from a
% 1st-kind Chebyshev or equispaced grid linspace(-1, 1, N), respectively. (In
% the latter case, a smooth interpolant is constructed using an adaptive
% Floater-Hormann scheme [Numer. Math. 107, 315-331 (2007)].). CHEBFUN(F, N) or
% CHEBFUN(F, N, 'chebkind', 2) is equivalent to CHEBFUN(feval(F, chebpts(N)).
%
% CHEBFUN(C, 'coeffs'), where C is an Nx1 matrix, constructs a CHEBFUN object
% representing the polynomial C(1) T_0(x) + ... + C(N) T_(N-1)(x),
% where T_K(x) denotes the K-th Chebyshev polynomial. This is equivalent to
% CHEBFUN({{[], C}}). C may also be an NxM matrix, as described below.
%
% CHEBFUN(F, ...), where F is an NxM matrix or an array-valued function handle,
% returns an "array-valued" CHEBFUN. For example,
%   CHEBFUN(rand(14, 2))
% or
%   CHEBFUN(@(x) [sin(x), cos(x)])
% Note that each column in an array-valued CHEBFUN object is discretized in the
% same way (i.e., the same breakpoint locations and the same underlying
% representation). For more details see ">> help quasimatrix". Note the
% difference between
%   CHEBFUN(@(x) [sin(x), cos(x)], [-1, 0, 1])
% and
%   CHEBFUN({@(x) sin(x), @(x) cos(x)}, [-1, 0, 1]).
% The former constructs an array-valued CHEBFUN with both columns defined on the
% domain [-1, 0, 1]. The latter defines a single column CHEBFUN which represents
% sin(x) in the interval [-1, 0) and cos(x) on the interval (0, 1]. 
%
% CHEBFUN({F1,...,Fk}, ENDS) constructs a piecewise smooth CHEBFUN which
% represents Fj on the interval [ENDS(j), END(j+1)]. Each entry Fj may be a
% string, function handle, or vector of doubles. For example
%   CHEBFUN({@(x) sin(x), @(x) cos(x)}, [-1, 0, 1])
%
% CHEBFUN(F, PREF) or CHEBFUN(F, [A, B], PREF) constructs a CHEBFUN object from
% F with the options determined by the CHEBFUNPREF object PREF. Construction
% time options may also be passed directly to the constructor in the form
% CHEBFUN(F, [A, B], PROP1, VAL1, PROP2, VAL2, ...). (See CHEBFUNPREF for
% details of the various preference options and their defaults.). In
% particular, CHEBFUN(F, 'splitting', 'on') allows the constructor to
% adaptively determine breakpoints to better represent piecewise smooth
% functions F. For example,
%   CHEBFUN(@(x) sign(x - .3), [-1, 1], 'splitting', 'on')
% CHEBFUN(F, 'extrapolate', 'on') prevents the constructor from evaluating the
% function F at the endpoints of the domain.
%
% If PROP/VAL and PREF inputs are mixed in a single constructor call, the
% preferences determined by the PROP/VAL inputs take priority over those
% determined by PREF.  At most one PREF input may be supplied to the
% constructor at any time.
%
% CHEBFUN(F, 'trunc', N) returns a smooth N-point CHEBFUN constructed by
% computing the first N Chebyshev coefficients from their integral form, rather
% than by interpolation at Chebyshev points.
%
% CHEBFUN(F, 'trig') constructs a CHEBFUN object representing a smooth and
% periodic function F on the interval [-1,1]. The resulting CHEBFUN is
% represented using a Fourier series. All operations done on F should preserve
% smoothness and periodicity, otherwise results are casted into chebfuns 
% represented by Chebyshev rather than Fourier series. Similar options
% as discussed above may be combined with the 'trig' flag, with exception to
% the 'chebkind' and 'splitting' flags.
%
% CHEBFUN(F, 'periodic') is the same as CHEBFUN(F, 'trig').
%
% CHEBFUN(F, ..., 'doubleLength') constructs a chebfun with twice the
% polynomial degree of that chosen by the Chebfun constructor algorithm.
% If LENGTH(F) == N, then LENGTH(CHEBFUN(F,'doubleLength')) == 2*N-1.
%
% CHEBFUN --UPDATE can be used to update to the latest stable release of CHEBFUN
% (obviously an internet connection is required!). CHEBFUN --UPDATE-DEVEL will
% update to the latest development release, but we recommend instead that you
% checkout from the Github repo https://github.com/chebfun/chebfun/. See
% CHEBFUN.UPDATE() for further details.
%
% See also CHEBFUNPREF, CHEBPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBFUN Class Description:
%
% The CHEBFUN class is for representations of piecewise functions on the
% interval [a,b].
%
% The CHEBFUN class is the main user interface. We do not expect users to
% directly invoke any objects below this level.
%
% A CHEBFUN object consists of a collection of FUN objects. There are two main
% tasks for the CHEBFUN constructor: (1) parse the user input, and (2) correctly
% piece together FUN objects to form a global approximation. If the input
% function is globally smooth then the resulting CHEBFUN contains a single FUN
% object. If the input is not smooth, or breakpoints are passed to the
% constructor, CHEBFUN must determine appropriate breakpoints and return a
% piecewise smooth CHEBFUN with multiple FUN objects.
%
% This is a user-level class, and all input arguments should be thoroughly
% sanity checked.
%
% Class diagram: [ADchebfun] <>-- [CHEBFUN] <>-- [<<fun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        % DOMAIN of definition of a CHEBFUN object. If K = length(F.DOMAIN) is
        % greater than 1 then the CHEBFUN is referred to as a "piecewise".
        % CHEBFUN. The first and last values of this vector define the left and
        % right endpoints of the domain, respectively. The other values give the
        % locations of the interior breakpoints that define the domains of the
        % individual FUN objects comprising the CHEBFUN. The entries in this
        % vector should be strictly increasing.
        domain              % (1x(K+1) double)

        % FUNS is a cell array containing the FUN objects that comprise a
        % piecewise CHEBFUN. The kth entry in this cell is the FUN defining
        % the representation used by the CHEBFUN object on the open interval
        % (F.DOMAIN(k), F.DOMAIN(k+1)). If M = size(f.funs, 2) is greater than
        % 1, then the CHEBFUN object is referred to as "array valued".
        funs                % (Kx1 cell array of FUN objects)
        
        % POINTVALUES Values of the function at the break points.
        pointValues = [];      % (1 x (K+1) double)

        % ISTRANSPOSED determines whether a (possibly array-valued) CHEBFUN F
        % should be interpreted as a collection of "column" CHEBFUN objects (if
        % F.ISTRANSPOSED == 0, the default), which are considered (infxM)
        % arrays, or "row" CHEBFUN objects (if F.ISTRANSPOSED == 1), which are
        % (Mxinf) arrays. This difference is only behavioral; the other
        % properties described above are _NOT_ stored differently if this flag
        % is set.)
        isTransposed = 0;   % (logical)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = chebfun(varargin)
            % The main CHEBFUN constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            if ( iscell(varargin{1}) && ...
                    all(cellfun(@(x) isa(x, 'fun'),  varargin{1})) )
                % Construct a CHEBFUN from a cell array of FUN objects.               
                %  Note, this is not affected by the input parser (see error
                %  message below) and must be _fast_ as it is done often.
                if ( nargin > 1 )
                    error('CHEBFUN:CHEBFUN:chebfun:nargin', ...
                     'Only one input is allowed when passing an array of FUNs.')
                end
                % Assign the cell to the .FUNS property:
                f.funs = varargin{1};
                % Collect the domains together:
                dom = cellfun(@(fun) get(fun, 'domain'), f.funs, ...
                    'uniformOutput', false);                
                f.domain = unique([dom{:}]);                
                % Update values at breakpoints (first row of f.pointValues):
                f.pointValues = chebfun.getValuesAtBreakpoints(f.funs, f.domain);
                return
            end
                       
            % Parse inputs:
            [op, dom, data, pref, flags] = parseInputs(varargin{:});
                        
            if ( flags.done )
                % An update was performed. Exit gracefully:
                throwAsCaller(MException('', ''))
            end
            
            if ( isa(op, 'chebfun') && flags.trunc )
                % Deal with the particular case when we're asked to truncate a
                % CHEBFUN:
                f = op;
                
            else
                % Construct from function_handle, numeric, or string input:
                
                % Call the main constructor:
                [f.funs, f.domain] = chebfun.constructor(op, dom, data, pref);
                
                if ( flags.doubleLength )
                    % Using the length of f.funs{1} is okay because the
                    % 'doubleLength' flag is mutually exclusive with 'splitting
                    % on'.
                    pref.techPrefs.fixedLength = 2*length(f.funs{1}) - 1;
                    [f.funs, f.domain] = chebfun.constructor(op, dom, data, pref);
                end

                % Update values at breakpoints (first row of f.pointValues):
                f.pointValues = chebfun.getValuesAtBreakpoints(f.funs, ...
                    f.domain, op);
                
                % Remove unnecessary breaks (but not those that were given):
                [ignored, index] = setdiff(f.domain, dom);
                f = merge(f, index(:).', pref);
               
            end

            if ( flags.trunc )
                % Truncate the CHEBFUN to the required length:
                f = truncate(f, flags.trunc);
            end

        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Hidden = true, Static = false )
        
        % Add breakpoints to the domain of a CHEBFUN.
        f = addBreaks(f, breaks, tol)
                
        % Add breaks at appropriate roots of a CHEBFUN.
        f = addBreaksAtRoots(f, tol)
        
        % Assign columns (or rows) of an array-valued CHEBFUN.
        f = assignColumns(f, colIdx, g)
        
        % Convert a CHEBFUN to another TECH.
        f = changeTech(f, newtech);

        % Deprecated function.
        f = define(f,s,v);
        
        % Supply a new definition for a CHEBFUN on a subinterval.
        f = defineInterval(f, subInt, g)
        
        % Supply new definition for a CHEBFUN at a point or set of points.
        f = definePoint(f, s, v)
        
        % Function used for deflation, wrapper for CHEBMATRIX/DEFLATIONFUN.
        out = deflationFun(Nu, u, r, p, alp, type)
        
        % Multiplication operator.
        M = diag(f)
        
        % Dimension check based on Matlab version
        out = dimCheck(f, g)

        % Useful information for DISPLAY.
        [name, data] = dispData(f)
        
        % Compare domains of two CHEBFUN objects.
        pass = domainCheck(f, g);        

        % Extract columns of an array-valued CHEBFUN object.
        f = extractColumns(f, columnIndex);
        
        % Deprecated function.
        varargin = fzero(varargout);
        
        % Get Delta functions within a CHEBFUN.
        [deltaMag, deltLoc] = getDeltaFunctions(f);
        
        % Get roots of a CHEBFUN and polish for use as breakpoints.        
        [rBreaks, rAll] = getRootsForBreaks(f, tol)
        
        % Returns true if numel(f) > 1
        out = isQuasi(f)
        
        % Number of columns (or rows) of a CHEBFUN quasimatrix.
        out = numColumns(f)
        
        % Solve a PDE with CHEBFUN.
        varargout = pdeSolve(varargin);
        
        % Obtain data used for plotting a CHEBFUN object:
        data = plotData(f, g, h)
        
        % Deprecated function.
        varargin = quad(varargout);

        % Sample a CHEBFUN on an "appropriate" grid.
        varargout = sample(f, n);
        
        % Reset pointValues property to the average of left and right limits.
        f = resetPointValues(f);
        
        % Remove all-zero layers of higher-order impulses.
        f = tidyImpulses(f)
        
        % Adjust nearby common break points in domains of CHEBFUN objects.
        [f, g, newBreaksLocF, newBreaksLocG] = tweakDomain(f, g, tol, pos)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = false )
        % Set small breakpoint values to zero.
        f = thresholdBreakpointValues(f);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Discrete cosine transform:
        y = dct(u, kind);
        
        % Discrete Legendre transform:
        y = dlt(u);
        
        % Discrete sine transform:
        y = dst(u, kind);
        
        % Inverse discrete cosine transform:
        u = idct(y, kind);
        
        % Inverse discrete Legendre transform:
        u = idlt(y, kind);
                
        % Inverse discrete sine transform:
        u = idst(y, kind);
        
        % Interpolate data:
        f = interp1(x, y, method, dom);
        
        % Inverse nonuniform fast Fourier transform: 
        [y, p] = inufft( varargin );

        % Compute Lagrange basis functions for a given set of points.
        f = lagrange(x, varargin);
        
        % Non-uniform discrete cosine transform:
        y = ndct(x, coeffs, th);
        
        % Non-uniform fast Fourier transform: 
        [y, p] = nufft( varargin );
        
        % Two-dimensional NUFFT: 
        f = nufft2( varargin ); 

        % ODE113 with CHEBFUN output.
        [t, y] = ode113(varargin);
        
        % ODE15S with CHEBFUN output.
        [t, y] = ode15s(varargin);
        
        % ODE45 with CHEBFUN output.
        [t, y] = ode45(varargin);

        % ODE78 with CHEBFUN output.
        [t, y] = ode78(varargin);

        % ODE89 with CHEBFUN output.
        [t, y] = ode89(varargin);        
        
        % Cubic Hermite interpolation:
        f = pchip(x, y, method);
        
        % Cubic spline interpolant:
        f = spline(x, y, d);
        
        % Update Chebfun source files:
        update(varargin)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Hidden = true, Static = true )

        % Convert a cell array of CHEBFUN objects to a quasimatrix.
        G = cell2quasi(F)
        
        % Determine values of CHEBFUN at breakpoints.
        vals = getValuesAtBreakpoints(funs, ends, op);
        
        % Which interval is a point in?
        out = whichInterval(dom, x, direction);

        % Like built-in UNION but with a tolerance for checking equality.
        C = tolUnion(A, B, tol);

        % Like built-in UNIQUE but with a tolerance for checking equality.
        C = tolUnique(A, tol);
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )

        % Main constructor.
        [funs, ends] = constructor(op, domain, data, pref);
        
        % Convert ODE solutions into CHEBFUN objects:
        [t, y] = odesol(sol, dom, opt);
        
        % Call one of the MATLAB ODE solvers and return a CHEBFUN
        [t, y] = constructODEsol(solver, odefun, tspan, uinit, varargin);
        
        % Parse inputs to PLOT. Extract 'lineWidth', etc.
        [lineStyle, pointStyle, jumpStyle, deltaStyle, out] = ...
            parsePlotStyle(varargin)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Class-related functions: private utilities for this m-file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op = str2op(op)
    % Convert string inputs to either numeric format or function_handles.
    sop = str2num(op); %#ok<ST2NM> % STR2DOUBLE doesn't support str2double('pi')
    if ( ~isempty(sop) )
        op = sop;
    else
        depVar = symvar(op);
        if ( numel(depVar) ~= 1 )
            error('CHEBFUN:CHEBFUN:str2op:indepvars', ...
                'Incorrect number of independent variables in string input.');
        end
        op = eval(['@(' depVar{:} ')', op]);
    end
end

function [op, dom, data, pref, flags] = parseInputs(op, varargin)
    
    % TODO: Should we 'data' structure to be passed to the constructor?
    % Currently, like in CHEBFUN/COMPOSE(), we don't have a use for this, but it
    % might be useful in the future.
    
    % Non-preferences that need to live beyond parseInputs.
    flags = struct();
    flags.done = false;         % No construction needs to take place.
    flags.doubleLength = false; % We will construct to double the length.
    flags.trunc = false;        % We will truncate the result.

    % Deal with string input options.
    if ( strncmp(op, '--', 2) )
        % An option has been passed to the constructor.
        if ( strcmpi(op, '--update') )
            chebfun.update();
        elseif ( strcmpi(op, '--update-devel') )
            chebfun.update('development');
        elseif ( strcmpi(op, '--version') )
            installDir = chebfunroot();
            fid = fopen(fullfile(installDir, 'Contents.m'), 'r');
            fgetl(fid);
            str = fgetl(fid);
            disp(['Chebfun ', str(3:end)]);
            fclose(fid);
        else
            error('CHEBFUN:parseInputs:unknown', ...
                'Unknow command %s.', op);
        end
        dom = [];
        data = struct();
        pref = [];
        flags.done = true;
        return
    end

    % Initialize data output.
    data.hscale = [];
    data.vscale = [];
    data.exponents = [];
    data.singType = [];
    
    args = varargin;

    % An op-only constructor call.
    if ( nargin == 1 )
        pref = chebfunpref();
    end

    % Try to parse out the domain which, if passed, is the second argument.
    domainWasPassed = false;
    if ( ~isempty(args) )
        if ( isnumeric(args{1}) && ...
                ((length(args{1}) >= 2) || isempty(args{1})) )
            dom = args{1};
            args(1) = [];
            domainWasPassed = true;
        elseif ( isa(args{1}, 'domain') )
            dom = double(args{1});
            args(1) = [];
            domainWasPassed = true;
        end
    end
    
    % A struct to hold any preferences supplied by keyword (name-value pair).
    keywordPrefs = struct();

    % Parse the remaining arguments.
    prefWasPassed = false;
    isPeriodic = false;
    vectorize = false;
    doVectorCheck = true;
    
    while ( ~isempty(args) )
        if ( isstruct(args{1}) || isa(args{1}, 'chebfunpref') )
            % Preference object input.  (Struct inputs not tied to a keyword
            % are interpreted as preference objects.)
            if ( ~prefWasPassed )
                pref = chebfunpref(args{1});
                prefWasPassed = true;
                args(1) = [];
            else
                error('CHEBFUN:CHEBFUN:parseInputs:twoPrefs', ...
                    'Multiple preference inputs are not allowed.');
            end
        elseif ( strcmpi(args{1}, 'equi') )
            % Enable FUNQUI when dealing with equispaced data.
            keywordPrefs.enableFunqui = true;
            args(1) = [];
        elseif ( strcmpi(args{1}, 'vectorize') || ...
                 strcmpi(args{1}, 'vectorise') )
            % Vectorize flag for function_handles.
            vectorize = true;
            args(1) = [];
        elseif ( strcmpi(args{1}, 'novectorcheck') )
            % Vector check for function_handles.
            doVectorCheck = false;
            args(1) = [];
        elseif ( strcmpi(args{1}, 'vectorcheck') )
            % Vector check for function_handles.
            doVectorCheck = strcmpi(args{2}, 'on');
            args(1:2) = [];            
        elseif ( strcmpi(args{1}, 'doublelength') )
            % Construct Chebfun twice as long as usually would be constructed.
            flags.doubleLength = true;
            args(1) = [];
        elseif ( strcmpi(args{1}, 'turbo') )
            % "turbo" flag for constructing "turbocharged" chebfuns.
            keywordPrefs.techPrefs.useTurbo = true;
            args(1) = [];
        elseif ( strcmpi(args{1}, 'coeffs') && isnumeric(op) )
            % Hack to support construction from coefficients.            
            op = {{[], op}};
            args(1) = [];
        elseif ( strcmpi(args{1}, 'coeffs') && iscell(op) )
            error('CHEBFUN:CHEBFUN:parseInputs:coeffcell', ...
                'Cannot construct CHEBFUN from a cell array of coefficients.');
        elseif ( any(strcmpi(args{1}, {'periodic', 'trig'})) )
            isPeriodic = true;
            args(1) = [];
        elseif ( strncmpi(args{1}, 'truncate', 5) )
            % Set the local truncation option.
            flags.trunc = args{2};
            args(1:2) = [];   
            % We split when truncation is selected. TODO: Why?
            keywordPrefs.splitting = true;
        elseif ( isnumeric(args{1}) && isscalar(args{1}) )
            % g = chebfun(@(x) f(x), N)
            keywordPrefs.techPrefs.fixedLength = args{1};
            args(1) = [];
        elseif ( strcmpi(args{1}, 'splitting') )
            keywordPrefs.splitting = strcmpi(args{2}, 'on');
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'minsamples') )
            % Translate "minsamples" --> "techPrefs.minSamples".
            keywordPrefs.techPrefs.minSamples = args{2};
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'blowup') )
            if ( strcmpi(args{2}, 'off') )
                % If 'blowup' is 'off'.
                keywordPrefs.blowup = 0;
            else
                % If 'blowup' is not 'off', set the singTypes.  (NB:  These
                % cells really need to store a left and right singType for each
                % domain subinterval, but we may not know the domain yet, so we
                % store just one cell for now and replicate it later, after
                % we've figured out the domain.)
                if ( (isnumeric(args{2}) && args{2} == 1 ) || ...
                        strcmpi(args{2}, 'on') )
                    % Translate "blowup" and flag "1" -->
                    % "blowup" and "poles only".
                    keywordPrefs.blowup = 1;
                    data.singType = {'pole'};
                elseif ( args{2} == 2 )
                    % Translate "blowup" and flag "2" -->
                    % "blowup" and "fractional singularity".
                    keywordPrefs.blowup = 1;
                    data.singType = {'sing'};
                else
                    error('CHEBFUN:CHEBFUN:parseInputs:badBlowupOption', ...
                        'Invalid value for ''blowup'' option.');
                end
            end
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'vscale') )
            % Store vscale types.
            data.vscale = args{2};
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'hscale') )
            % Store vscale types.
            data.hscale = args{2};
            args(1:2) = [];            
        elseif ( strcmpi(args{1}, 'singType') )
            % Store singularity types.
            data.singType = args{2};
            args(1:2) = [];            
        elseif ( strcmpi(args{1}, 'exps') || strcmpi(args{1}, 'exponents') )
            % Store exponents.
            data.exponents = args{2};
            args(1:2) = [];
        elseif ( any(strcmpi(args{1}, 'chebkind')) )
            % Translate "chebkind" and "kind" --> "tech.@chebtech".
            if ( (isnumeric(args{2}) && (args{2} == 1)) || ...
                     (ischar(args{2}) && strncmpi(args{2}, '1st', 1)) )
                keywordPrefs.tech = @chebtech1;
            elseif ( (isnumeric(args{2}) && (args{2} == 2)) || ...
                     (ischar(args{2}) && strncmpi(args{2}, '2nd', 1)) )
                keywordPrefs.tech = @chebtech2;
            else
                error('CHEBFUN:CHEBFUN:parseInputs:badChebkind', ...
                    'Invalid value for ''chebkind'' option.');
            end
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'resampling') )
            % Translate "resampling" --> "techPrefs.refinementFunction".
            if ( strcmpi(args{2}, 'on') )
                keywordPrefs.techPrefs.refinementFunction = 'resampling';
            elseif ( strcmpi(args{2}, 'off') )
                keywordPrefs.techPrefs.refinementFunction = 'nested';
            end
            args(1:2) = [];
        elseif ( any(strcmpi(args{1}, 'eps')) )
            % Translate "eps" --> "techPrefs.chebfuneps".
            keywordPrefs.techPrefs.chebfuneps = args{2};
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'maxdegree') )
            % Translate "maxdegree" --> "techPrefs.maxLength".
            keywordPrefs.techPrefs.maxLength = args{2};
            args(1:2) = [];
        elseif ( any(strcmpi(args{1}, {'splitLength', 'splitdegree'})) )
            % Translate "splitdegree" --> "splitPrefs.splitLength".
            keywordPrefs.splitPrefs.splitLength = args{2};
            args(1:2) = [];
        elseif ( strcmpi(args{1}, 'splitMaxLength') )
            % Translate "splitMaxLength" --> "splitPrefs.splitMaxLength".
            keywordPrefs.splitPrefs.splitMaxLength = args{2};
            args(1:2) = [];
        elseif ( ischar(args{1}) )
            % Update these preferences:
            if ( length(args) < 2 )
                error('CHEBFUN:CHEBFUN:parseInputs:noPrefValue', ...
                    ['Value for ''' args{1} ''' preference was not supplied.']);
            end
            keywordPrefs.(args{1}) = args{2};
            args(1:2) = [];
        else
            if ( isnumeric(args{1}) )
                error('CHEBFUN:CHEBFUN:parseInputs:badInputNumeric', ...
                    ['Could not parse input argument sequence.\n' ...
                     '(Perhaps the construction domain is not the second ' ...
                     'argument?)']);
            else
                error('CHEBFUN:CHEBFUN:parseInputs:badInput', ...
                    'Could not parse input argument sequence.');
            end
        end
    end
    
    % Construction from equispaced data requires the number of points to be
    % specified
    if ( ~isnumeric(op) && isfield(keywordPrefs, 'enableFunqui') && ...
            (~isfield(keywordPrefs, 'techPrefs') || ...
            (isfield(keywordPrefs, 'techPrefs') && ...
            ~isfield(keywordPrefs.techPrefs,'fixedLength'))) )
        error('CHEBFUN:CHEBFUN:parseInputs:equi', ...
            '''equi'' flag requires the number of points to be specified.');
    end
    
    % It doesn't make sense to construct from values and coeffs at the same
    % time.
    if ( iscell(op) && iscell(op{1}) && isfield(keywordPrefs, 'tech') && ...
            ~isempty(keywordPrefs.tech) )
        error('CHEBFUN:CHEBFUN:parseInputs:coeffschebkind', ...
            [' ''coeffs'' and ''chebkind'' should not be ' ...
            'specified simultaneously.']);
    end
    
    % Override preferences supplied via a preference object with those supplied
    % via keyword.
    if ( prefWasPassed )
        pref = chebfunpref(pref, keywordPrefs);
    else
        pref = chebfunpref(keywordPrefs);
    end

    % Use the domain of the chebfun that was passed if none was supplied.
    if ( ~domainWasPassed || isempty(dom) )
        if ( isa(op, 'chebfun') )
            dom = [ op.domain(1) op.domain(end) ];
        else
            dom = pref.domain;
        end
    end
    numIntervals = numel(dom) - 1;

    % Error if 'doubleLength' and 'splitting on' are both passed:
    % This combination is not supported.
    if ( pref.splitting && flags.doubleLength )
        error('CHEBFUN:CHEBFUN:parseInputs:doubleLengthSplitting', ...
            'doubleLength not supported with splitting on.')
    end

    % Error if 'doubleLength' is used on a domain with breakpoints:
    % This combination is not supported.
    if ( (length(dom) > 2) && flags.doubleLength )
        error('CHEBFUN:CHEBFUN:parseInputs:doubleLengthBreakpoints', ...
            'doubleLength not supported on domains with breakpoints.')
    end

    % Deal with the 'periodic' or 'trig' flag:
    if ( isPeriodic )
        % Translate 'periodic' or 'trig'.
        pref.tech = @trigtech;
        pref.splitting = false;
        pref.enableFunqui = false;
        if ( numel(dom) > 2 )
            error('CHEBFUN:parseInputs:periodic', ...
                '''periodic'' or ''trig'' option is only supported for smooth domains.');
        end
    end

    % Parse the OP (handle the vectorize flag, etc.).
    if ( iscell(op) )
        for k = 1:numel(op)
            op{k} = parseOp(op{k});
        end
    else
        op = parseOp(op);
    end

    function op = parseOp(op)
        % Convert string input to function_handle:
        if ( ischar(op) )
            op = str2op(op);
        end
        if ( doVectorCheck && isa(op, 'function_handle') )
            op = vectorCheck(op, dom, vectorize);
        end
        if ( isa(op, 'chebfun') )
            if ( op.isTransposed )
                op = op';
            end
            op = @(x) feval(op, x);
        end
        if ( isa(op, 'function_handle') && pref.enableFunqui )
            if ( isfield(pref.techPrefs, 'fixedLength') && ...
                 ~isnan(pref.techPrefs.fixedLength) )
                x = linspace(dom(1), dom(end), pref.techPrefs.fixedLength).';
                op = feval(op, x);
                pref.techPrefs.fixedLength = NaN;
            end
        end
    end

    % Enable singularity detection if we have exponents or singTypes:
    if ( any(data.exponents) || ~isempty(data.singType) )
        pref.blowup = true;
    end
    % Sort out the singularity types:
    if ( numel(data.singType) == 1 )
        singType = data.singType{1};
        data.singType = cell(1, 2*numIntervals);
        for j = 1:2*numIntervals
            data.singType{j} = singType;
        end
    elseif ( ~isempty(data.singType) && ...
            (numel(data.singType) ~= 2*numIntervals) )
        % If the number of exponents supplied by user isn't equal to twice the
        % the number of the FUNs, throw an error message:
        error('CHEBFUN:CHEBFUN:parseInputs:badExponents', ...
            'The number of the exponents is inappropriate.');
    end
    % Sort out the exponents:
    if ( ~isempty(data.exponents) )
        exps = data.exponents;
        nExps = numel(exps);
        if ( nExps == 1 )
            % If only one exponent is supplied, assume the exponent at other
            % breakpoints are exactly same.
            exps = exps*ones(1, 2*numIntervals);
        elseif ( nExps == 2 )
            % If the exponents are only supplied at endpoints of the entire
            % domain, then pad zeros at the interior breakpoints.
            exps = [exps(1) zeros(1, 2*(numIntervals-1)) exps(2)];
        elseif ( nExps == numIntervals + 1 )
            % If only one exponent is supplied for each interior breakpoint,
            % then we assume that the singularity take the same order on each
            % side.
            exps = exps(ceil(1:0.5:nExps - 0.5));
        elseif( nExps ~= 2*numIntervals )
            % The number of exponents supplied by user makes no sense.
            error('CHEBFUN:CHEBFUN:chebfun:parseInputs', ...
                'Invalid length for vector of exponents.');
        end
        data.exponents = exps;
    end
    
    % Ensure DOM is a double (i.e., not a domain object).
    dom = double(dom);

end

function op = vectorCheck(op, dom, vectorize)
%VECTORCHECK   Try to determine whether op is vectorized. 
%   It's impossible to cover all eventualities without being too expensive. 
%   We do the best we can.

y = dom([1 end]); y = y(:);

% Make a slightly narrower domain to evaluate on. (Endpoints can be tricky).
% This used to be fixed at 0.01. But this can cause troubles at very narrow
% domains, where 1.01*y(1) might actually be larger than y(end)!
y = y + [1;-1].*diff(y)/200;

if ( vectorize )
    op = vec(op, y(1));
end

try
    % Evaluate a vector of (near the) endpoints
    v = op(y);
    
    % Get the size of the output:
    sv = size(v);
    sy = size(y);
    
    % Check the sizes:
    if ( sv(1) == sy(1) )
        % Here things seem OK! 
        
        % However, we may possibly be fooled if we have an array-valued function
        % whose number of columns equals the number of test points(i.e., 2) or 
        % something unvectorized like sin(x)/x (no dot). We choose one 
        % additional point as a final check:
        if ( sv(2) == sy(1) )
            v = op(y(1));
            if ( size(v, 1) > 1 )
                op = @(x) op(x).';
                warning('CHEBFUN:CHEBFUN:vectorCheck:transpose', ...
                    ['Chebfun input should return a COLUMN array.\n', ...
                     'Attempting to transpose.'])
            elseif ( size(v, 2) ~= sv(2) )
                % It doesn't really matter what this error message is as it will
                % be caught in the try-catch statement.
                error('CHEBFUN:CHEBFUN:vectorCheck:numColumns', ...
                    'Number of columns increases with length(x).');
                
            end
                
        end
        
    elseif ( all( sv == 1 ) )
        % The operator always returns a scalar:
        op = @(x) repmat(op(x), length(x), 1);
        
    elseif ( any(sv == sy(1)) )
        
        if ( any(sv) == 1 )
            % We check to see if we have something like @(x) [1 1].
            v = op(y(1)); % Should evaluate to a scalar, unless array-valued.
            if ( all(size(v) == sv) )
                % Scalar expand:
                op = @(x) repmat(op(x), length(x), 1);
                return
            end
        end
        
        if (size(op(y), 1) > 1)
            % Try and transpose:
            op = @(x) op(x).';
            warning('CHEBFUN:CHEBFUN:vectorCheck:transpose', ...
                    ['Chebfun input should return a COLUMN array.\n', ...
                    'Attempting to transpose.'])
        else
            % Try and vectorize:
            op = vectorCheck(op, dom, 1);
        end

    elseif ( any(sv == 1) )
        % The operator always returns a scalar:
        op = @(x) repmat(op(x), length(x), 1);
        
    end

catch ME
    % The above didn't work. :(
    
    if ( vectorize )
        % We've already tried vectorizing, so we've failed.
        rethrow(ME)
        
    else
        % Try vectorizing. (This is now done silently.)
        op = vectorCheck(op, dom, 1);
    end
    
end

end

function g = vec(op, y)
%VEC  Vectorize a function or string expression.
%   VEC(OP, Y), if OP is a function handle or anonymous function, returns a
%   function that returns vector outputs for vector inputs by wrapping F inside
%   a loop. Y, serving as a testing point, is a point in the domain where OP is
%   defined and is used to determine if OP is array-valued or not.
    
    % Check to see if OP is array-valued:
    opy = op(y);
    if ( any(size(opy) > 1) )
        % Use the array-valued wrapper:
        g = @loopWrapperArray;    
    else
        % It's not array-valued. Use the scalar wrapper:
        g = @loopWrapperScalar;
    end
    
    % Nested functions:
    
    % Scalar case:
    function v = loopWrapperScalar(x)
        v = zeros(size(x));
        for j = 1:numel(v)
            v(j) = op(x(j));
        end
    end

    % Array-valued case:
    function v = loopWrapperArray(x)
        numCol = size(op(x(1)), 2);
        numRow = size(x, 1);
        v = zeros(numRow, numCol);
        for j = 1:numRow
            v(j,:) = op(x(j));
        end
    end

end
