clear;
close all;
clc;

%% Test Algorithm: Kolmanovsky and Gilbert, Theory and Computation of Disturbance Invariant Sets for Discrete-Time Linear Systems

% x+ = A*x + B*u; y = C*x + D*u
m = 0.5; M = 0.5; l = 1.4; g = 10;
A = [0 1 0 0; 0 0 -m*g/M 0; 0 0 0 1; 0 0 (M+m)*g/(M*l) 0];
B = [0 1/M 0 -1/(M*l)]';
C = eye(size(A));
D = zeros(length(A),size(B,2));
% discretize
T = 0.1;
sysc = ss(A,B,C,D);
sysd = c2d(sysc,T);
A = sysd.A;
B = sysd.B;

context = Context();
context.disturbance.PHI = [0 0 1 0]';

specifics.sys.A = A;
specifics.sys.B = B;
specifics.sys.C = zeros(1,length(A)); % y = u
specifics.sys.D = eye(1); % y = u
specifics.sys.K = [0.5451 1.8357 27.2815 8.6552];
specifics.sys.Ak = A + B*specifics.sys.K;

% set polytope for Y (||u|| < 1/2)
Y = Polyhedron([1; -1],[1/2; 1/2]);
specifics.Ybound.polytope.Y = Y;
Msize = 2;

% Set polytope for disturbance
dd = [0 0.001 0.003 0.00475 0.0048];
f1 = figure;
hold on;
f2 = figure;
hold on;


for i = 1:length(dd)
    d = dd(i);
    W = Polyhedron('lb',-d,'ub',d);
    specifics.unc.polytope.W = W;
    gen = processUncertainty.RobustSetGenerator(context,specifics);
    
    tic;
    [Oinf, t, Nt] = gen.compute_MaxRPI();
    Teval = toc;
    fprintf("**** Disturbance magnitude: %.4f **** Time evaluation: %.2f ****\n",dd(i),Teval)
    if ~isEmptySet(Oinf)
        P1 = slice(Oinf,[3 4],[0 0]); % fix dimension 3 and 4 at 0
        figure(f1)
        xlabel("$x_1$ [m]"); ylabel("$x_2$ [m/s]")
        plot(P1)
        P2 = slice(Oinf,[1 2],[0 0]); % fix dimension 1 and 2 at 0
        figure(f2)
        xlabel("$x_3$ [rad]"); ylabel("$x_4$ [rad/s]")
        plot(P2)
    end
end