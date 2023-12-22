clear;
close all;
clc;

%% Test Algorithm: Rakovic, Invariant Approximations of the Minimal Robust Positively Invariant Set
context = Context();
specifics.sys.A = [1 1; 0 1];
specifics.sys.B = [1 1]';
K = -[1.17 1.03];
specifics.sys.Ak = specifics.sys.A + specifics.sys.B*K;
center = [0; 0];
generators = diag([1 1]);
specifics.unc.zonotope.W = conZono(center,generators);
specifics.unc.polytope.W = Polyhedron('lb',[-1 -1],'ub',[1 1]);

gen = processUncertainty.RobustSetGenerator(context, specifics);

figure;
hold on;
eps = 5e-5;
tic
[Finf, s_inf, alfa] = gen.compute_minRPI(eps);
Tpoly = toc;
plot(Finf)
for s = 10:-1:1
    Fs = gen.incomplete_minRPI(s);
    plot(Fs)
end
annotation('textarrow',[0.45 0.4],[0.5 0.42],'String','s=1','FontSize',12);
annotation('textarrow',[0.7 0.6],[0.75 0.82],'String','s=2','FontSize',12);
n_steps = 5;
xlabel("$x_1$"); ylabel("$x_2$")
tic;
Finf_zono = gen.compute_zonotope_minRPI(n_steps);
Tzone = toc;