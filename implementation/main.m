clear;
close all;
clc;
%% Plot Settings
set(0,'defaultLineLineWidth', 1.5)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Define parameters for context
%rng("shuffle")
rng(1)
addpath(genpath('third-party'));
context = Context();

% ----- PARAMETERS FOR DYNAMIC SIMULATION ----- %
parameters.J_T = [17023.3 397.1 -2171.4;
                  397.1 124825.7 344.2;
                  -2171.4 344.2 129112.2]; % kg m^2
parameters.J_C = diag([362, 604, 560]);
%parameters.m_T = 7827.867; % kg
m_C = 850;
parameters.m_C = m_C; % kg
parameters.OM = 1.0454e-3; % rad/s from (a=7144.8km)
parameters.pBerth = [0 0 -5.5]'; % m -> from sum pieces
%parameters.pBerth = [4 0 0]'; % m -> from sum pieces

parameters.u_lim = 100*ones(3,1);
% Chaser dimensions
% x = 2.38m , y = 1.5m, z = 1.695m
% Tumbler dimension
% Main: x = 10 2.75 2.075
% Boom: 3m -> connect panel to main body
% Solar panel section 1.3333 4.97

% Ellipsoid keep-out zone definition
parameters.ellipse.r = [17 8 6];
parameters.ellipse.P = diag(1./parameters.ellipse.r.^2);
parameters.ellipse.xc = [1.5 0 0.75]';

context.parameters = parameters;
%% Define Reference and MPC settings for context
initialChaserState = [40 0 0 0 0 0]';
context.ref.initialChaserState = initialChaserState;
context.real.initialChaserState = initialChaserState;
maneuverTime = 300;
context.ref.maneuverTime = maneuverTime;
context.mpc.N = 20;
context.mpc.T = 0.5;
context.mpc.R = 1*eye(3);
%context.mpc.Q = diag([1e5 1e-4 1e-4 1e7 1e8 1e8]); % true
context.mpc.Q = diag([1e4 1e4 1e4 1e-3 1e-3 1e-3]);
Ttube = 10; context.ref.Ttube = Ttube;

% Define uncertainty generation
context.unc.nTrials = 500;
context.unc.flagDynamics = false;
context.unc.dJ = 100*ones(3);

context.real.J_T = parameters.J_T + context.unc.dJ.*utils.randnsym(3);
context.disturbance.PSI = eye(2);

% context.unc.dalfa = 0.0001745329;
% context.unc.dtheta = 0.0001745329;
context.unc.baseOmega = 0.0005*diag([2,1,2]);
context.unc.baseQ = 0.001*diag([1,1,1,1]);
% context.unc.baseOmega = 0.005*randn(3);
% context.unc.baseQ = 0.01*randn(4);


%% Define mpc variables for context
% A = zeros(6);
% A(1:3,4:6) = eye(3);
[~,~,m_C,OM] = dynamics.set_parameters(parameters);
A = zeros(6);
A(1:3,4:6) = eye(3);
A(4:6,1:3) = [3*OM^2 0 0;  0 0 0; 0 0 -OM^2];
A(4:6,4:6) = [0 2*OM 0; -2*OM 0 0; 0 0 0];

B = zeros(6,3);
B(4:6,1:3) = 1/m_C*eye(3);

C = eye(6);
D = zeros(6,3);

sys = ss(A,B,C,D);
T = context.mpc.T;
sysd = c2d(sys,T);
A = sysd.A;
B = sysd.B;
context.mpc.sysd = sysd;

Q = 1e3*eye(6);
R = 1*eye(3);
% Q = blkdiag(eye(3),inv(parameters.J_T));
% R = 1e-3*eye(3);
[K,P,e] = dlqr(sysd.A,sysd.B,Q,R);
K = -K;
% context.mpc.Q = Q;
% context.mpc.R = R;
context.mpc.P = P;
context.mpc.K = K;
Ak = context.mpc.sysd.A + context.mpc.sysd.B*context.mpc.K;
AkX = Ak([1 4],[1 4]);
AkY = Ak([2 5],[2 5]);
AkZ = Ak([3 6],[3 6]);

% Check behaviour of lqr
% sysClose = ss(Ak,B,C,D,T);
% step(sysClose)
%% Instantiate Reference Generator

%Instantiate reference generator (with best knowledge at t=0)
qT_init =  [0.3826834, 0, 0, 0.9238795 ]'; % 45° rotation around x axis
wT_init = [0.0174533 0.0349066 0.0174533]'; % 2deg/s on y and z 0.0174533 0.0349066

base = context.unc.baseQ;
covQ = base.'*base;
qT_real = mvnrnd(qT_init,covQ,1)'; % qT_real = [0.382814540459003 9.20233639674008e-05 -0.00191834239742808 0.925496145793015]
base = context.unc.baseOmega;
covOmega = base.'*base;
wT_real = mvnrnd(wT_init,covOmega,1)'; % wT_real = [0.0187137391698965	0.0345187856757413	0.0181421389368183]

initialTumblerState = [wT_init; qT_init];
context.ref.initialTumblerState = initialTumblerState;

realOrientation.tumblerState = [wT_real; qT_real];
context.real.initialTumblerState = [wT_real; qT_real];


trajectoryGenerator = motionPlanning.ReferenceGenerator(context);
%%
solutionRealInertial = visualizeRealTraj(trajectoryGenerator.solutionTumblerFrame,context);
f = figure(1);
hold on;
visualization.plotSolutionTrajectoryGeneration(trajectoryGenerator)
figure(f);
plot3(solutionRealInertial.xRefI(1,:),solutionRealInertial.xRefI(2,:),solutionRealInertial.xRefI(3,:))
%% Assign bounds
t = 0;
simulator = Simulator(context);
% 0) Do dynamic step
outputEstimator = simulator.reset(); % first step -> no control

% 1) Generate reference 0->N
% outputRefGenerator = trajectoryGenerator.generateSafeReferenceTumblerFrame(t);
% outputRefGeneratorIn = trajectoryGenerator.generateDebugReferenceInertialFrame(t);
% outputRefGenerator = tumblerToInertial(outputRefGenerator, outputEstimator);
outputRefGenerator = trajectoryGenerator.generateForBound(t);

% 2) Generate bounds
mode = "sigmasigma";
boundGenerator = processUncertainty.MonteCarloBound(context);
outputBoundGenerator = boundGenerator.disturbanceBound(outputRefGenerator, outputEstimator, mode);
outputBoundGenerator = max(outputBoundGenerator,[],2);

%% Prepare Robust set generators

% 3) Generate RPI
rpiGens = [];
Xbound = [100 100 100 5 5 5]';
Ubound = 100*ones(3,1);
selector = [1 4; 2 5; 3 6];
specifics = struct([]);
for i = 1:3
    sel = selector(i,:);
    specifics(i).sys.Ak = Ak(sel,sel);
    specifics(i).sys.A = A(sel,sel);
    specifics(i).sys.B = B(sel,i);
    specifics(i).sys.K = K(i,sel);
    specifics(i).sys.Q = Q(sel,sel);
    specifics(i).sys.R = R(i,i);
    specifics(i).sys.P = P(sel,sel);
    specifics(i).sys.C = [eye(2);zeros(1,2)];
    specifics(i).sys.D = [zeros(2,1);1];
    specifics(i).unc.zonotope.W = conZono(zeros(2,1),diag(outputBoundGenerator(sel)));
    specifics(i).unc.polytope.W = Polyhedron('lb',-outputBoundGenerator(sel),'ub',outputBoundGenerator(sel));
    specifics(i).Xbound.zonotope.W = conZono(zeros(2,1),diag(Xbound(sel)));
    specifics(i).Xbound.polytope.W = Polyhedron('lb',-Xbound(sel),'ub',Xbound(sel));
    specifics(i).Ubound.zonotope.W = conZono(0,diag(Ubound(i)));
    specifics(i).Ubound.polytope.W = Polyhedron('lb',-Ubound(i),'ub',Ubound(i));
    specifics(i).Ybound.polytope.Y = Polyhedron('lb',[-Xbound(sel); -Ubound(i)],'ub',[Xbound(sel); Ubound(i)]);
    rpiGens = [rpiGens, processUncertainty.RobustSetGenerator(context, specifics(i))];
end
%%
%[outputRefVec,outputStateEstVec] = utils.demux(outputRefGenerator, outputEstimator);
tic;
for i = 1:3
    outputRobustSet{i} = rpiGens(i).tightenConstraint();
    %[uLaw,sol{i}] = ctr(i).solve(outputRefVec(i), outputStateEstVec(i), outputRobustSet);
end
T = toc;
%% Plot sets
ylabels = ["$\dot{x} [m/s]$","$\dot{y} [m/s]$","$\dot{z}$ [m/s]"];
xlabels = ["x [m]","y [m]","z [m]"];
for i = 1:3
    figure;
    plot(outputRobustSet{i}.MaxRPI,'color',[1,0,0],'alpha',0.5)
    hold on;
    plot(outputRobustSet{i}.minRPI,'color',[0.7 0.7 0.7],'alpha',0.5)
    legend(["MRPI","mRPI"])
    xlabel(xlabels(i)); ylabel(ylabels(i));
    figure;
    plot(outputRobustSet{i}.minRPI,'color',[0.7 0.7 0.7],'alpha',0.5)
    hold on
    xlabel(xlabels(i)); ylabel(ylabels(i));
%     for j = 10:-1:1
%         plot(rpiGens(i).incomplete_minRPI(j))
%     end
    %plot(outputRobustSet{i}.minRPI)
end
%%
upState = [];
upAction = [];
for i = 1:3
    Xtight = outputRobustSet{i}.Xtight;
    Vtight = outputRobustSet{i}.Vtight;
    % Obtain lower and upper limits on elements
    Xt = abs(Xtight.V(1,:)');
    Vt = abs(Vtight.V(1,:)');
    upState = [upState, Xt];
    upAction = [upAction; Vt];
end
upState = [upState(1,:), upState(2,:)]'; lwState = -upState; lwAction = -upAction;
tightConstraints.upState = upState;
tightConstraints.upAction = upAction;
%% Casadi with class
import casadi.*

robustCon = tubeMPCcontrol.RobustCasadi(context, outputRobustSet, tightConstraints);
t = 0;
simulator = Simulator(context);
dtEst = 0.1;
outputEst = simulator.reset();
N = context.mpc.N; nStates = 6; nActions = 3;
Tcomp = 3; % 3s to compute bounds and RPI
TtubeNow = 0;

actionRecord = [];
mpcStateRecord = [];
mpcActionRecord = [];
stateRecord = [];
TcompTube = [];
Tcomp = [];
while true
    outputRefGenerator = trajectoryGenerator.generateSafeReferenceTumblerFrame(t);
    if isempty(outputRefGenerator.time)
        fprintf("Terminated")
        break;
    end
    outputRefGenerator = tumblerToInertial(outputRefGenerator, outputEst);
    if TtubeNow + Tcomp >= Ttube
        fprintf("Computing new tube at t=%f\n",t)
        tic;
        outputRefBound = trajectoryGenerator.generateForBound(t);
        outputBoundGenerator = boundGenerator.disturbanceBound(outputRefBound, outputEst, mode); % bonus
        outputBoundGenerator = max(outputBoundGenerator,[],2);
        T1 = toc;
        tic;
        [outputRobustSet, tightConstraints] = genRobustSet(Ak,A,B,K,Q,R,P,outputBoundGenerator,context); % bonus
        T2 = toc;
        robustCon = tubeMPCcontrol.RobustCasadi(context, outputRobustSet, tightConstraints); % bonus
        TtubeNow = 0;
        TcompTube = [TcompTube, [T1; T2]];
    end
    if size(outputRefGenerator.xRefI,2) <= context.mpc.N
        context.mpc.N = context.mpc.N - 1;
        N = context.mpc.N;
        robustCon = tubeMPCcontrol.RobustCasadi(context, outputRobustSet, tightConstraints); % bonus
    end

    args.x0 = [reshape(outputRefGenerator.xRefI,[(N+1)*nStates 1]); reshape(outputRefGenerator.vRefI(:,1:N),[N*nActions 1])];
    args.p = [outputEst.chaserState; reshape(outputRefGenerator.xRefI,[(N+1)*nStates 1]); reshape(outputRefGenerator.vRefI(:,1:N),[N*nActions 1])];

    tic;
    [uLaw, xsol, usol] = robustCon.solve(outputRefGenerator, outputEst);
    Tcomp = [Tcomp, toc];

    for i = 1:context.mpc.T/dtEst
        outputEst = simulator.stepUncertain(dtEst,uLaw);
    end
    mpcStateRecord = [mpcStateRecord, xsol(:,1)];
    mpcActionRecord = [mpcActionRecord, usol(:,1)];
    actionRecord = [actionRecord, outputEst.uAction];
    stateRecord = [stateRecord, outputEst.chaserState];
    xpoint = outputEst.chaserState;
    t = t + context.mpc.T;
    TtubeNow = TtubeNow + context.mpc.T;
    figure(f)
    plot3(xpoint(1),xpoint(2),xpoint(3),'k.','MarkerSize',5)
    pause(0.0000000001);
    % [ts,us] = stairs(Info.Topt,Info.Uopt);
    % plot(ts,us,"r-",Info.Topt,Info.Yopt,"b--")
    % legend("MV","OV")
    if t >= context.ref.maneuverTime
        break
    end
end
% figure;
% stairs(1:length(actionRecord),actionRecord')
%% Plot the solution
figure; hold on;
tt = 0:0.5:299.5;
plot(tt,stateRecord(1:3,:),'LineWidth',2)
tt = linspace(0,300,1000);
plot(tt,solutionRealInertial.xRefI(1:3,:),'--')
plot(trajectoryGenerator.solutionInertialFrame.t,trajectoryGenerator.solutionInertialFrame.z(1:3,:),'-.')
legend(["x_{mpc}","y_{mpc}","z_{mpc}","x_{real}","y_{real}","z_{real}","x_{nom}","y_{nom}","z_{nom}"])
xlabel("Time [s]"); ylabel("Position [m]"); grid on;
figure; hold on;
tt = 0:0.5:299.5;
plot(tt,stateRecord(4:6,:))
tt = linspace(0,300,1000);
plot(tt,solutionRealInertial.xRefI(4:6,:),'--')
legend(["x_{mpc}","y_{mpc}","z_{mpc}","x_{real}","y_{real}","z_{real}"])
xlabel("Time [s]"); ylabel("Velocity [m/s]"); grid on;
%% Plot the error with respect to the real trajectory
figure;
tt = linspace(0,300,1000);
tt2 = 0:0.5:299.5;
s = interp1(tt,solutionRealInertial.xRefI(1:3,:).',tt2).';
plot(tt2, vecnorm(s - stateRecord(1:3,:),2,1))
%% Plot of the actions
figure; hold on;
tt = 0:0.5:299.5;
plot(tt,mpcActionRecord)
figure; hold on;
tt = 0:0.5:299.5;
plot(tt,abs(actionRecord - mpcActionRecord))
%% Help Functions

function solutionRealInertial = visualizeRealTraj(solutionTumblerFrame,context)

    tt = solutionTumblerFrame.t;
    zRefT = solutionTumblerFrame.z;
    vRefT = solutionTumblerFrame.v;
    x0 = context.real.initialTumblerState;
    param_real = context.realParameters();

    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,stateRotTrue] = ode45(@(t,x) dynamics.odetumbler(t,x,param_real), ...
                    tt,x0, options);
    stateRotTrue = stateRotTrue';

    for i = 1:length(tt)
        q = stateRotTrue(4:7,i);
        zRefI(1:3,i) = quat.rotate(zRefT(1:3,i),q);
        zRefI(4:6,i) = quat.rotate(zRefT(4:6,i),q);
        vRefI(:,i) = quat.rotate(vRefT(:,i),q);
    end

    solutionRealInertial.xRefI = zRefI;
    solutionRealInertial.vRefI = vRefI;


end

function [outputRobustSet, tightConstraints] = genRobustSet(Ak,A,B,K,Q,R,P,outputBoundGenerator,context)

% 3) Generate RPI
rpiGens = [];
Xbound = [100 100 100 5 5 5]';
Ubound = 100*ones(3,1);
selector = [1 4; 2 5; 3 6];
specifics = struct([]);
for i = 1:3
    sel = selector(i,:);
    specifics(i).sys.Ak = Ak(sel,sel);
    specifics(i).sys.A = A(sel,sel);
    specifics(i).sys.B = B(sel,i);
    specifics(i).sys.K = K(i,sel);
    specifics(i).sys.Q = Q(sel,sel);
    specifics(i).sys.R = R(i,i);
    specifics(i).sys.P = P(sel,sel);
    specifics(i).sys.C = [eye(2);zeros(1,2)];
    specifics(i).sys.D = [zeros(2,1);1];
    specifics(i).unc.zonotope.W = conZono(zeros(2,1),diag(outputBoundGenerator(sel)));
    specifics(i).unc.polytope.W = Polyhedron('lb',-outputBoundGenerator(sel),'ub',outputBoundGenerator(sel));
    specifics(i).Xbound.zonotope.W = conZono(zeros(2,1),diag(Xbound(sel)));
    specifics(i).Xbound.polytope.W = Polyhedron('lb',-Xbound(sel),'ub',Xbound(sel));
    specifics(i).Ubound.zonotope.W = conZono(0,diag(Ubound(i)));
    specifics(i).Ubound.polytope.W = Polyhedron('lb',-Ubound(i),'ub',Ubound(i));
    specifics(i).Ybound.polytope.Y = Polyhedron('lb',[-Xbound(sel); -Ubound(i)],'ub',[Xbound(sel); Ubound(i)]);
    rpiGens = [rpiGens, processUncertainty.RobustSetGenerator(context, specifics(i))];
end
for i = 1:3
    outputRobustSet{i} = rpiGens(i).tightenConstraint();
    %[uLaw,sol{i}] = ctr(i).solve(outputRefVec(i), outputStateEstVec(i), outputRobustSet);
end
upState = [];
upAction = [];
for i = 1:3
    Xtight = outputRobustSet{i}.Xtight;
    Vtight = outputRobustSet{i}.Vtight;
    % Obtain lower and upper limits on elements
    Xt = abs(Xtight.V(1,:)');
    Vt = abs(Vtight.V(1,:)');
    upState = [upState, Xt];
    upAction = [upAction; Vt];
end
upState = [upState(1,:), upState(2,:)]'; lwState = -upState; lwAction = -upAction;
tightConstraints.upState = upState;
tightConstraints.upAction = upAction;
end