clear;
close all;
clc;

%%
wBounds = combinations({[-1, 1], [-2 2], [-5 5]});

%%

A = [1 0.2 -1; 0 1 -0.2; 0 0 0.6];
B = [0 0 0.6]';
sysd.A = A;
sysd.B = B;
Q = diag([1 1 0.1]);
R = 0.1;
[K,~] = dlqr(A,B,Q,R);
K = -K;

X.cX = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
X.d = ones(3,1);
V.cV = [0 0 0 1]';
V.d = 1;

N = 10;
wlim = 0.001*[5, 5, 5]';
[Xtight,Vtight] = tightenConstraint(sysd,K,X,V,wlim,N);

%%

A = [1 1; 0 1];
B = [0; 1];
sysd.A = A;
sysd.B = B;

K = [-0.4 -1.2];
Ak = A + B*K;
nState = 2;

wlim = [-0.1, 0.1]; % limit of an edge
wBounds = combinations({wlim, wlim}); % find all vertexes
wBounds = wBounds';
wLow = [-0.1; -0.1];
wHigh = [0.1; 0.1];
fun = @(x) norm(x,Inf);

% convert box bounds in form:
% c'*z <= d 
% c'_x*x + c'_u*u <= d
cX = [eye(2); 0, 0]; % I think we can ignore these if they are just eye
cU = [0 0 1]';
d = ones(3,1);

X.cX = cX;
X.d = ones(2,1);
V.cV = cU;
V.d = 1;

% probably instead use always the same vertex, trying various sequence of
% vertexes lead to a too big number of computations
nDist = 2;
nCst = 3;
Nmax = 26;
alfa = zeros(1,Nmax);
thN = zeros(nCst,Nmax);
xLim = repmat(d,[1, Nmax]);
f = zeros(nCst,Nmax*nDist);
for N = 0:Nmax-1
    alfa(N+1) = max([maxOverSet(fun,Ak^N*wBounds)/maxOverSet(fun,wBounds), maxOverSet(fun,K*Ak^N*wBounds)/maxOverSet(fun,K*wBounds)]);
    f(:,nDist*N+1:nDist*(N+1)) = cX*Ak^N + cU*K*Ak^N;
%     th = -Inf*ones(3,1);
%     for vtx = 1:size(wBounds,2)
%         th_tmp = zeros(3,1);
%         for i = 0:N-1
%             th_tmp = th_tmp + [Ak^i*wBounds(:,vtx);K*Ak^i*wBounds(:,vtx)];
%         end
%         th = max(th,th_tmp);
%     end
    %xLim(:,N+1) = d - 1/(1 - alfa(N+1))*th;
    %thN(:,N+1) = th;
end
opts = optimoptions(@linprog,"Display","off");
Z1 = zeros(nCst,Nmax);
for N = 0:Nmax-1
    if alfa(N+1) < 1
        llim = repmat(wLow,[N, 1]);
        ulim = repmat(wHigh,[N, 1]);
        for s = 1:nCst
            c = -f(s,1:nDist*N);
            [w,th] = linprog(c,[],[],[],[],llim, ulim, opts);
            th = -th;
            thN(s,N+1) = th;
            xLim(s,N+1) = xLim(s,N+1) - 1/(1 - alfa(N+1))*th;
        end
        [Xtight,Vtight,Z] = tightenConstraint(sysd,K,X,V,[0.1; 0.1],N);
        Z1(:,N+1) = Z;
%         xLim(:,N+1)
%         Xtight
%         Vtight
    end
end

idx = 0:25;
figure(1)
semilogy(idx,alfa)
hold on;
figure(2);
plot(idx(alfa < 1),xLim(:,alfa < 1))
yticks(-1:0.25:1)
%%






function c = combinations(v)
 elements = v; %cell array with N vectors to combine
 combinations = cell(1, numel(elements)); %set up the varargout result
 [combinations{:}] = ndgrid(elements{:});
 combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
 c = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.

end

function m = maxOverSet(fun, set)
    m = -Inf;
    for i = 1:size(set,1)
        t = fun(set(:,i));
        if t > m
            m = t;
        end
    end
end


function [Xtight,Vtight,Z] = tightenConstraint(sysd,K,X,V,w,N)
A = sysd.A;
B = sysd.B;

Ak = A + B*K;

w = [w, -w];
wLow = min(w,[],2);
wHigh = max(w,[],2);
w = mat2cell(w,ones(size(w,1),1),size(w,2));

wBounds = combinations(w);
wBounds = wBounds';
fun = @(x) norm(x,Inf);

% convert box bounds in form:
% c'*z <= d 
% c'_x*x + c'_u*u <= d
cX = X.cX; % I think we can ignore these if they are just eye
cU = V.cV;
d = [X.d; V.d];

% probably instead use always the same vertex, trying various sequence of
% vertexes lead to a too big number of computations
nState = length(A);
nCst = length(d);
xLim = zeros(nCst,1);
f = zeros(nCst,N*nState);
alfa = max([maxOverSet(fun,Ak^N*wBounds)/maxOverSet(fun,wBounds), maxOverSet(fun,K*Ak^N*wBounds)/maxOverSet(fun,K*wBounds)]);

for i = 0:N-1
    f(:,nState*i+1:nState*(i+1)) = cX*Ak^i + cU*K*Ak^i;
end
opts = optimoptions(@linprog,"Display","off");
Z = zeros(nCst,1);
if alfa < 1
    for s = 1:nCst
        c = -f(s,:);
        llim = repmat(wLow,[N, 1]);
        ulim = repmat(wHigh,[N, 1]);
        [w,th] = linprog(c,[],[],[],[],llim, ulim, opts);
        th = -th;
        Z(s) = 1/(1 - alfa)*th;
        xLim(s) = d(s) - 1/(1 - alfa)*th;
    end
else
    error("Alfa must be within (0,1)")
end

Xtight = xLim(1:length(X.d));
Vtight = xLim(length(X.d)+1:end);
end