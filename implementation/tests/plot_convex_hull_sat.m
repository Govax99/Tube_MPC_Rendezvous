clear;
close all;
%%


[X,Y,Z] = sphere;
points = [];
for i = 1:length(X)
    for j = 1:length(X)
        points = [points, [X(i,j); Y(i,j); Z(i,j)]];
    end
end
C = 5*points;
C = C + [-6; 0; 1];

[X,Y,Z] = sphere;
points = [];
for i = 1:length(X)
    for j = 1:length(X)
        points = [points, [X(i,j); Y(i,j); Z(i,j)]];
    end
end
C2 = 6*points;
C2 = C2 + [0; 0; 1];

[X,Y,Z] = sphere;
points = [];
for i = 1:length(X)
    for j = 1:length(X)
        points = [points, [X(i,j); Y(i,j); Z(i,j)]];
    end
end
C3 = 6*points;
C3 = C3 + [8; 0; 1.5];

[X,Y,Z] = ellipsoid(2.5,0,0.5,15,6.75,4.75);
points = [];
for i = 1:length(X)
    for j = 1:length(X)
        points = [points, [X(i,j); Y(i,j); Z(i,j)]];
    end
end
C4 = points;

% Read mesh and compute normal vector
TR = stlread('ENVISAT.stl'); % contains "Mesh" structure
[X,~,J] = unique(TR.Points,'rows');
F = J(TR.ConnectivityList).';
T = X.';
T(3,:) = T(3,:)+1;
T = 5*T;

R1 = [  1.0000000,  0.0000000,  0.0000000;
   0.0000000,  0.0000000, -1.0000000;
   0.0000000,  1.0000000,  0.0000000 ];
R2 = [  0.0000000,  1.0000000,  0.0000000;
  -1.0000000,  0.0000000,  0.0000000;
   0.0000000,  0.0000000,  1.0000000 ];
T = R2*R1*T;
DT = delaunayTriangulation(T');
%%

sat = Polyhedron(DT.Points);
[X,Y,Z] = sphere(4);
r = 1.5*sqrt(2);
X = X * r;
Y= Y * r;
Z = Z * r;
sph = [];
for i = 1:length(X)
    for j = 1:length(Y)
        sph = [sph, [X(i,j); Y(i,j); Z(i,j)]];
    end
end

sph = unique(sph.','rows','stable').';
sph = Polyhedron(sph');

new = sph + sat;

%%

% VISUALISE RESULTS ONLY IN MATLAB
if(exist('OCTAVE_VERSION', 'builtin') == 0)
       % .. create new figure
       figure('units','centimeters', 'WindowStyle','normal', 'color','w',...
       'Position',[0 8.5 9 6],'defaultAxesColorOrder',parula,...
       'Renderer','opengl') 
       % .. adjust properties
       axis equal tight off; hold all; 
       % .. display ellipsoid and convex hull
       DT = delaunayTriangulation(T');
       [K,~] = convexHull(DT);
       trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),...
              'EdgeColor','none','FaceColor',[.4 1 .9 ],...
              'FaceLighting','flat' )

       DT = delaunayTriangulation(C4');
       [K,~] = convexHull(DT);
       trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),...
              'EdgeColor','none','FaceColor',[.4 1 .8 ],...
              'FaceLighting','flat' , 'FaceAlpha', 0.5)

       % ... adjust point of view   
       view(42,21)
       % ... add light
       light('Position',[5 -10 20],'Style','local'); 
end