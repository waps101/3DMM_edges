clear FV

%% ADJUSTABLE PARAMETERS

% Number of model dimensions to use
ndims = 50;

% Number of standard deviations that each parameter is allowed to deviate
% from the mean
numsd = 0.5;

%% 

addpath('BoundaryVertices','RamananDetector','Functions','SOP');

% YOU MUST set this to the base directory of the Basel Face Model
BFMbasedir = '/Users/williamsmith/Documents/arnaud/dbs/BaselFaceModel/';

% Load morphable model
load(strcat(BFMbasedir,'01_MorphableModel.mat'));

% Load edge structure to speed up finding occluding boundaries

load('BFMedgestruct.mat')

%% Setup basic parameters:

landmarks=[8333,7301,6011,9365,10655,8320,8311,8286,8275,5959,4675,3642,4922,3631,2088,27391,39839,40091,40351,6713,10603,11641,12673,11244,12661,14472,27778,41804,41578,41310,9806,8345,7442,6936,5392,7335,7851,8354,9248,9398,11326,9399,9129,9406,9128,8890,8367,7858,7580,7471,8374,23002,32033]';
FV.faces = tl;

%% Landmark Detection

im = imread('pie.png');
xp= RamananDetector(im);

%% Fit model

[b,R_est,t_est,s_est] = FitEdges( im,xp,landmarks,shapePC,shapeMU,shapeEV,Ef,Ev,tl,ndims,numsd );

%% Display fitted model

FV.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,53490)';
figure; patch(FV, 'FaceColor', [1 1 1],'EdgeColor', 'none', 'FaceLighting', 'phong'); axis equal; axis off;
light('Position',[0 0 1],'Style','infinite');
patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong','AmbientStrength',0,'DiffuseStrength',1,'SpecularStrength',0,'BackFaceLighting','lit');

 %% Texture sampling
im = double(im)./255;
pts=FV.vertices';
rotpts=R_est*pts;
pts2D = [s_est.*(rotpts(1,:)+t_est(1)); s_est.*(rotpts(2,:)+t_est(2))];

[rows,cols]=meshgrid(1:size(im,2),1:size(im,1));

for col=1:3
    FV.facevertexcdata(:,col) = interp2(rows,cols,im(:,:,col),pts2D(1,:),size(im,1)+1-pts2D(2,:));
end
nonVisibleVertices= setdiff(1:53490, visiblevertices(FV, R_est));

for col=1:3
    FV2.facevertexcdata(nonVisibleVertices,col) = NaN;
end

figure; patch(FV, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceLighting', 'phong'); axis equal; axis off;
