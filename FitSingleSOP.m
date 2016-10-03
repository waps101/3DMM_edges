function [b,R,t,s] = FitSingleSOP( xp,shapePC,shapeMU,shapeEV,ndims,landmarks,w1,w2,nlandmarks )
%FITSINGLESOP Fit morphable model to single image landmarks
%   Given 2D locations of landmark vertices, fit morphable model under
%   scaled orthographic projection
%
% Inputs:
% xp - 2 by nlandmarks matrix of 2D landmark positions
% shapePC, shapeMU, shapeEV - 3D morphable model
% ndims - number of model dimensions to fit
% landmarks - vector of landmark vertex indices
%
% Outputs:
% b - shape parameter vector
% R - 3 by 3 rotation matrix
% t - 2 by 1 2D translation vector
% s - scale

% Make sure xp is 2 by nlandmarks
if size(xp,1)>size(xp,2)
    xp = xp';
end

% Parameters
numsd= 3; % Number of standard deviations for hyperbox constraint
niter = 10; % Number of iterations of coordinate ascent

nverts = size(shapePC,1)/3;

% Subselect landmark vertices from mean and principal component vectors
sortedfps = reshape(1:nverts*3,3,nverts)';
fps_sel = sortedfps(landmarks,1:3)';
fps_sel = fps_sel(:);
sortedfps = fps_sel;
shapePC = double(shapePC(sortedfps,1:ndims));
shapeMU = double(shapeMU(sortedfps,:));
shapeEV = double(shapeEV(1:ndims));

% Obtain initial estimate of pose parameters using mean shape
disp('Initialising camera parameters...')
x = reshape(shapeMU,3,size(shapePC,1)/3);
[ R,t,s ] = EstimateSOPwithRefinement( xp,x );

% Obtain initial shape estimate by solving linear system
disp('Initialising shape parameters...')
b = EstimateShapeSingleSOP( shapePC,shapeMU,shapeEV,R,t,s,numsd,ndims,xp );

% Perform co-ordinate ascent, alternatively solving for pose and shape
% whilst fixing the other
fprintf('Coordinate ascent... Iteration: ');
for k=1:niter
    
    x = reshape(shapePC*b+shapeMU,3,size(shapePC,1)/3);
    [ R,t,s ] = EstimateSOPwithRefinement( xp,x );
    
    b = EstimateShapeSingleSOP( shapePC,shapeMU,shapeEV,R,t,s,numsd,ndims,xp );
%    disp(['Residual fitting error: ' num2str(norm(A*b-h)^2)])
    fprintf('%d ',k);
end
fprintf('\n');

% Perform non-linear bundle adjustment, simultaneously optimising shape and
% pose parameters
disp('Performing bundle adjustment...')
if nargin==7
[ b,R,t,s ] = BundleAdjustSingleSOP( xp,b,R,t,s,shapePC,shapeMU,shapeEV,numsd,w1 );
elseif nargin==9
[ b,R,t,s ] = BundleAdjustSingleSOP( xp,b,R,t,s,shapePC,shapeMU,shapeEV,numsd,w1,w2,nlandmarks );    
end
disp('Done.')
end

