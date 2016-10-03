function [ R,t,s,b ] = optimiseHardEdgeCost( b0,x_landmarks,shapeEV,shapeMU,shapePC,R0,t0,s0,r,c,landmarks,options,faces,runtoconvergence )
%OPTIMISEHARDEDGECOST Nonlinear LS optimisation of hard edge cost
%   Inputs
%      b0 - initial estimate of shape parameters
%      x_landmarks - 2 by k matrix of 2D positions of landmarks
%      shapeEV - standard deviations of shape model parameters
%      shapeMU - mean shape
%      shapePC - shape principal components
%      R0,t0,s0 - initial pose estimate
%      edgeim - edge image
%      landmarks - model indices of vertices corresponding to landmarks
%      options
%          occludingVertices - (optional) fixed occluding boundary verts
%          Ef,Ev - (optional) edge structures of mesh
%          w1 - edge weight
%          w2 - landmark weight
%          Note: w1 + w2 <= 1
%      faces - mesh face structure
%
% Updated Sep 2016

numsd = 3;
ndims = length(b0);
nedges = length(options.occludingVertices);
nlandmarks = length(landmarks);
J = zeros(2*nedges+2*nlandmarks+ndims,6+ndims);
J(1:2*nedges+2*nlandmarks,:)=1;
J(2*nedges+2*nlandmarks+1:end,7:end)=eye(ndims);

if runtoconvergence
    optoptions = optimoptions('lsqnonlin','Display','iter','JacobPattern',J);
else
    optoptions = optimoptions('lsqnonlin','Display','iter','MaxIter',5,'JacobPattern',J);
end


% Rescale shape params to be in units of standard deviations
b0 = b0./shapeEV(1:ndims);
b0 = [Rts2vec(R0,t0,s0); b0];

% With no re-scaling of shape params
%LB=[-inf -inf -inf -inf -inf 0 -numsd.*shapeEV(1:ndims)'];
%UB=[inf inf inf inf inf inf numsd.*shapeEV(1:ndims)'];
% With re-scaling of shape params 
LB=[-inf -inf -inf -inf -inf 0 -numsd.*ones(1,ndims)];
UB=[inf inf inf inf inf inf numsd.*ones(1,ndims)];


b = lsqnonlin(@(b) hardEdgeCost(b,x_landmarks,shapeEV(1:ndims),shapeMU,shapePC(:,1:ndims),r,c,landmarks,options,faces), b0, LB, UB, optoptions);
[R,t,s]=vec2Rts(b(1:6));

b = b(7:end);
b = b.*shapeEV(1:ndims);

end

function residuals = hardEdgeCost(b,x_landmarks,shapeEV,shapeMU,shapePC,r,c,landmarks,options,faces)

[R,t,s]=vec2Rts(b(1:6));
b = b(7:end);
b = b.*shapeEV(1:length(b));

X = reshape(shapePC*b+shapeMU,3,size(shapePC,1)/3);

if ~isfield(options,'occludingVertices')
    % We have not been supplied with fixed occluding boundary vertices
    % So recompute at every iteration (slow and objective discontinuous)
    FV.faces = faces;
    FV.vertices = X';
    [ options.occludingVertices ] = occludingBoundaryVertices( FV,options.Ef,options.Ev,R );
end

% Scale weights to be invariant to number of landmarks/edge vertices/number
% of dimensions
nedges = length(options.occludingVertices);
nlandmarks = length(landmarks);
w3 = 1-options.w1-options.w2;
w1 = options.w1./(2*nedges);
w2 = options.w2./(2*nlandmarks);
w3 = w3./(length(shapeEV));

% Compute position of projected occluding boundary vertices
x_edge = R*X(:,options.occludingVertices);
x_edge = x_edge(1:2,:);
x_edge(1,:)=s.*(x_edge(1,:)+t(1));
x_edge(2,:)=s.*(x_edge(2,:)+t(2));
% Find edge correspondences
[idx,d] = knnsearch([c r],x_edge');
edgeResiduals = double([c(idx)'-x_edge(1,:) r(idx)'-x_edge(2,:)]);

% Compute position of projected landmark vertices
x_landmarks2 = R*X(:,landmarks);
x_landmarks2 = x_landmarks2(1:2,:);
x_landmarks2(1,:)=s.*(x_landmarks2(1,:)+t(1));
x_landmarks2(2,:)=s.*(x_landmarks2(2,:)+t(2));
landmarkResiduals = double([x_landmarks2(1,:)-x_landmarks(1,:) x_landmarks2(2,:)-x_landmarks(2,:)]);

priorResiduals = double(b./shapeEV(1:length(b)));
residuals = [w1.*edgeResiduals w2.*landmarkResiduals w3.*priorResiduals'];

end