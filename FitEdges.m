function [b,R,t,s] = FitEdges( im,xp,landmarks,shapePC,shapeMU,shapeEV,Ef,Ev,tri,ndims,w_prior, w_edges, w_landmarks,niter )
%FITEDGES Perform morphable model fitting using edges
%   This function initialises by fitting to a sparse set of landmarks, then
%   iteratively fits to edges by finding nearest neighbour between
%   projected model edges and image edges
%
% Inputs:
%   im        - Image to fit to
%   xp        - 2 by nlandmarks matrix containing locations of 2D landmarks
%   landmarks - Indices of vertices corresponding to points in xp
%   shapePC   - 3DMM principal components
%   shapeMU   - 3DMM average shape
%   shapeEV   - standard deviations of each 3DMM dimension
%   Ef        - nedges by 2 matrix storing faces adjacent to each edge
%   Ev        - nedges by 2 matrix storing vertices adjacent to each edge
%   tri       - Face structure of 3DMM mesh
%   ndims     - Number of model dimensions to use
%   numsd     - Hard limit on number of standard deviations each parameter
%               may deviate
%
% Outputs:
%   b         - Estimated shape parameters
%   R,t,s     - Estimated rotation, translation and scale of face

FV.faces = tri;
% PARAMETERS
% Number of model dimensions used for initial landmark fit
%ndims = 40;
% Number of edge fitting iterations (may want to replace with a convergence
% test)
%niter = 1;
% Proportion of nearest-neighbour edge matches used at each iteration, may
% be more sensible to threshold on distance to nearest neighbour
percentile = 1;

% Perform initial landmark-only fit
%w_prior = 0.7;
[b,R,t,s] = FitSingleSOP( xp,shapePC,shapeMU,shapeEV,ndims,landmarks,w_prior );
%b=zeros(20,1);
FV.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,size(shapePC,1)/3)';

% Display input image
%figure; imshow(im);

% Extract image edges
edges = edge(rgb2gray(im),'canny',0.15);
%figure; imshow(edges);

% Store pixel locations of edges
[r,c]=find(edges);
r = size(edges,1)+1-r;
%bold=b;
for iter=1:niter
    % Compute vertices lying on occluding boundary
    [ occludingVertices ] = occludingBoundaryVertices( FV,Ef,Ev,R );
        
    % Project occluding boundary vertices
    x2 = R*FV.vertices(occludingVertices,:)';
    x2 = x2(1:2,:);
    x2(1,:) = x2(1,:)+t(1);
    x2(2,:) = x2(2,:)+t(2);
    x2 = x2.*s;
    
    % Find edge correspondences
    [idx,d] = knnsearch([c r],x2'); %its correct, searching based on columns dont change!

    % Filter edge matches - probably want to do something better here
    sortedd=sort(d);
    threshold = sortedd(round(percentile*length(sortedd)));
    idx = idx(d<threshold);
    occludingVertices = occludingVertices(d<threshold);

    [b,R,t,s] = FitSingleSOP( [c(idx)' xp(1,:); r(idx)' xp(2,:)],shapePC,shapeMU,shapeEV,ndims,[occludingVertices; landmarks],w_edges,w_landmarks,length(landmarks));
    % Note: this is completely refitting from scratch. We could just
    % re-start the nonlinear optimisation using previous estimates as
    % initialisation but this doesn't seem to work well.
    %disp(num2str(norm(b-bold)));
    %differb(iter)=norm(b-bold);
    FV.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,size(shapePC,1)/3)';
    %bold=b;
end

end

