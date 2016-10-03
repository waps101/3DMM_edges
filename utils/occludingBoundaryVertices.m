function [ occludingVertices ] = occludingBoundaryVertices( FV,Ef,Ev,R )
%OCCLUDINGBOUNDARYVERTICES Return vertices in FV that lie on boundary edges
%   FV - mesh structure containing vertices and faces
%   Ef - faces adjacent to each edge in the mesh
%   Ev - vertices adjacent to each face in the mesh
%
% This runs during edge fitting. As shape/pose changes, so do the vertices
% that lie on the occluding boundary. So these probably want to be updated
% every iteration.

if nargin>3
    FV.vertices = (R*FV.vertices')';
end

% Compute face normals
fn = facenormals( FV.vertices,FV.faces );

% Edges with a zero index lie on the mesh boundary, i.e. they are only
% adjacent to one face
boundaryEdges = Ef(:,1)==0;

% To avoid matlab indexing problems, replace zero index with an arbitrary
% value (these edges are going to be ignored anyway)
idx1 = Ef(:,1);
idx1(boundaryEdges)=1;

% Compute the occluding edges as those where the two adjacent face normals
% differ in the sign of their Z component
occludingEdges = (sign(fn(idx1,3)) ~= sign(fn(Ef(:,2),3))) & ~boundaryEdges;

% Select the vertices lying at the end of the occluding edges and remove
% duplicates
occludingVertices = Ev(occludingEdges,:);
occludingVertices = unique(occludingVertices(:));

% Remove vertices from occluding boundary list that are not visible
occludingVertices = intersect(occludingVertices,visiblevertices(FV, eye(3)));

end

