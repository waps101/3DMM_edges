function Ef = meshFaceEdges( faces,edges )
%MESHFACEEDGES Compute faces adjacent to each edge in the mesh
%   faces - nverts by 3 matrix of mesh faces
%   edges - nedges by 2 matrix containing vertices adjacent to each edge
%
% This function is slow! But it only needs to be run once for a morphable
% model and the edge-face list can then be saved

nedges = size(edges,1);

faces = sort(faces,2);
edges = sort(edges,2);

disp('      ');
for i=1:nedges
    idx = find(((faces(:,1)==edges(i,1)) & ( (faces(:,2)==edges(i,2)) | (faces(:,3)==edges(i,2)) )) | ((faces(:,2)==edges(i,1)) & (faces(:,3)==edges(i,2))));
    if length(idx)==1
        idx = [0 idx];
    end
    Ef(i,:)=[idx(1) idx(2)];
    fprintf('\b\b\b\b\b\b%05.2f%%',i/nedges*100);
end

end
