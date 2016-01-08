function fn = facenormals( V,F )
%FACENORMALS Compute face normals of triangular mesh
%   V - nverts by 3 matrix containing vertex positions
%   F - ntri by 3 matrix containing triangle vertex indices

% Get the triangle vertices
v1      = F(:, 1);
v2      = F(:, 2);
v3      = F(:, 3);

% Compute the edge vectors
e1s = V(v2, :) - V(v1, :);
e2s = V(v3, :) - V(v1, :);

% Compute cross products between edge vectors
fn    = cross(e1s, e2s, 2);

% Normalise face normals to unit length
norms = sqrt(fn(:,1).^2+fn(:,2).^2+fn(:,3).^2);
fn=fn./repmat(norms,[1 3]);

end

