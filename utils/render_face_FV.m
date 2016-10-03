function im = render_face_FV(FV, T, oglp, background)

% Get the rendering parameters
width           = oglp.width;               % width of the image plane
height          = oglp.height;              % height of the image plane
i_amb_light     = oglp.i_amb_light;   	% ambient light intensity
i_dir_light     = oglp.i_dir_light;    % directed light intensity
d_dir_light     = [0; 0; 1];       % directed light direction
shininess       = 1;           % shininess
specularity     = 0;      % specularity
do_cast_shadows = false;     % flag for casting shadows
%shadow_buf_size = oglp.sbufsize;            % buffer size for casting shadows

% Get the extrinsic transformation matrix
M = T(1: 3, :); % the output does not need to be in homogeneous coordinates

% Get the vertices
V           = FV.vertices;
Nvertices   = size(FV.vertices, 1);

% Compute the transformed vertices
V(:, 4)	= 1;        % use homogeneous coordinates for input
V2   	= V * M.';	% the vertices are transposed

clear V

% Store the vertex depths for z-buffering
Z = V2(:, 3);

% Compute the projected vertices in the image plane
UV(:, 1)    = V2(:, 1) ;	% orthographic projection for x
UV(:, 2)    = V2(:, 2) ;	% orthographic projection for y

% Transform to the pixel plane (the axes remain switched)
%UV(:, 1)    = width - UV(:, 1);
%UV          = UV + 0.5;

% Get the triangle vertices
v1      = FV.faces(:, 1);
v2      = FV.faces(:, 2);
v3      = FV.faces(:, 3);
Nfaces  = size(FV.faces, 1);

% Compute bounding boxes for the projected triangles
x       = [UV(v1, 1), UV(v2, 1), UV(v3, 1)];
y       = [UV(v1, 2), UV(v2, 2), UV(v3, 2)];
minx    = ceil (min(x, [], 2));
maxx    = floor(max(x, [], 2));
miny    = ceil (min(y, [], 2));
maxy    = floor(max(y, [], 2));

clear x y

% Frustum culling
minx    = max(1,        minx);
maxx    = min(width,    maxx);
miny    = max(1,        miny);
maxy    = min(height,   maxy);

% Construct the pixel grid (can speed up by precomputing if shared among the images)
[rows, cols] = meshgrid(1: width, 1: height);

% Initialize the depth-, face- and weight-buffers
zbuffer     = -inf(height, width);
fbuffer     = zeros(height, width);
wbuffer1    = NaN(height, width);
wbuffer2    = NaN(height, width);
wbuffer3    = NaN(height, width);

% For each triangle (can speed up by comparing the triangle depths to the z-buffer and priorly sorting the triangles by increasing depth)
for i = 1: Nfaces
    
    % If some pixels lie in the bounding box
    if minx(i) <= maxx(i) && miny(i) <= maxy(i)
        
        % Get the pixels lying in the bounding box
        px = rows(miny(i): maxy(i), minx(i): maxx(i));
        py = cols(miny(i): maxy(i), minx(i): maxx(i));
        px = px(:);
        py = py(:);
        
        % Compute the edge vectors
        e0 = UV(v1(i), :);
        e1 = UV(v2(i), :) - e0;
        e2 = UV(v3(i), :) - e0;
        
        % Compute the barycentric coordinates (can speed up by first computing and testing a solely)
        det     = e1(1) * e2(2) - e1(2) * e2(1);
        tmpx    = px - e0(1);
        tmpy    = py - e0(2);
        a       = (tmpx * e2(2) - tmpy * e2(1)) / det;
        b       = (tmpy * e1(1) - tmpx * e1(2)) / det;
        
        % Test whether the pixels lie in the triangle
        test = a >= 0 & b >= 0 & a + b <= 1;
        
        % If some pixels lie in the triangle
        if any(test)
            
            % Get the pixels lying in the triangle
            px = px(test);
            py = py(test);
            
            % Interpolate the triangle depth for each pixel
            w2 = a(test);
            w3 = b(test);
            w1 = 1 - w2 - w3;
            pz = Z(v1(i)) * w1 + Z(v2(i)) * w2 + Z(v3(i)) * w3;
            
            % For each pixel lying in the triangle
            for j = 1: length(pz)
                
                % Frustum culling
%                if pz(j) <= -near && pz(j) >= -far
                    
                    % Update the depth-, face- and weight-buffers
                    if pz(j) > zbuffer(py(j), px(j))
                        zbuffer(py(j), px(j))   = pz(j);
                        fbuffer(py(j), px(j))   = i;
                        wbuffer1(py(j), px(j))  = w1(j);
                        wbuffer2(py(j), px(j))  = w2(j);
                        wbuffer3(py(j), px(j))  = w3(j);
                    end
                    
%                end
                
            end
            
        end
        
    end
    
end

clear UV Z zbuffer px py pz minx maxx miny maxy

% Get the vertices to render
test    = fbuffer ~= 0;
f       = unique(fbuffer(test));
v       = unique([v1(f); v2(f); v3(f)]);
f       = find(any(ismember(FV.faces, v), 2));
Nfaces  = length(f);
 
% Compute the edge vectors
e1s = V2(v2(f), :) - V2(v1(f), :);
e2s = V2(v3(f), :) - V2(v1(f), :);
e3s = V2(v2(f), :) - V2(v3(f), :);

clear V2

% Normalize the edge vectors
e1s_norm = e1s ./ repmat(sqrt(sum(e1s.^2, 2)), 1, 3);
e2s_norm = e2s ./ repmat(sqrt(sum(e2s.^2, 2)), 1, 3);
e3s_norm = e3s ./ repmat(sqrt(sum(e3s.^2, 2)), 1, 3);

% Compute the angles
angles(:, 1) = acos(sum(e1s_norm .* e2s_norm, 2));
angles(:, 2) = acos(sum(e3s_norm .* e1s_norm, 2));
angles(:, 3) = pi - (angles(:, 1) + angles(:, 2));

% Compute the triangle weighted normals
triangle_normals    = cross(e1s, e3s, 2);
w1_triangle_normals = triangle_normals .* repmat(angles(:, 1), 1, 3);
w2_triangle_normals = triangle_normals .* repmat(angles(:, 2), 1, 3);
w3_triangle_normals = triangle_normals .* repmat(angles(:, 3), 1, 3);

clear e1s e2s e3s e1s_norm e2s_norm e3s_norm angles triangle_normals

% Initialize the vertex normals
normals = zeros(Nvertices, 3);

% Update the vertex normals
for i = 1: Nfaces
    normals(v1(f(i)), :) = normals(v1(f(i)), :) + w1_triangle_normals(i, :);
    normals(v2(f(i)), :) = normals(v2(f(i)), :) + w2_triangle_normals(i, :);
    normals(v3(f(i)), :) = normals(v3(f(i)), :) + w3_triangle_normals(i, :);
end

clear w1_triangle_normals w2_triangle_normals w3_triangle_normals

% Normalize the vertex normals
%normals = normals(v, :);
normals = normals ./ repmat(sqrt(sum(normals.^2, 2)), 1, 3);

% Get the texture
texture1 = zeros(Nvertices,1);
texture2 = zeros(Nvertices,1);
texture3 = zeros(Nvertices,1);

H = [0; 0; 1] + d_dir_light;
H = H / norm(H, 2);

% Compute the dot products between the vertex normals and the directed light source
NdotL = d_dir_light(1) * normals(:, 1) + d_dir_light(2) * normals(:, 2) + d_dir_light(3) * normals(:, 3);

% Compute the dot products between the vertex normals and the directed light source
diffuse = max(0,normals(:, 3));

diffuse1 = i_dir_light(1) *  NdotL;
diffuse2 = i_dir_light(2) *  NdotL;
diffuse3 = i_dir_light(3) *  NdotL;
% Render the texture
texture1 =  diffuse1;
texture2 =  diffuse2;
texture3 =  diffuse3;


% Render the texture
% texture1 = diffuse;
% texture2 = diffuse;
% texture3 = diffuse;

clear ambient1 diffuse1 specular1 ambient2 diffuse2 specular2 ambient3 diffuse3 specular3 v

% Initialize the image
im1 = NaN(height, width);
im2 = NaN(height, width);
im3 = NaN(height, width);


% Rasterize the image
v1          = v1(fbuffer(test));
v2          = v2(fbuffer(test));
v3          = v3(fbuffer(test));
w1          = wbuffer1(test);
w2          = wbuffer2(test);
w3          = wbuffer3(test);
im1(test)   = w1 .* texture1(v1) + w2 .* texture1(v2) + w3 .* texture1(v3);
im2(test)   = w1 .* texture2(v1) + w2 .* texture2(v2) + w3 .* texture2(v3);
im3(test)   = w1 .* texture3(v1) + w2 .* texture3(v2) + w3 .* texture3(v3);

clear v1 v2 v3 w1 w2 w3 

% Reshape the image
im = cat(3, im1, im2, im3);

clear im1 im2 im3

%im = flipud(fliplr(im));
im = flipud(im);
test = repmat(flipud(test),[1 1 3]);

im(~test)=background(~test);

end
