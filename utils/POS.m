function [ R,t,s ] = POS( xp,x )
%POS Estimate scaled orthographic projection parameters from 2D-3D
%correspondences
%   xp - 2*n matrix of 2D feature point positions
%   x  - 3*n matrix of 3D feature point positions
%   R, s, t - rotation, scale, translate parameters
%
% The algorithm is loosely based on the first iteration of the POSIT
% algorithm, although I enforce a valid rotation matrix via SVD. This
% function only provides a linear estimate and the reprojection error can
% be improved by subsequently refining the pose estimates using a
% non-linear optimiser. This is achieved by RefineOrthoCam.
%
% Note that SOP has 6 DOF: 3 for rotation, 2 for 2D translation, 1 for
% scale so 3 2D-3D correspondences should be enough. However, our linear
% solution does not enforce any constraints on the 8 estimated parameters
% and so a minimum of 4 points are required.

npts = size(xp,2);

% Build linear system of equations in 8 unknowns of projection matrix
A = zeros(2*npts,8);

A(1:2:2*npts-1,1:3)=x';
A(1:2:2*npts-1,4)=1;

A(2:2:2*npts,5:7)=x';
A(2:2:2*npts,8)=1;

b = reshape(xp,2*npts,1);

% Solve linear system
k = (A\b);

% Extract params from recovered vector
R1 = k(1:3);
R2 = k(5:7);
sTx = k(4);
sTy = k(8);
s = (norm(R1)+norm(R2))/2;
r1 = R1./norm(R1);
r2 = R2./norm(R2);
r3 = cross(r1,r2);
R = [r1'; r2'; r3'];
% Set R to closest orthogonal matrix to estimated rotation matrix
[U,S,V]=svd([r1'; r2'; r3']);
R = U*V';
% Determinant of R must = 1
if (det(R)<0)
    U(3,:)=-U(3,:);
    R = U*V';
end
% Re-estimate scale and translation based on R

% Centroid of rotated 3D points
% Rx = R(1:2,:)*x;
% c3d = mean(Rx,2);
% c2d = mean(xp,2);
% c2d-c3d
% 
% a(1,:) = Rx(1,:)-c3d(1);
% a(2,:) = Rx(2,:)-c3d(2);
% b(1,:) = x(1,:)-c2d(1);
% b(2,:) = x(2,:)-c2d(2);
% s2 = a(:)'\b(:);
% t2 = (c2d-c3d)./s;

% Remove scale from translations
t(1,1)=sTx/s;
t(2,1)=sTy/s;

end

