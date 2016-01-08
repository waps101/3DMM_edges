function b = EstimateShapeSingleSOP( shapePC,shapeMU,shapeEV,R,t,s,numsd,ndims,xp )
%ESTIMATESHAPESINGLESOP Compute optimal shape parameters given R,t,s
%   Solve linear system to find best shape parameters to fit 2D landmarks
%   given estimated rotation, translation and scale

P = reshape(shapePC,3,size(shapePC,1)/3,ndims);
mu = reshape(shapeMU,3,size(shapeMU,1)/3);
A = zeros(2*size(P,2),ndims);
% Loop over all observed points in all images
for j=1:size(P,2)
    % Each observed point adds two equations to the system
    A(2*j-1,:)=s*R(1,1).*squeeze(P(1,j,:))';
    A(2*j-1,:)=A(2*j-1,:)+s*R(1,2).*squeeze(P(2,j,:))';
    A(2*j-1,:)=A(2*j-1,:)+s*R(1,3).*squeeze(P(3,j,:))';
    
    A(2*j,:)=s*R(2,1).*squeeze(P(1,j,:))';
    A(2*j,:)=A(2*j,:)+s*R(2,2).*squeeze(P(2,j,:))';
    A(2*j,:)=A(2*j,:)+s*R(2,3).*squeeze(P(3,j,:))';
    
    h(2*j-1,1)=xp(1,j)- s * (R(1,1).*mu(1,j)+R(1,2).*mu(2,j)+R(1,3).*mu(3,j)+t(1));
    h(2*j,1)=xp(2,j)- s * (R(2,1).*mu(1,j)+R(2,2).*mu(2,j)+R(2,3).*mu(3,j)+t(2));
end
% Unconstrained linear solution:
%b = A\h;

% Hyper-box constrained linear solution:
options = optimset('LargeScale','off','Algorithm', 'active-set', 'Display','off');
C = [eye(ndims); -eye(ndims)];    
d = [numsd.*shapeEV; numsd.*shapeEV];
b = lsqlin(A,h,C,d,[],[],[],[],[],options);

end

