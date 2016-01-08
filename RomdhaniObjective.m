function dists = RomdhaniObjective(b,shapeMU,shapePC,shapeEV,occludingVertices,edges,DT)
%ROMDHANIOBJECTIVE Summary of this function goes here
%   Detailed explanation goes here
%b=b.*shapeEV;
r = b(1:3)';
t = b(4:5)';
s = b(6);

X = reshape(shapePC*(b(7:end).*shapeEV)+shapeMU,3,size(shapePC,1)/3);

if (norm(r)==0) | (real(r)~=r)
    R = eye(3);
else
    R = vrrotvec2mat([r./norm(r) norm(r)]);
end

X = X(:,occludingVertices);
x2=R*X;
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));
%x2(2,:) = size(edges,1)+1-x2(2,:);
plot(x2(1,:),size(edges,1)+1-x2(2,:),'.r');
drawnow

dists = double(interp2(DT,x2(1,:),x2(2,:),'cubic')); %linear
% dists = sum(dists.^2);
end

