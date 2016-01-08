function dists = RomdhaniObjectivePose( Rts,b,shapeMU,shapePC,occludingVertices,DT,landmarks,xp )
%ROMDHANIOBJECTIVEPOSE Summary of this function goes here
%   Detailed explanation goes here
Rts(4:5)=Rts(4:5).*1e+7;
r = Rts(1:3)';
t = Rts(4:5)';
s = Rts(6);
X = reshape(shapePC*b+shapeMU,3,size(shapePC,1)/3);

if (norm(r)==0) | (real(r)~=r)
    R = eye(3);
else
    R = vrrotvec2mat([r./norm(r) norm(r)]);
end

x2 = R*X(:,occludingVertices);
x2 = x2(1:2,:);
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));
%x2(2,:) = size(edges,1)+1-x2(2,:);
plot(x2(1,:),size(DT,1)+1-x2(2,:),'.b');
drawnow

dists = double(interp2(DT,x2(1,:),x2(2,:),'cubic')); %linear
dists=mean(dists(~isnan(dists)));

%landmarks
weight2=1e-3;
x2 = R*X(:,landmarks);
x2 = x2(1:2,:);
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));
dists(end+1) = weight2.*mean( (xp(1,:)-x2(1,:)).^2 + (xp(2,:)-x2(2,:)).^2 );
dists=sum(dists).*1000000;





end

