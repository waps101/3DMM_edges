function dists = RomdhaniObjectiveShape(b,xp,shapeEV,shapeMU,shapePC,R,t,s,occludingVertices,DT,landmarks )
%ROMDHANIOBJECTIVESHAPE Summary of this function goes here
%   Detailed explanation goes here

X = reshape(shapePC*b+shapeMU,3,size(shapePC,1)/3);

x2 = R*X(:,occludingVertices);
x2 = x2(1:2,:);
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));

plot(x2(1,:),size(DT,1)+1-x2(2,:),'.r');
drawnow


dists = double(interp2(DT,x2(1,:),x2(2,:),'cubic')); %linear
%normalization
%dists=mean(dists(~isnan(dists)));
dists=mean(dists(~isnan(dists))).*10;
%dists = mean(dists.^2);
%34
%22
weight=1e-3;
weight2=1e-3;


dists(end+1)= weight.*sum(b.^2./shapeEV.^2);

%landmarks
%weight2=1e-3;

x2 = R*X(:,landmarks);
x2 = x2(1:2,:);
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));
dists(end+1) = weight2.*mean( (xp(1,:)-x2(1,:)).^2 + (xp(2,:)-x2(2,:)).^2 );
%dists=sum(dists).*10000;
dists=sum(dists).*1000000;
%[idx,d] = knnsearch(xp',x2');
%dists(end+1)= sum(d.^2);
%D = pdist2(xp',double(x2'));
%dists(end+1) = sqrt(sum((xp - x2) .^ 2));
end
