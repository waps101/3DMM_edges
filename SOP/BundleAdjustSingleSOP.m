function [ b,R,t,s ] = BundleAdjustSingleSOP( x,b0,R0,t0,s0,shapePC,shapeMU,shapeEV,numsd )
%ORTHOCAMNONLIN Summary of this function goes here
%   Detailed explanation goes here

b0 = [Rts2b(double(R0),double(t0),double(s0)); b0];


options = optimset('TolFun',1e-20,'Display','off'); %,'LargeScale','off');
%options = optimset('LargeScale','off','Display','off'); 
%options = optimset('TolX',1e-4,'MaxIter',1e4,'MaxFunEvals',1e4,'TolFun',1e-4,'Algorithm','trust-region-reflective'); 

LB=[-Inf -Inf -Inf -Inf -Inf 0 -numsd.*shapeEV']';
UB=[Inf Inf Inf Inf Inf Inf numsd.*shapeEV']';

%sum(BANonlinerrfun( b0,double(x),shapePC,shapeMU).^2)
b = lsqnonlin(@(b) BANonlinerrfun( b,double(x),shapePC,shapeMU),double(b0),LB,UB,options);
%sum(BANonlinerrfun( b,double(x),shapePC,shapeMU).^2)

[R,t,s]=b2Rts(b(1:6));
b = b(7:end);

end

function result = BANonlinerrfun( b,x,shapePC,shapeMU )
%ORTHOCAMNONLINERRFUN Summary of this function goes here
%   Detailed explanation goes here

r = b(1:3)';
t = b(4:5);
s = b(6);

X = reshape(shapePC*b(7:end)+shapeMU,3,size(shapePC,1)/3);

if (norm(r)==0) | (real(r)~=r)
    R = eye(3);
else
    R = vrrotvec2mat([r./norm(r) norm(r)]);
end

x2=R*X;
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));

result = [(x(1,:)-x2(1,:)) (x(2,:)-x2(2,:))];

end
