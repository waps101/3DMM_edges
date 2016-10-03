function [ b,R,t,s ] = BundleAdjustSingleSOP( x,b0,R0,t0,s0,shapePC,shapeMU,shapeEV,numsd,w1,w2,nlandmarks )
%ORTHOCAMNONLIN Summary of this function goes here
%   w2 and nlandmarks are optional

b0 = [Rts2vec(double(R0),double(t0),double(s0)); b0];


options = optimset('Display','iter'); 

LB=[-inf -inf -inf -inf -inf 0 -numsd.*shapeEV'];
UB=[inf inf inf inf inf inf numsd.*shapeEV'];

if nargin==10
    b = lsqnonlin(@(b) BANonlinerrfun( b,double(x),shapePC,shapeMU,shapeEV,w1 ),double(b0),LB,UB,options);
elseif nargin==12
    b = lsqnonlin(@(b) BANonlinerrfun( b,double(x),shapePC,shapeMU,shapeEV,w1,w2,nlandmarks ),double(b0),LB,UB,options);
end

[R,t,s]=vec2Rts(b(1:6));
b = b(7:end);

end

function result = BANonlinerrfun( b,x,shapePC,shapeMU,shapeEV,w1,w2,nlandmarks )
%ORTHOCAMNONLINERRFUN Summary of this function goes here
%   Detailed explanation goes here

[R,t,s]=vec2Rts(b(1:6));

X = reshape(shapePC*b(7:end)+shapeMU,3,size(shapePC,1)/3);

x2=R*X;
x2(1,:)=s.*(x2(1,:)+t(1));
x2(2,:)=s.*(x2(2,:)+t(2));

if nargin==6
    w=w1;
    w1=(1-w)/(2*size(x2,2));
    w2 = w/(length(shapeEV));

    result = w1.*[(x(1,:)-x2(1,:)) (x(2,:)-x2(2,:))];

    result = [result, (w2.*(b(7:end)./(shapeEV)))'];
elseif nargin==8
    w3 = 1-w1-w2;
    nedges = size(x,2)-nlandmarks;
    w1 = w1./(2*nedges);
    w2 = w2./(2*nlandmarks);
    w3 = w3./(length(shapeEV));
    residuals_edges = w1.*[(x(1,1:nedges)-x2(1,1:nedges)) (x(2,1:nedges)-x2(2,1:nedges))];
    residuals_landmarks = w2.*[(x(1,nedges+1:end)-x2(1,nedges+1:end)) (x(2,nedges+1:end)-x2(2,nedges+1:end))];
    residuals_prior = (w3.*(b(7:end)./(shapeEV)))';
    result = [residuals_edges residuals_landmarks residuals_prior];
end
%w2=1000;
%result(end+1)= (w2.*norm(b(7:end)./(shapeEV).^2).^2);


end
