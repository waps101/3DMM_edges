function [ b ] = RefineOrthoCam( x,X,b0 )
%ORTHOCAMNONLIN Summary of this function goes here
%   Detailed explanation goes here

% Initialise b
if nargin<3
    b0 = [0 0 0 0 0 1];
end

options = optimset('TolFun',1e-20,'Display','off'); %,'LargeScale','off');

%sum(OrthoCamNonlinerrfun( b0,x,X ).^2)
b = lsqnonlin(@(b) OrthoCamNonlinerrfun( b,double(x),double(X) ),double(b0),[-Inf -Inf -Inf -Inf -Inf 0],[Inf Inf Inf Inf Inf Inf],options);
%sum(OrthoCamNonlinerrfun( b,x,X ).^2)

end

function result = OrthoCamNonlinerrfun( b,x,X )
%ORTHOCAMNONLINERRFUN Summary of this function goes here
%   Detailed explanation goes here

r = b(1:3)';
t = b(4:5);
s = b(6);


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