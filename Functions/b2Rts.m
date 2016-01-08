
function [ R,t,s ] = b2Rts( b )
%B2RST Summary of this function goes here
%   Detailed explanation goes here

r = b(1:3)';
t = b(4:5)';
s = b(6);
if (norm(r)==0) | (real(r)~=r)
    R = eye(3);
else
    R = vrrotvec2mat([r./norm(r) norm(r)]);
end

end

