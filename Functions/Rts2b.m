function b = Rts2b( R,t,s )
%RST2B Summary of this function goes here
%   Detailed explanation goes here

r = vrrotmat2vec(R);
b(1:3,1)=(r(1:3)./norm(r(1:3))).*r(4);
b(4:5,1)=t;
b(6,1)=s;

end

