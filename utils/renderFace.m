function [ renderim ] = renderFace( FV,im,R,t,s,istexture )
%RENDERFACE Orthographic projection of the mesh on the image
%   Inputs  
%   FV - face mesh
%   im - face image
%   R,t,s - pose parameteres
%   istexture - control for FV has facevertexcdata field.
%       true; renders the face with texture.
%       false; renders only the face shape.
%
% Updated Sep 2016

if(~isa(im,'double'))
im=double(im)./255;
end

oglp.height=size(im,1);
oglp.width=size(im,2);
oglp.i_amb_light = [1 1 1];
oglp.i_dir_light = [1 1 1];

Rr = R;
Rr(4,4)=1;
Sr = eye(4).*s;
Tr = eye(4);
Tr(1:2,4)=t;
T = Tr*Sr*Rr;
clear Tr Sr Rr

if(istexture)
FV.facevertexcdata = faceTexture(FV,R,t,s,im);
oglp.i_dir_light = [0 0 0];
renderim = render_face_back(FV, T, oglp, zeros(size(im)));
else
renderim = render_face_FV(FV, T, oglp, im);
end

end

