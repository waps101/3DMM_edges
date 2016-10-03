function texture = faceTexture( FV,R,t,s,im)
%FACETEXTURE Summary of this function goes here
%   Detailed explanation goes here

rotpts = R*FV.vertices';
pts2D =[s.*(rotpts(1,:)+t(1)); s.*(rotpts(2,:)+t(2))];
[rows,cols]=meshgrid(1:size(im,2),1:size(im,1));

FV.facevertexcdata = zeros(0);

for col=1:3%size(im,1)+1-
FV.facevertexcdata(:,col) = interp2(rows,cols,im(:,:,col),pts2D(1,:),size(im,1)+1-pts2D(2,:));
end

%non visible texture set to zero
nonVisibleVertices = setdiff(1:length(FV.vertices), visiblevertices(FV, R));
for col=1:3
FV.facevertexcdata(nonVisibleVertices,col) = NaN; 
end

texture = FV.facevertexcdata;
end

