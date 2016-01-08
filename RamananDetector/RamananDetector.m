function [ xp ] = RamananDetector( im )
%RAMANANDETECTOR Summary of this function goes here
%   Detailed explanation goes here
load('RamananModel.mat')
model.interval = 5;
model.thresh = min(-0.65, model.thresh);
posemap = 90:-15:-90;
%clf; imagesc(im); axis image; axis off; drawnow;
    bs = detect(im, model, model.thresh);
    bs = clipboxes(im, bs);
    bs = nms_face(bs,0.3);
%figure,showboxes(im, bs(1),posemap),title('Highest scoring detection');
%text((bs.xy(:,1)+bs.xy(:,3))/2,(bs.xy(:,2)+bs.xy(:,4))/2,cellstr(num2str([1:68]')))
listX= (bs.xy(:,1)+bs.xy(:,3))/2;
listY= (bs.xy(:,2)+bs.xy(:,4))/2;
list = [listX'; listY'];
im = double(im)./255;
for i=1:51
seqX(i)=listX(i)';
seqY(i)=listY(i)';
end
seqX(52)=listX(56)';seqX(53)=listX(64)';
seqY(52)=listY(56)';seqY(53)=listY(64)';
seq = [seqX; seqY];
xp=seq;
xp(2,:) = size(im,1)+1-xp(2,:);

end

