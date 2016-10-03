function [ bs ] = LandmarkDetector( im )
%LANDMARKDETECTOR Summary of this function goes here
%   Detailed explanation goes here
load('ZhuRamananModel.mat');
%posemap = 90:-15:-90;
%clf; imagesc(im); axis image; axis off; drawnow;
    bs = detect(im, model, model.thresh);
    bs = clipboxes(im, bs);
    bs = nms_face(bs,0.3);
% figure,showboxes(im, bs(1),posemap),title('Highest scoring detection');
% text((bs.xy(:,1)+bs.xy(:,3))/2,(bs.xy(:,2)+bs.xy(:,4))/2,cellstr(num2str([1:68]')))

end


