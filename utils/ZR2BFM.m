function [ xp,landmarkidx ] = ZR2BFM( bs,im )
%ZR2BFM Convert Zhu and Ramamnan output to BFM vertex correspondence
%   Converts the detected landmarks from the Zhu and Ramamnan detector into
%   references into the BFM and 2D positions in the format used by our
%   code. There is only support for 68 landmark detections at present.

% No support for profile poses yet
if size(bs.xy,1)==68
    % Frontal pose
    
    listX= (bs.xy(:,1)+bs.xy(:,3))/2;
    listY= (bs.xy(:,2)+bs.xy(:,4))/2;
    list = [listX'; listY'];
    
    landlist= [1:51, 56, 64];
    xp=[list(1,landlist);list(2,landlist)];
    xp(2,:) = size(im,1)+1-xp(2,:);
    
    landmarkidx = [8333,7301,6011,9365,10655,8320,8311,8286,8275,5959,4675,3642,4922,3631,2088,27391,39839,40091,40351,6713,10603,11641,12673,11244,12661,14472,27778,41804,41578,41310,9806,8345,7442,6936,5392,7335,7851,8354,9248,9398,11326,9399,9129,9406,9128,8890,8367,7858,7580,7471,8374,23002,32033]';    
end

end

