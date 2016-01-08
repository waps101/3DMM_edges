function [ data ] = createRenderings( basedir )
%CREATERENDERINGS Summary of this function goes here
%   Detailed explanation goes here
filePattern = sprintf('%s/*_0_0.anchors', basedir);
baseFileNames = dir(filePattern);
        for i=1:size(baseFileNames,1)
            [landmarks, xp] = read_anchors([basedir baseFileNames(i).name]);
            data(i).id=i;
            data(i).im=imread([basedir baseFileNames(i).name(1:length(baseFileNames(1).name)-7) 'png'], 'BackgroundColor',[1 1 1]);
            data(i).landmarks = landmarks;
            data(i).xp = xp;
        end

end

