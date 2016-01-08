function [landmarks,xp] = read_anchors( filename )
%READ_ANCHORS Summary of this function goes here
%   Detailed explanation goes here
fp = fopen(filename,'r');
if fp == -1
    fclose all;
    error(['Cannot open file ' filename]);
end

tempstr = ' ';
key = 'p';
index= 1;

while ( tempstr ~= -1)
tempstr = fgets(fp);
    if( ~isempty(strfind(tempstr,key)) && isempty(strfind(tempstr,'-')) )
        list = sscanf(tempstr, ['%c %d %d %d %d %d %d'], inf);
        % 4 landmark 5 6 pos
        landmarks(index,:)=[list(4)];
        xp(:,index)=[list(5),list(6)];
        index=index+1;
    end
end

end

