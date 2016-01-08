function DT = GetDT( im,k )

%figure; imshow(Gmag>0.1)
%thold = 0.01:0.01:0.5;

%thold = [0.1 0.2 0.3 0.4 0.5];
%thold = [0.5 0.6 0.7 0.8 1];
thold = 0.01:0.01:0.15;
%sigma = [0.6];
sigma = [0 0.5 1 2];

for j=1:length(sigma)
    
    if (sigma(j)~=0)
        g=fspecial('gaussian',[10 10],sigma(j));
        [Gmag, Gdir] = imgradient(rgb2gray(im2double(imfilter(im,g))));
    else
        [Gmag, Gdir] = imgradient(rgb2gray(im2double(im)));
    end
  
    Gmagnonmax=nonmax(Gmag, Gdir);
   % figure; imshow(Gmagnonmax);
for i=1:length(thold)
    edgeim=Gmagnonmax>thold(i);
    % non maxima supp
%    Gdir=-Gdir;

  %  Gdir(Gdir < 0) = 360 + Gdir(Gdir < 0);
  
  %  figure; imshow(edgeim);
%     Gdir(Gdir<0)=-Gdir(Gdir<0)+90;
%     edgeim=nonmax(edgeim,Gdir*(pi/180));
   
    % distance transform
    E(:,:,i,j)= bwdist(edgeim);
    % smoothing edges
    E(:,:,i,j)=E(:,:,i,j)./(E(:,:,i,j)+k);
end
end

E = reshape(E,size(im,1),size(im,2),length(thold)*length(sigma));

DT = mean(E,3);

end



