%% Setup basic parameters:

landmarks=[8333,7301,6011,9365,10655,8320,8311,8286,8275,5959,4675,3642,4922,3631,2088,27391,39839,40091,40351,6713,10603,11641,12673,11244,12661,14472,27778,41804,41578,41310,9806,8345,7442,6936,5392,7335,7851,8354,9248,9398,11326,9399,9129,9406,9128,8890,8367,7858,7580,7471,8374,23002,32033]';
%landmarks=[20896, 2088, 5394, 8320, 14472, 11328, 4646]';
%load('01_MorphableModel.mat')
load('BFMedgestruct.mat')
[tri,pts]=plyread('C:\Users\Dell\Documents\MATLAB\faces\00001_20061015_00418_neutral_face05.ply','tri');
FV.vertices=pts;
FV.faces=tri;

%data=createRenderings('C:\Users\Dell\Documents\MATLAB\face-release1.0-basic\renderings\');

% figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light; axis equal; axis off
% hold on
% plot3(FV.vertices(landmarks,1),FV.vertices(landmarks,2),FV.vertices(landmarks,3),'xr')


%% Landmark Detection
%im = imread('images/9.jpg');
%im = imread('renderings/00001_20061015_00418_neutral_0_0.png');
im = imread('pie.png');
xp= RamananDetector(im);
%i=6;

%im = data(i).im;
%xp =[data(i).xp(2,:); data(i).xp(1,:)];
%xp =data(i).xp;
%xp(2,:) = size(im,1)+1-xp(2,:);
%landmarks= data(i).landmarks;


%% Display fitted model using xp :
ndims = 75; %20Number of model dimensions
numsd=2;    %Number of standard deviations

[b,R_est,t_est,s_est] = FitEdges( im,xp,landmarks,shapePC,shapeMU,shapeEV,Ef,Ev,tri );
%[b,R_est,t_est,s_est] = FitSingleSOP( xp,shapePC,shapeMU,shapeEV,ndims,numsd,landmarks);
FV2.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,53490)';
FV2.faces = tri;
FV2.facevertexcdata = min(1,max(0,double(reshape(texMU,3,53490)')./255));
figure; patch(FV2, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceLighting', 'phong'); axis equal; axis off;

%% Plot projected 3D landmarks
%xp(2,:) = size(im,1)+1-xp(2,:);
% figure; imshow(im);
% hold
% plot(xp(1,:),size(im,1)+1-xp(2,:),'xb');
% 
% pts=FV2.vertices(landmarks,:)';
% rotpts=R_est*pts;
% xp2 = [s_est.*(rotpts(1,:)+t_est(1)); s_est.*(rotpts(2,:)+t_est(2))];
% xp2(2,:) = size(im,1)+1-xp2(2,:);
% plot(xp2(1,:),xp2(2,:),'xr');
% 
% 
% error = mean(sqrt(sum((xp2-xp).^2,1)));
% title(['Error rate: ',num2str(error)]);
% xlabel({
%      [' Number of standard deviations: ',num2str(numsd)];
%      [' Number of model dimensions: ', num2str(ndims)] }');
% hold off;

 %% Texture sampling
im = double(im)./255;
pts=FV2.vertices';
rotpts=R_est*pts;
pts2D = [s_est.*(rotpts(1,:)+t_est(1)); s_est.*(rotpts(2,:)+t_est(2))];
%figure; plot(pts2D(1,:),pts2D(2,:),'.')

[rows,cols]=meshgrid(1:size(im,2),1:size(im,1));

for col=1:3
    FV2.facevertexcdata(:,col) = interp2(rows,cols,im(:,:,col),pts2D(1,:),size(im,1)+1-pts2D(2,:));
end
% 
% % for col= 1:53490
% % Vertices(col) = col;
% % end
% 
 nonVisibleVertices= setdiff(1:53490, visiblevertices(FV2, R_est));
% 
for col=1:3
FV2.facevertexcdata(nonVisibleVertices,col) = NaN;
end
% 
% 
% 
% %intersect(occludingVertices,visiblevertices(FV, eye(3)));
% %FV2.facevertexcdata(visiblevertices(FV2, eye(3)),2) = NaN;
% %FV2.facevertexcdata(visiblevertices(FV2, R_est),2) = NaN;
% 
% %intersect(FV2.facevertexcdata,visiblevertices(FV2, R_est));
% 
% %burada remove ediyoz col yerine landmark pointleri yazdýkmý tamam
% % for col=1:
% % FV2.facevertexcdata(col,2) = NaN;
% % end
% 
% 
 figure; patch(FV2, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceLighting', 'phong'); axis equal; axis off;
