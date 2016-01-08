%% Setup basic parameters:

landmarks=[8333,7301,6011,9365,10655,8320,8311,8286,8275,5959,4675,3642,4922,3631,2088,27391,39839,40091,40351,6713,10603,11641,12673,11244,12661,14472,27778,41804,41578,41310,9806,8345,7442,6936,5392,7335,7851,8354,9248,9398,11326,9399,9129,9406,9128,8890,8367,7858,7580,7471,8374,23002,32033]';
%load('01_MorphableModel.mat')
load('BFMedgestruct.mat')
[tri,pts]=plyread('C:\Users\Dell\Documents\MATLAB\faces\00001_20061015_00418_neutral_face05.ply','tri');
FV.vertices=pts;
FV.faces=tri;

data=createRenderings('C:\Users\Dell\Documents\MATLAB\face-release1.0-basic\renderings\');

%% Landmark Detection
im = imread('pie.png');
xp= RamananDetector(im);

%% Sample Selection
%i=8;
%im = data(i).im;
%xp =data(i).xp;
%xp(2,:) = size(im,1)+1-xp(2,:);
%landmarks= data(i).landmarks;

%% Initialization
ndims=5;
[b,R,t,s] = FitSingleSOP( xp,shapePC,shapeMU,shapeEV,ndims,landmarks);
FV.vertices=reshape(shapePC(:,1:ndims)*b+shapeMU,3,size(shapePC,1)/3)';
figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light; axis equal; axis off
% FV.vertices=reshape(shapePC(:,1:ndims)*b+shapeMU,3,size(shapePC,1)/3)';
% figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light; axis equal; axis off
%alterb=zeros(26,1);
%figure; imshow(im);
%edges = imresize(double(edge(imresize(rgb2gray(im),0.5))),[size(im,1) size(im,2)])>0.5;
%DT=bwdist(edges);

% b = [Rts2b(double(R),double(t),double(s)); b];

%% Distance Transform
DT = GetDT( im,60);
tolerance = [1e-1 1e-2 1e-3 1e-4 1e-5];
fina = [1e-3 1e-3 1e-3 1e-4 1e-4];
 figure; imshow(DT,[]);
% hold on;
for i=1:5

%alterb(:,i)=b;
[ occludingVertices ] = occludingBoundaryVertices( FV,Ef,Ev,R );

figure; imshow(DT,[]);
hold on;

% Get interpolation
%options = optimset('Display','iter');
%options = optimset('Display','iter','TolFun',1e-9);
%options = optimoptions('fmincon','Display','iter','TolFun',1e-20,'FinDiffRelStep',1e-3);

%b2 = lsqnonlin(@(b) RomdhaniObjective( b,shapeMU,shapePC(:,1:ndims),shapeEV(1:ndims),R,t,s,occludingVertices,edges,DT),b,-3.*double(shapeEV),3.*double(shapeEV),options);
%b2 = fminunc(@(b) RomdhaniObjective( b,shapeMU,shapePC(:,1:ndims),shapeEV(1:ndims),R,t,s,occludingVertices,edges,DT),b,options);
%b = lsqnonlin(@(b) RomdhaniObjective( b,shapeMU,shapePC(:,1:ndims),occludingVertices,edges,DT ),double(b0),-3.*double(shapeEV),3.*double(shapeEV),options);

%TolFun: Termination tolerance on the function value. 1e-20 default: 1e-6 
%FinDiffRelStep: step size factor for finite differences. 1e-3 -1 -9
%options = optimoptions('TolFun',1e-6,'Display','iter','FinDiffRelStep',1e-9); %,'LargeScale','off');
%options = optimoptions('fmincon','Display','iter','TolFun',1e-20,'FinDiffRelStep',1e-1);
%b = fminunc(@(b) RomdhaniObjectiveShape(b,shapeMU,shapePC(:,1:ndims),R,t,s,occludingVertices,edges,DT),double(b),options);
%b = fmincon(@(b) RomdhaniObjectiveShape(b,shapeMU,shapePC(:,1:ndims),R,t,s,occludingVertices,edges,DT),double(b),[],[],[],[],-5.*double(shapeEV),5.*double(shapeEV),[],options);
%options = optimset('Display','iter','TolFun',1e-3,'FinDiffRelStep',1e-3);
%options = optimset('Display','iter','MaxIter',30,'FinDiffType','central');

%options = optimset('Display','iter','TolFun',1e-9,'FinDiffRelStep',1e-6); %'TolFun',1e-6,'FinDiffRelStep',1e-3,
%b = lsqnonlin(@(b) RomdhaniObjectiveShape( b,xp,shapeEV(1:ndims),shapeMU,shapePC(:,1:ndims),R,t,s,occludingVertices,DT,landmarks ),b,[],[],options);
options = optimoptions('fmincon','Display','iter','TolFun',tolerance(i),'FinDiffRelStep',fina(i));
%options = optimoptions('fmincon','Display','iter','TolFun',1e-6,'TolX',1e-6,'FinDiffRelStep',1e-9);
%options = optimoptions('fmincon','Display','iter');
b = fmincon(@(b) RomdhaniObjectiveShape( b,xp,shapeEV(1:ndims),shapeMU,shapePC(:,1:ndims),R,t,s,occludingVertices,DT,landmarks ),double(b),[],[],[],[],-3.*double(shapeEV(1:ndims)),3.*double(shapeEV(1:ndims)),[],options);

%options = optimoptions('fminunc','Display','iter','TolFun',1e-6,'FinDiffRelStep',1e-1);
%b = fminunc(@(b) RomdhaniObjectiveShape( b,xp,shapeEV(1:ndims),shapeMU,shapePC(:,1:ndims),R,t,s,occludingVertices,DT,landmarks ),double(b),options);




n = 2;
pause(n);

Rts= Rts2b(double(R),double(t),double(s));
Rts(4:5)=Rts(4:5)./1e+7;
%Rtsoptions = optimoptions('fminunc','Display','iter','TolFun',1e-6,'FinDiffRelStep',1e-1);
Rtsoptions = optimoptions('fminunc','Display','iter','TolFun',tolerance(i),'FinDiffRelStep',fina(i));
Rts = fminunc(@(Rts) RomdhaniObjectivePose( Rts,b,shapeMU,shapePC(:,1:ndims),occludingVertices,DT,landmarks,xp ),Rts,Rtsoptions);
%Rts = lsqnonlin(@(Rts) RomdhaniObjectivePose( Rts,b,shapeMU,shapePC(:,1:ndims),occludingVertices,DT ),Rts,[],[],Rtsoptions);
Rts(4:5)=Rts(4:5).*1e+7;
[R,t,s]=b2Rts(Rts);

% options = optimset('TolFun',1e-6,'Display','iter'); %,'LargeScale','off'); %Sum Dists cancelled
% b(7:end)=b(7:end)./shapeEV(1:ndims);
% b = lsqnonlin(@(b) RomdhaniObjective( b,shapeMU,shapePC(:,1:ndims),shapeEV(1:ndims),occludingVertices,edges,DT ),b,[],[],options);
% %b = fminunc(@(b) RomdhaniObjective( b,shapeMU,shapePC(:,1:ndims),shapeEV(1:ndims),occludingVertices,edges,DT ),b,options);
% b(7:end)=b(7:end).*shapeEV(1:ndims);
% [R,t,s]=b2Rts(b(1:6));
% b = b(7:end);

FV.vertices=reshape(shapePC(:,1:ndims)*b+shapeMU,3,size(shapePC,1)/3)';
%figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light; axis equal; axis off
end

% %% TEST
 figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light; axis equal; axis off
% FV2.faces=tri;
% FV2.vertices=reshape(shapePC(:,1:ndims)*b2+shapeMU,3,size(shapePC,1)/3)';
% figure; patch(FV2, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light; axis equal; axis off


%% Landmark Difference
% figure;
% lndmrkpts=FV.vertices(occludingVertices,:)';
% plot(lndmrkpts(1,:),lndmrkpts(2,:),'.r');
% hold on;
% romdhanipts=FV2.vertices(occludingVertices,:)';
% plot(romdhanipts(1,:),romdhanipts(2,:),'.b');

%% Texture
%im = imread('images/1.jpg');
im = imread('pie.png');
im = double(im)./255;
rotpts=R*FV.vertices';
romdhani2d = [s.*(rotpts(1,:)+t(1)); s.*(rotpts(2,:)+t(2))];
[rows,cols]=meshgrid(1:size(im,2),1:size(im,1));
FV.facevertexcdata = zeros(0);
for col=1:3
    FV.facevertexcdata(:,col) = interp2(rows,cols,im(:,:,col),romdhani2d(1,:),size(im,1)+1-romdhani2d(2,:));
end
nonVisibleVertices= setdiff(1:53490, visiblevertices(FV, R));
for col=1:3
FV.facevertexcdata(nonVisibleVertices,col) = NaN;
end
figure; patch(FV, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceLighting', 'phong'); axis equal; axis off;
