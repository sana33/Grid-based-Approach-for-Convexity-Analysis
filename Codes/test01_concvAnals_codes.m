%%%%%%%%%%%%% Efficacy test

close all
Eps = .05; q = 2; t = 2; grdSmpRat = .5;
[cncvStat,~] = concvAnals(ring,Eps,q,t,grdSmpRat)

close all
Eps = .05; q = 2; t = 2; grdSmpRat = 1;
[cncvStat,~] = concvAnals(ringHalf,Eps,q,t,grdSmpRat)

close all
Eps = .05; q = 3; t = 2; grdSmpRat = .5;
[cncvStat,~] = concvAnals(ring3D_half,Eps,q,t,grdSmpRat)

close all
Eps = .02; q = 2; t = 2; grdSmpRat = .5;
[cncvStat,~] = concvAnals(ring3D_half,Eps,q,t,grdSmpRat)

close all
Eps = .05; q = 2; t = 2; grdSmpRat = .5;
[cncvStat,~] = concvAnals(ring3D_half,Eps,q,t,grdSmpRat)


%%%%%%%%%%%%%%%% Single point concavity test

close all
Eps = .05; q = 2; t = 2; grdSmpRat = 1;
[cncvStat,~] = concvAnals([1 1],Eps,q,t,grdSmpRat)


%%%%%%%%%%%%%%%% test on grid accuracy

% close all
Eps = .005; q = 2; t = 2; grdSmpRat = 1;
[cncvStat,~] = concvAnals(rh1,Eps,q,t,grdSmpRat)

Eps = .05; q = 2; t = 2; grdSmpRat = 1;
[cncvStat,~] = concvAnals(rh1,Eps,q,t,grdSmpRat)

Eps = .2; q = 2; t = 2; grdSmpRat = 1;
[cncvStat,~] = concvAnals(rh1,Eps,q,t,grdSmpRat)

Eps = .5; q = 2; t = 2; grdSmpRat = 1;
[cncvStat,~] = concvAnals(rh1,Eps,q,t,grdSmpRat)


%%%%%%%%%%%%%%%% test on random sampling

% close all
Eps = .05; q = 2; t = 2; grdSmpRat = .05;
[cncvStat,~] = concvAnals(ring,Eps,q,t,grdSmpRat)

Eps = .05; q = 2; t = 2; grdSmpRat = .1;
[cncvStat,~] = concvAnals(ring,Eps,q,t,grdSmpRat)

Eps = .05; q = 2; t = 2; grdSmpRat = .15;
[cncvStat,~] = concvAnals(ring,Eps,q,t,grdSmpRat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% After disregarding Dimensionality Reduction and Random Sampling

close all
tStart = tic;
Eps = .05; grdSmpRat = 1;
[cncvStat,~] = concvAnals(ring,Eps,grdSmpRat)
toc(tStart)

close all
tStart = tic;
Eps = .05; grdSmpRat = .5;
[cncvStat,~] = concvAnals(ringHalf,Eps,grdSmpRat)
toc(tStart)

%% Evaluating Yeast data
close all

% Evaluating the whole data without the 6th(pox) feature
[~,score,~] = pca(YeastMod);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;

epsilon = .6; MinPts = 20; sizLim = 1e3;
[idx,~] = DBSCAN(YeastMod,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = .6; grdSmpRat = .5;
[cncvStat,~] = concvAnals(YeastMod,Eps,grdSmpRat)
toc(tStart)


% Evaluating the outliers containing data
W = X(y==0,:);
[~,score,~] = pca(W);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;

epsilon = .6; MinPts = 20; sizLim = 1e3;
[idx,~] = DBSCAN(W,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = .6; grdSmpRat = .5;
[cncvStat,~] = concvAnals(W,Eps,grdSmpRat)
toc(tStart)

% Evaluating single classes, disregarding outliers containing data
% Evaluating yeast, CYT class with 463 instances and 8 attributes
[~,score,~] = pca(yeast_CYT);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;
title('yeast, CYT class (463 by 8)');

epsilon = .6; MinPts = 10; sizLim = 1e3;
[idx,~] = DBSCAN(yeast_CYT,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = .6; grdSmpRat = .5;
[cncvStat,~] = concvAnals(yeast_CYT,Eps,grdSmpRat)
toc(tStart)

% Evaluating yeast, NUC class with 429 instances and 8 attributes
[~,score,~] = pca(yeast_NUC);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;
title('yeast, NUC class (429 by 8)');

epsilon = .6; MinPts = 10; sizLim = 1e3;
[idx,~] = DBSCAN(yeast_NUC,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = .6; grdSmpRat = .5;
[cncvStat,~] = concvAnals(yeast_NUC,Eps,grdSmpRat)
toc(tStart)

% Evaluating yeast, MIT class with 244 instances and 8 attributes
[~,score,~] = pca(yeast_MIT);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;
title('yeast, MIT class (244 by 8)');

epsilon = .6; MinPts = 10; sizLim = 1e3;
[idx,~] = DBSCAN(yeast_MIT,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = .6; grdSmpRat = .5;
[cncvStat,~] = concvAnals(yeast_MIT,Eps,grdSmpRat)
toc(tStart)

%% Evaluating Mammography
close all

% Evaluating the outliers containing data
W = X(y==0,:);
[~,score,~] = pca(W);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;

epsilon = 16; MinPts = 20; sizLim = 1e3;
[idx,~] = DBSCAN(W,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = 16; grdSmpRat = .5;
[cncvStat,~] = concvAnals(W,Eps,grdSmpRat)
toc(tStart)

% Evaluating single classes, disregarding outliers containing data
% Evaluating mamography, non-calcification class with 10923 instances and 6 attributes
[~,score,~] = pca(mamography_1);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;
title('mamography, non-calcification class (10923 by 6)');

epsilon = 16; MinPts = 10; sizLim = 1e3;
[idx,~] = DBSCAN(mamography_1,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = 16; grdSmpRat = .5;
[cncvStat,~] = concvAnals(mamography_1,Eps,grdSmpRat)
toc(tStart)

% Evaluating mamography, calcification class with 260 instances and 6 attributes
[~,score,~] = pca(mamography_2);
figure; gscatter(score(:,1),score(:,2),ones(size(score,1),1),'b'); grid on;
title('mamography, calcification class (260 by 6)');

epsilon = 13; MinPts = 10; sizLim = 1e3;
[idx,~] = DBSCAN(mamography_2,epsilon,MinPts,sizLim);
unique(idx)
sum(idx==0)

tStart = tic;
Eps = 13; grdSmpRat = .5;
[cncvStat,~] = concvAnals(mamography_2,Eps,grdSmpRat)
toc(tStart)



