function [cncvStat,concvEvdCands] = concvAnals(X,Eps,eta)
% convxAnals - Analyzes the convexity of a density-based cluster.
% 
% This function is an implementation of the paper named "A Grid-based
% Approach for Convexity Analysis of a Density-based Cluster" by
% Sayyed-Ahmad Naghavi-Nozad et al
% 
% Syntax:  [cncvStat,concvEvdCands] = convxAnals(X,Eps,eta)
% 
% Inputs:
%     X - Input dataset
%     Eps - Grid accuracy
%     eta - Random sampling rate for the grid
% 
% Outputs:
%     cncvStat - Concavity status of a cluster
%     concvEvdCands - Candidate points for proving the concavity
% 
% Example: 
%     [cncvStat,~] = convxAnals(ring,.05,.5);
% 
% Other m-files required: none
% Subfunctions: yes
% MAT-files required: yes
% 
% 
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Dept. of Computer Engineering and
% Information Technology
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ceit.aut.ac.ir/~sann_cv/
% June 2019; Last revision: 31-May-2019

%------------- BEGIN CODE --------------

%% Phase 0: Initializing the grid
p = size(X,2);
Xmin = min(X,[],1);
Xmax = max(X,[],1);

Xgv = cell(0);
sizVec = zeros(1,p);
for c1 = 1:p
    Xgv{c1} = Xmin(c1)-2*Eps:Eps:Xmax(c1)+2*Eps;
    sizVec(c1) = numel(Xgv{c1});
end

[gridPts,gridSampLinIdx,~] = cartesProd(Xgv,sizVec,eta);

%% Detecting non-neighboring grid points
[~,GXdist] = knnsearch(X,gridPts);
clstGridEpsNghbIdxLog = GXdist<=Eps;
nonClstGridEpsNghbIdx = gridSampLinIdx(~clstGridEpsNghbIdxLog);

%% Finding marginal cluster points
[nonClstGridNghbIdx] = gridNghbFind(nonClstGridEpsNghbIdx,sizVec);
margGridPtsIdx = setdiff(nonClstGridNghbIdx,nonClstGridEpsNghbIdx);
[margGridPts,~] = gridLinIdx2Pts(Xgv,sizVec,margGridPtsIdx);

[margPtsIdx,mGXdist] = knnsearch(X,margGridPts);
margGridPtsIdxLog = mGXdist<=Eps;
margGridPtsIdx = margGridPtsIdx(margGridPtsIdxLog);
margGridPts = margGridPts(margGridPtsIdxLog,:);

margPtsIdx = unique(margPtsIdx);
margPts = X(margPtsIdx,:);
margPtsNo = size(margPts,1);

%% Analyzing midpoint convexity
cncvStat = 0;
concvEvdCands = [];
for c1 = 1:margPtsNo
    means = (margPts(c1,:)+margPts(c1+1:margPtsNo,:))./2;
    [~,muXdist] = knnsearch(X,means);
    concvEvdCandsIdx = find(muXdist>Eps);
    
    if ~isempty(concvEvdCandsIdx)
        cncvStat = 1;
        concvEvdCands = means(concvEvdCandsIdx,:);
        break;
    end
end

%% Plotting 2-D results
if p==2
    figure;
    
    subplot(2,3,1);
    gscatter(X(:,1),X(:,2),ones(size(X,1),1),'b'); grid on; hold on;
    xlim([Xmin(1)-2*Eps Xmax(1)+2*Eps]); ylim([Xmin(2)-2*Eps Xmax(2)+2*Eps]);
    title('(a) Ring-Shaped Cluster'); legend('off');
    
    subplot(2,3,2);
    gscatter(X(:,1),X(:,2),ones(size(X,1),1),'b'); grid on; hold on;
    gscatter(gridPts(:,1),gridPts(:,2),zeros(size(gridPts,1),1),'r');
    xlim([Xmin(1)-2*Eps Xmax(1)+2*Eps]); ylim([Xmin(2)-2*Eps Xmax(2)+2*Eps]);
    title('(b) Sampled Grid Points'); legend('off');
    
    subplot(2,3,3);
    gscatter(X(:,1),X(:,2),ones(size(X,1),1),'b'); grid on; hold on;
%     gscatter(gridPts(clstGridEpsNghbIdxLog,1),gridPts(clstGridEpsNghbIdxLog,2),zeros(sum(clstGridEpsNghbIdxLog),1),'r');
    gscatter(gridPts(~clstGridEpsNghbIdxLog,1),gridPts(~clstGridEpsNghbIdxLog,2),zeros(sum(~clstGridEpsNghbIdxLog),1),'g');
    xlim([Xmin(1)-2*Eps Xmax(1)+2*Eps]); ylim([Xmin(2)-2*Eps Xmax(2)+2*Eps]);
    title('(c) Non-Neighboring Grid Points'); legend('off');
    
    subplot(2,3,4);
    gscatter(X(:,1),X(:,2),ones(size(X,1),1),'b'); grid on; hold on;
    gscatter(margGridPts(:,1),margGridPts(:,2),zeros(length(margGridPtsIdx),1),'k');
    xlim([Xmin(1)-2*Eps Xmax(1)+2*Eps]); ylim([Xmin(2)-2*Eps Xmax(2)+2*Eps]);
    title('(d) Marginal Grid Points'); legend('off');
    
    subplot(2,3,5);
    gscatter(X(:,1),X(:,2),ones(size(X,1),1),'b'); grid on; hold on;
    gscatter(X(margPtsIdx,1),X(margPtsIdx,2),zeros(length(margPtsIdx),1),'k');
    xlim([Xmin(1)-2*Eps Xmax(1)+2*Eps]); ylim([Xmin(2)-2*Eps Xmax(2)+2*Eps]);
    title('(e) Marginal Cluster Points'); legend('off');
    
    if ~isempty(concvEvdCandsIdx)
        subplot(2,3,6);
        gscatter(X(:,1),X(:,2),ones(size(X,1),1),'b'); grid on; hold on;
        plot(margPts([c1,c1+concvEvdCandsIdx(1)],1),margPts([c1,c1+concvEvdCandsIdx(1)],2),'^k','MarkerFaceColor','k','LineStyle','none');
        plot(concvEvdCands(1,1),concvEvdCands(1,2),'sr','MarkerFaceColor','r','LineStyle','none');
        fh = @(x,y) (x-concvEvdCands(1,1)).^2+(y-concvEvdCands(1,2)).^2-Eps^2; ezplot(fh,[Xmin(1) Xmax(1)]);
        xlim([Xmin(1)-2*Eps Xmax(1)+2*Eps]); ylim([Xmin(2)-2*Eps Xmax(2)+2*Eps]);
        title('(f) Midpoint Convexity Analysis'); legend('off');
    end
end

end


%% Subfunctions Here!

function [gridPts,gridSampLinIdx,gridSampSubIdx] = cartesProd(Xgv,sizVec,eta)
    
    grdCnt = prod(sizVec);
    gridSampLinIdx = 1:floor(1/eta):grdCnt;
%     gridSampLinIdx = randperm(grdCnt,ceil(eta*grdCnt));
    
    [gridPts,gridSampSubIdx] = gridLinIdx2Pts(Xgv,sizVec,gridSampLinIdx);

end

function [gridPts,subIdx] = gridLinIdx2Pts(Xgv,sizVec,linIdx)

p = numel(sizVec);
[subIdx{1:p}] = ind2sub(sizVec,linIdx);

for c1 = 1:p
    gridPts(:,c1) = Xgv{c1}(subIdx{c1});
end

end

function [nghbLinInd] = gridNghbFind(linInd,sizVec)

p = length(sizVec);
[SUB{1:p}] = ind2sub(sizVec,linInd);
subs = zeros(numel(linInd),p);
for c1 = 1:p
    subs(:,c1) = SUB{c1}(:);
end

nghbSub = [];
for c1 = 1:p
    nghbSub = [nghbSub
        subs(:,1:c1-1) subs(:,c1)-1 subs(:,c1+1:end)
        subs(:,1:c1-1) subs(:,c1)+1 subs(:,c1+1:end)];
end
nghbSub = nghbSub(all(nghbSub>=1 & nghbSub<=sizVec,2),:);
nghbSub = unique(nghbSub,'rows');

nghbLinInd = SUB2IND(sizVec,nghbSub);

end

function [linIdx] = SUB2IND(sizVec,nghbSub)

prodVec = cumprod(sizVec);
linIdx = sum((nghbSub(:,2:end)-1).*prodVec(1:end-1),2)+nghbSub(:,1);

end

%------------- END OF CODE --------------
