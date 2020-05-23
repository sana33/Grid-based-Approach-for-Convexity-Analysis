% 
% Copyright (c) 2018
% 
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% 
% Developer: Sayyed Ahmad Naghavi Nozad (SANN)
% 
% Contact Info: ahmad.naghavi.aut@gmail.com
% 

function [idx,isNoise] = DBSCAN(X,epsilon,MinPts,sizLim)

n = size(X,1);
idx = sparse(n,1);
isNoise = sparse(n,1);
clustNo = 0;
neighbMat_Maker();
idxCore = find(sum(neighbMatX,2)>=MinPts);
corPntNo = numel(idxCore);
visited = sparse(n,1);

for d1 = 1:corPntNo
    if ~visited(idxCore(d1))
        clustNo = clustNo+1;
        clust = idxCore(d1);
        idxCorConsd = [];
        while true
            clustTemp = clust;
            clust = union(clust,find(any(neighbMatX(intersect(clust,setdiff(idxCore,idxCorConsd)),:),1)));
            if isempty(setdiff(clust,clustTemp))
                break
            end
            idxCorConsd = intersect(idxCore,clustTemp);
        end
        idx(clust) = clustNo;
        visited(clust) = 1;
    end
end

isNoise(idx==0) = 1;
idx = full(idx);

% Nested function for computing the neighborhood matrix
    function neighbMat_Maker()
        
        neighbMatX = logical([]);
        
        cntLim = ceil(n/sizLim);
        
        for c1 = 1:cntLim
            %     c1
            if c1 ~= cntLim
                Y1 = X((c1-1)*sizLim+1:c1*sizLim,:);
            else
                Y1 = X((c1-1)*sizLim+1:end,:);
            end
            neighbMatY = logical([]);
            
            for c2 = 1:cntLim
                %         c2
                if c2 ~= cntLim
                    indVec = (c2-1)*sizLim+1:c2*sizLim;
                else
                    indVec = (c2-1)*sizLim+1:n;
                end
                Y2 = X(indVec,:);
                neighbMatY = [neighbMatY pdist2(Y1,Y2)<=epsilon];
            end
            
            neighbMatX = [neighbMatX; neighbMatY];
        end
        
    end

end



