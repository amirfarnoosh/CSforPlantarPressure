function [err,w,pv] = gmmError(pd, gmm, locs,v)
    if(nargin < 3)
        pdLocs = pd;
    else
        pdLocs = locs;
    end
    
    [gmm,w] = adjustWeights(pdLocs, gmm,v);
    pv = getPointValues(pd, gmm,v,w);
    pv(pv<0)=0;
    err = abs(pv) - pd(:,3);
    
    if(any(isnan(err)))
        [gmm.w]
    end
    err = err.*err;
end


function pv = getPointValues(pd, gmm,v,w)
    pv = zeros(size(pd,1), numel(gmm));
    
    for idx=1:numel(gmm)
        pv(:,idx) = mvnpdf(pd(:,1:2), gmm(idx).mu, gmm(idx).sigma);
    end
    
    weights = w;%[gmm.w];
    weights(isnan(weights)) = 0;
    
    pv = pv * v* weights;
end

% use least squared to determine ideal weights
function [gmm,w] = adjustWeights(pd, gmm,v)
    A = zeros(size(pd,1), numel(gmm));
    
    for idx=1:numel(gmm)
        A(:,idx) = mvnpdf(pd(:,1:2), gmm(idx).mu, gmm(idx).sigma);
    end

    x = pd(:,3);
    [w flag] = lsqr(A*v,x);
    
%     for idx = 1:numel(gmm)
%         gmm(idx).w = w(idx);
%     end

end
