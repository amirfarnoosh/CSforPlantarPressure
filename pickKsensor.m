
function [regionIdxs]=pickKsensor(subi, legi, k,cfg)

expname=cfg.subjects(subi,1).src;

if legi==1
    lrFoot='left';
elseif legi==2
    lrFoot='right';
end

% region image
imfile = ['Mask\','Mask-', expname,'-', lrFoot,'-', 'out', '.png'];
regionData = imread(imfile);

%=== Iterate through regions ===

for region=1:numel(cfg.regions)
    rpos = find(regionData==region);
    rnd=randperm(length(rpos));
    regionIdxs(region, :)=rpos(rnd(1:k));
    
end
%---------------------------