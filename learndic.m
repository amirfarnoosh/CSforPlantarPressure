function DCT=learndic(data,subi,legi,dict)

warning('off');
if isempty(data)
    load('Dataset\Data');
end

% subi %subject number
% legi %leg number 1 for left and 2 for right

% dict.type='learned'; % 'learned', 'DCT' and 'DWT'
% dict.name='DCTtn3'; %if type is 'learned' load this dictionary
% dict.learn=0; %set to 1 to learn dictionary
% dict.iter=5;
% dict.S=2;
% dict.blkSize=[20,37];
% dict.slidingDis=1;
% dict.disp=0;

D=data.datatrain{subi,legi};

s_level=dict.S;
blkSize=dict.blkSize;
slidingDis=dict.slidingDis;

blocks=[];
for i=1:size(D,2)
    blocks=[blocks,my_im2col(reshape(D(:,i),20,37),blkSize,slidingDis)];
end
blocks(:,~any(blocks,1))=[];
%size(blocks)

%initialize dictionary
DD=50;%20*37*1 %dic second dimension

%learning
param.errorFlag=dict.errorFlag;
param.errorGoal = dict.errorGoal;
param.totalerr=dict.totalerr;
param.L=s_level;
param.preserveDCAtom=0;
param.K =DD;
param.numIteration = dict.iter;
param.displayProgress=1;
%param.initialDictionary = DCT;
param.InitializationMethod =  'DataElements';
param.OMPFlag = 1; %solve with selection matrix
param.OMP = dict.DCTm; %solve with this selection matrix

if dict.nn==1
    [DCT,output] = KSVD_NN(blocks,param);
elseif dict.nn==0
    [DCT,output] = KSVD(blocks,param);
end

if dict.disp==1
    j1=20;
    img=[];
    imgt=[];
    cnt=0;
    for i=1:DD
        img=[img,reshape(DCT(:,i),[20,37])];
        cnt=cnt+1;
        if cnt==j1
            imgt=[imgt;img];
            img=[];
            cnt=0;
        end
    end
    imshow(imgt,[]);
end

if dict.nn==1
    save([dict.dir,'DCTt_' num2str(dict.S) '_' num2str(dict.iter) 'nn_sub' num2str(subi) num2str(legi)],'DCT');
elseif dict.nn==0
    save([dict.dir,'DCTt_' num2str(dict.S) '_' num2str(dict.iter) '_sub' num2str(subi) num2str(legi)],'DCT');
end