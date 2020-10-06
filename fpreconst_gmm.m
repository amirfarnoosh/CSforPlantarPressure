function [LSE,dist_err,RImg,W]=fpreconst_gmm(data,subi,legi,method,kbest,K,varargin)

nVararg=length(varargin);
%varargin=[nframe,prog_flag,blkSize]
fig_flag=0;

if isempty(data)
    load('Dataset\Data');
end

% subi %subject number
% legi %leg number 1 for left and 2 for right

% K=5; %number of random sensors;
% kbest.K=1; %number of sensors in each pre-defined region
% kbest.mode %'peaks' & 'rnd'
% method.name %'Laplace', 'nearest', 'linear', 'natural'

if nVararg<3
    blkSize=[20,37]; %image size or block size
else
    blkSize=varargin{3};
end
if nVararg<2
    prog_flag=0; %if 1, shows progress
else
    prog_flag=varargin{2};
end
if nVararg<1
    nframe=size(data.datatest{subi,legi},2); %number of test frames to consider
elseif isempty(varargin{1})
    nframe=size(data.datatest{subi,legi},2);
else
    nframe=varargin{1};
end

%binary selection matrix formation for dictionary

sx=data.sx(:,(subi-1)*2+legi);
sy=data.sy(:,(subi-1)*2+legi);


% regionIdxs=round(reshape([data.gmm{subi,legi}.mu],2,40));
% regionIdxs=sub2ind([44,52],regionIdxs(1,:)',regionIdxs(2,:)');
% 
% if kbest.mode=='rnd'
%     [regionIdxs1]=pickKsensor(subi,legi,kbest.K,data.cfg);
%     regionIdxs1=regionIdxs1(:);
% elseif kbest.mode=='peaks'
%     regionIdxs1=data.regionIdxs{subi,legi}(:,1:kbest.K);
%     regionIdxs1=regionIdxs1(:);
% end
% 
% org_img=zeros([44,52]);
% org_img(regionIdxs1(:))=1;
% norg_img=org_img(sx(1):sx(2),sy(1):sy(2));
% indxs1=find(norg_img==1)';
% 
% sen_indx=[];
% 
% org_img=zeros([44,52]);
% org_img(regionIdxs(:))=1;
% norg_img=org_img(sx(1):sx(2),sy(1):sy(2));
% indxs=find(norg_img==1)';
% sen_indx=[sen_indx,indxs];
% rnd=randperm(length(sen_indx));
% sen_indx=sen_indx(rnd);
% 
% sen_indx=[indxs1,sen_indx];
% sen_indx=sen_indx(1:kbest.K*4+K);


if kbest.mode=='rnd'
    [regionIdxs]=pickKsensor(subi,legi,kbest.K,data.cfg);
    regionIdxs=regionIdxs(:);
elseif kbest.mode=='peaks'
    regionIdxs=data.regionIdxs{subi,legi}(:,1:kbest.K);
    regionIdxs=regionIdxs(:);
end

org_img=zeros([44,52]);
org_img(regionIdxs(:))=1;
norg_img=org_img(sx(1):sx(2),sy(1):sy(2));
indxs=find(norg_img==1)';

sen_indx=[];
sen_indx=[sen_indx,indxs];


zrm=find(data.zeromask(:,(subi-1)*2+legi)==1)';
for i=1:length(sen_indx)
    zrm(zrm==sen_indx(i))=[];
end
rnd=randperm(length(zrm));
rndd=zrm(rnd(1:K));
sen_indx=[sen_indx,rndd];


zrmask=data.zeromask(:,(subi-1)*2+legi);

%sen_indx=find(zrmask==1)'; %consider commenting this
%error per distance%
dist_err=[];
% [xd,yd]=ind2sub([20,37],(1:20*37)');
% [xs,ys]=ind2sub([20,37],sen_indx);
% dist=repmat([xd(zrmask~=0),yd(zrmask~=0)],1,1,kbest.K*4+K)-...
%     repmat(reshape([xs;ys],1,2,kbest.K*4+K),length(find(zrmask~=0)),1,1);
% dist=(min(sum(dist.^2,2),[],3)).^0.5;
%%%%%%%%%%%%%%%

%reconstruction phase

if prog_flag==1
    textprogressbar('Reconstructing: ');
end

lsey=0;
RImg=[];
cnt=0;
W=[];
for icnt=1:1:nframe
    
    I=data.datatest{subi,legi}(:,icnt);
    Iot=zeros(20,37);
    if I==0
        %cnt=cnt+1;
        %continue
    end
    
   
    zrmask_org=zeros([44,52]);
    zrmask_org(sx(1):sx(2),sy(1):sy(2))=reshape(zrmask,20,37);
    [x,y]=find(zrmask_org==1);
    pd(:,1)=x;pd(:,2)=y;pd(:,3)=I(zrmask==1);
    
    sen_img=zeros([20,37]);
    sen_img(sen_indx)=1;
    sen_org=zeros([44,52]);
    sen_org(sx(1):sx(2),sy(1):sy(2))=sen_img;
    [x,y]=find(sen_org==1);
    locs(:,1)=x;locs(:,2)=y;locs(:,3)=I(sen_img==1);
    
    [lse,w,vp] = gmmError(pd, data.gmm{subi,legi}, locs,data.V{subi,legi}(:,1:5));
    W=[W,w];
    Iot(zrmask==1)=vp;
    
    if fig_flag==1
        fig=figure('Visible','off');
        set(fig,'WindowStyle','docked');
        
        surf(reshape(I,20,37));
        saveas(gcf,['Images\img_a' num2str(icnt) '.fig']);
        surf(Iot);
        saveas(gcf,['Images\img_gmm' num2str(icnt) '.fig']);
            
%         %map=jet(225);%colormap('jet');
%         subplot(2,1,1);
%         surf(reshape(I,20,37));
%         zlim([0,250]);
%         %colormap(map);
%         caxis([0,180]);
%         set(gcf,'color','black');
%         ax=gca;
%         ax.Color='black';
%         view([13.43 76.93])
%         subplot(2,1,2);
%         pic=surf(Iot);
%         zlim([0,250]);
%         %colormap(map);
%         caxis([0,180]);
%         set(gcf,'color','black');
%         ax=gca;
%         ax.Color='black';
%         view([13.43 76.93]);
%         fig=gcf;
%         fig.InvertHardcopy = 'off';
%         saveas(pic,['Images\img' num2str(icnt) '.jpg']);
    end
    if kbest.option==1
        RImg=[RImg;I(zrmask~=0),Iot(zrmask~=0)];
    end
        
    lsey=lsey+lse;
    
    if prog_flag==1
        textprogressbar(icnt/nframe*100);
    end

end

%err=(lsey/(nframe-cnt)).^0.5;
%dist_err=[dist,err];
LSE=(sum(lsey)/(nframe-cnt)/length(find(data.zeromask(:,(subi-1)*2+legi)==1)))^0.5;

if prog_flag==1
    textprogressbar('done');
end