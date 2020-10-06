function [LSE,dist_err,RImg]=fpreconst_interp(data,subi,legi,method,kbest,K,varargin)

dist_err=0;
nVararg=length(varargin);
%varargin=[nframe,prog_flag,blkSize]

if isempty(data)
    load('Dataset\Data');
end

% subi %subject number
% legi %leg number 1 for left and 2 for right

% K=5; %number of random sensors;
% kbest.K=1; %number of sensors in each pre-defined region
% kbest.mode %'peaks' & 'rnd'
% method.name %'Laplace', 'nearest', 'linear', 'natural'
fig_flag=0;

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

if kbest.mode=='rnd'
    [regionIdxs]=pickKsensor(subi,legi,kbest.K,data.cfg);
    regionIdxs=regionIdxs(:);
elseif kbest.mode=='peaks'
    regionIdxs=data.regionIdxs{subi,legi}(:,1:kbest.K);
    regionIdxs=regionIdxs(:);
end

sen_indx=[];

org_img=zeros([44,52]);
org_img(regionIdxs(:))=1;
norg_img=org_img(sx(1):sx(2),sy(1):sy(2));
indxs=find(norg_img==1)';
sen_indx=[sen_indx,indxs];

if K~=0
    zrm=find(data.zeromask(:,(subi-1)*2+legi)==1)';
    for i=1:length(indxs)
        zrm(zrm==indxs(i))=[];
    end
    rnd=randperm(length(zrm));
    rndd=zrm(rnd(1:K));
    sen_indx=[sen_indx,rndd];
end

zrmask=data.zeromask(:,(subi-1)*2+legi);
zrm=find(zrmask==1);
index=[sen_indx,find(zrmask==0)'];
[x,y]=ind2sub([20,37],index);
[xq,yq]=ind2sub([20,37],find(zrmask==1)');

%laplace matrix
sen_mask=zeros(20,37);
sen_mask(index)=1;
lap=zeros(20*37,20*37);
for j=1:37
    for i=1:20
            lap((j-1)*20+i,(j-1)*20+i)=1;
            if sen_mask(i,j)==0
                lap((j-1)*20+i,(j-2)*20+i)=-1/4;
                lap((j-1)*20+i,(j)*20+i)=-1/4;
                lap((j-1)*20+i,(j-1)*20+i-1)=-1/4;
                lap((j-1)*20+i,(j-1)*20+i+1)=-1/4;
            end
    end
end
lap=lap^-1;
%error per distance%
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
for icnt=1:1:nframe
    
    I=data.datatest{subi,legi}(:,icnt);
    I=reshape(I,20,37);
    Iot=zeros(20,37);
    if I==0
        %cnt=cnt+1;
        %continue
    end
    
    if method.name=='Laplace'
    
    b=zeros(20*37,1);
    b(index)=I(index);
    vq=lap*b;
    Iot=reshape(vq,20,37);
    else
    vq = griddata(x,y,I(index),xq,yq,char(method.name));
    Iot(zrm)=vq;
    Iot=abs(Iot);
    end
    
    if fig_flag==1
        fig=figure('Visible','off');
        set(fig,'WindowStyle','docked');
        
        %surf(reshape(I,20,37));
        %saveas(gcf,['Images\img_a' num2str(icnt) '.fig']);
        surf(Iot);
        saveas(gcf,['Images\img_lap' num2str(icnt) '.fig']);
    end
%     if fig_flag==1
%         fig=figure('Visible','off');
%         set(fig,'WindowStyle','docked');
%         map=jet(180);%colormap('jet');
%         subplot(2,1,1);
%         imshow(I,map);
%         subplot(2,1,2);
%         pic=imshow(Iot,map);
%         saveas(pic,['Images\img' num2str(icnt) '.jpg']);
%     end
    
    if kbest.option==1
        RImg=[RImg;I(zrmask~=0),Iot(zrmask~=0)];
    end
    lse=(I-Iot).^2;
    lsey=lsey+lse;
    
    if prog_flag==1
        textprogressbar(icnt/nframe*100);
    end

end

% err=(lsey/(nframe-cnt)).^0.5;
% dist_err=[dist,err];
LSE=(sum(sum(lsey))/(nframe-cnt)/length(find(data.zeromask(:,(subi-1)*2+legi)==1)))^0.5;

if prog_flag==1
    textprogressbar('done');
end