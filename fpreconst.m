function [LSE,mu_coh,dist_err,RImg]=fpreconst(data,subi,legi,dict,method,kbest,K,varargin)

nVararg=length(varargin);
%varargin=[nframe,prog_flag,zr_flag,fig_flag,blkSize,slidingDis]

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
% ict.disp=0;

% method.name='OMP';
% method.S=2; %level of sparsity for OMP algorithm
% method.p=1;

% K=5; %number of random sensors;
% kbest.K=1; %number of sensors in each pre-defined region
% kbest.mode %'peaks' & 'rnd'
% kbest.sur %surrounding points to consider

if nVararg<6
    slidingDis=1; %sliding distance
else
    slidingDis=varargin{6};
end
if nVararg<5
    blkSize=[20,37]; %image size or block size
else
    blkSize=varargin{5};
end
if nVararg<4
    fig_flag=0;  %save real vs reconstructed image
else
    fig_flag=varargin{4};
end
if nVararg<3
    zr_flag=0; %if 1, include zero values in reconstruction phase
else
    zr_flag=varargin{3};
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

%fig_flag=1;
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
if zr_flag==1
    rnds=find(data.zeromask(:,(subi-1)*2+legi)==0)';
    sen_indx=[sen_indx,rnds];
end

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

%error per distance%
zrmask=data.zeromask(:,(subi-1)*2+legi);
[xd,yd]=ind2sub([20,37],(1:20*37)');
[xs,ys]=ind2sub([20,37],sen_indx);
dist=repmat([xd(zrmask~=0),yd(zrmask~=0)],1,1,kbest.K*4+K)-...
    repmat(reshape([xs;ys],1,2,kbest.K*4+K),length(find(zrmask~=0)),1,1);
dist=(min(sum(dist.^2,2),[],3)).^0.5;
%%%%%%%%%%%%%%%

I_temp=zeros(20*37,1);
I_temp(sen_indx,1)=1;
I_temp=reshape(I_temp,20,37);
blocks2=my_im2col(I_temp,blkSize,slidingDis);

for i=1:size(blocks2,2)
    [ix,iy]=find(reshape(blocks2(:,i),blkSize)==1);
    len_ix(i)=length(ix);
    for j=1:length(ix)
        ixt=ix(j)-ceil((kbest.sur(1)-1)/2):ix(j)-ceil((kbest.sur(1)-1)/2)+kbest.sur(1)-1;
        if min(ixt)<1
            ixt=ixt-min(ixt)+1;
        elseif max(ixt)>blkSize(1)
            ixt=ixt-(max(ixt)-blkSize(1));
        end
        iyt=iy(j)-ceil((kbest.sur(2)-1)/2):iy(j)-ceil((kbest.sur(2)-1)/2)+kbest.sur(2)-1;
        if min(iyt)<1
            iyt=iyt-min(iyt)+1;
        elseif max(iyt)>blkSize(2)
            iyt=iyt-(max(iyt)-blkSize(2));
        end
        ixy=(iyt-1)*blkSize(1)+ixt;
        DCTmt=zeros(1,blkSize(1)*blkSize(2));
        DCTmt(ixy)=1/(kbest.sur(1)*kbest.sur(2));
        DCTm(j,:,i)=DCTmt;
    end
end
dict.DCTm=DCTm;

%dictionary formation

if dict.type=='learned'
    %load or learn a dictionary
    if dict.learn==0
        load(dict.name)
    elseif dict.learn==1
        DCT=learndic(data,subi,legi,dict);
    end
elseif dict.type=='DCT'
    %DCT dictionary
    DCT1=dct(eye(blkSize(1)));
    DCT2=dct(eye(blkSize(2)));
    DCT=kron(DCT1,DCT2)';
elseif dict.type=='DWT'
    %DWT dictionary
    Iden=eye(blkSize(1));
    DWT1=zeros(blkSize(1),blkSize(1));
    for i=1:blkSize(1)
        DWT1(:,i)=wavedec(Iden(:,i),nextpow2(blkSize(1)),'haar');
    end

    Iden=eye(blkSize(2));
    DWT2=zeros(blkSize(2),blkSize(2));
    for i=1:blkSize(2)
        DWT2(:,i)=wavedec(Iden(:,i),nextpow2(blkSize(2)),'haar');
    end
    
    DCT=kron(DWT1,DWT2)';
end

%reconstruction phase

if prog_flag==1
    textprogressbar('Reconstructing: ');
end


%rndm=(1+randn(size(DCTm(1:len_ix(1),:,1),1),1)*0.05); %randomness
cnt=0;
lsey=0;
RImg=[];
for icnt=1:1:nframe
        
    I=data.datatest{subi,legi}(:,icnt);
    I=reshape(I,20,37);
    [m,n]=size(I);

    blocks1=my_im2col(I,blkSize,slidingDis);

    i1=ceil((m-blkSize(1))/slidingDis)+1;
    j1=ceil((n-blkSize(2))/slidingDis)+1;

    Iot=zeros(m,n);
    Icnt=zeros(m,n);
    for i=1:size(blocks1,2)
        block=DCTm(1:len_ix(i),:,i)*blocks1;
        %block=block.*rndm; %randomness
        if block==0
            Io=zeros(blkSize(1),blkSize(2));
            %cnt=cnt+1;
        else
            if method.name=='FOCUSS'
                Io=RFOCUSS(DCTm(1:len_ix(i),:,i)*DCT,block,method.p,1,1e-3); %change p
            elseif method.name=='OMP'
                Io=ORMP(DCTm(1:len_ix(i),:,i)*DCT,block,method.S);
            elseif method.name=='LASSO'
                [Io,fit]=lasso(DCTm(1:len_ix(i),:,i)*DCT,block,'NumLambda',20,'Alpha',method.p);
                [~,i_mse]=min(fit.MSE);
                Io=Io(:,i_mse);
            end
     
            Io=DCT*Io;
            Io=reshape(Io,blkSize);
        end
    cl=ceil(i/i1);
    rw=i-(cl-1)*i1;
    if rw~=i1
        indx=(rw-1)*slidingDis+1:(rw-1)*slidingDis+blkSize(1);
    else
        indx=m-blkSize(1)+1:m;
    end
    if cl~=j1
        indy=(cl-1)*slidingDis+1:(cl-1)*slidingDis+blkSize(2);
    else
        indy=n-blkSize(2)+1:n;
    end

    Iot(indx,indy)=Iot(indx,indy)+Io;
    Icnt(indx,indy)=Icnt(indx,indy)+1;

    end
    Iot=Iot./Icnt;
    Iot(data.zeromask(:,(subi-1)*2+legi)==0)=0;
    %Iot(Iot<0)=0;
    Iot=abs(Iot);

    
    if fig_flag==1
        fig=figure('Visible','off');
        set(fig,'WindowStyle','docked');
        
        %surf(reshape(I,20,37));
        %saveas(gcf,['Images\img_a' num2str(icnt) '.fig']);
        surf(Iot);
        saveas(gcf,['Images\img_sp4' num2str(icnt) '.fig']);
    end
    
%     if fig_flag==1
%         %fig=figure('Visible','off');
%         %set(fig,'WindowStyle','docked');
%         map=jet(180);%colormap('jet');
%         subplot(2,1,1);
%         imshow(I,map);
%         subplot(2,1,2);
%         pic=imshow(Iot,map);
%         saveas(pic,['Images\img' num2str(icnt) '.fig']);
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

mu_coh=corr(DCTm(1:len_ix(1),:,1)*DCT);
mu_coh(isnan(mu_coh))=1;
mu_coh=mu_coh(:);
mu_coh(1:blkSize(1)*blkSize(2)+1:end)=[];
maxcoh=mean(abs(mu_coh));

err=(lsey(zrmask~=0)/(nframe-cnt)).^0.5;
dist_err=[dist,err];
LSE=(sum(sum(lsey))/(nframe-cnt)/length(find(data.zeromask(:,(subi-1)*2+legi)==1)))^0.5;

if prog_flag==1
    textprogressbar('done');
end