function [img,y,dist_terr,ReImg,W]=testallfunc(data,subi,legi,dict,method,kbest)

W=0;
% subi=3; %subject number
% legi=1; %leg number 1 for left and 2 for right

dict.type=string('learned'); % 'learned', 'DCT' and 'DWT'
%dict.learn=0; %set to 1 to learn dictionary
% dict.iter=5;
% dict.S=2;
%dict.name='DCTt_2_5nn_sub21';
%dict.dir='C:\Users\Amir\Documents\MATLAB\dics\';
if dict.nn==1
    dict.name=[dict.dir,'DCTt_',num2str(dict.S),'_',num2str(dict.iter),'nn_sub',num2str(subi),num2str(legi)]; %if type is 'learned' load this dictionary
elseif dict.nn==0
    dict.name=[dict.dir,'DCTt_',num2str(dict.S),'_',num2str(dict.iter),'_sub',num2str(subi),num2str(legi)];
end

dict.blkSize=[20,37];
dict.slidingDis=1;
dict.disp=1;

method.name=string('OMP');
%method.S=2; %level of sparsity for OMP algorithm
method.p=1;

% K=5; %number of random sensors;
kbest.K=1; %number of sensors in each pre-defined region
kbest.mode=string('peaks'); %'peaks' & 'rnd'
%kbest.sur=[1,1]; %surrounding points about each sensor

kvalue=[0 1 5 7 10 14 20 29 42];% 60 85 122 200];
cnt=0;
dist_terr=cell(1,length(kvalue));
ReImg=cell(1,length(kvalue));
for K=kvalue
    cnt=cnt+1;
    for iter=1:1
        
        [LSE(iter,cnt),mu_coh(:,cnt),dist_err,RImg]=fpreconst(data,subi,legi,dict,method,kbest,K,90,1,[],1);  %change to this for sparse methods
        %also uncomment end of code for mutual coherence
        %[LSE(iter,cnt),dist_err,RImg]=fpreconst_interp(data,subi,legi,method,kbest,K,90,1);  %change to this for interpolation methods
        %[LSE(iter,cnt),dist_err,RImg,W]=fpreconst_gmm(data,subi,legi,method,kbest,K,90,1);
        dict.learn=0;
        dist_terr{cnt}=[dist_terr{cnt};dist_err];
        ReImg{cnt}=[ReImg{cnt};RImg];
    end
end

x=1:length(kvalue);
y=mean(LSE,1);
err=std(LSE,[],1);

img.hold=0;
img.title=['Reconstruction Using OMP for SUB#',num2str(subi),num2str(legi)];
img.xlabel='Number of Sensors';
img.ylabel='RMS Error';
img.xlim=[1 length(kvalue)+1];
img.xtick=1:length(kvalue)+1;

for i=1:length(kvalue)
    img.xticklabel{1,i}=4*kbest.K+kvalue(i);
end

my_plot(x,y,err,img)

%%%%%%%%%distance plot%%%%%%%%
if kbest.option==1
nbin=20;
for iter=1:length(kvalue)
    x=zeros(1,nbin);
    y=zeros(1,nbin);
    [pres,index]=sort(dist_terr{1,iter}(:,1));
    err=dist_terr{1,iter}(index,2);
    zrmax=find(pres==0,1,'last');
    x(1)=0;y(1)=mean(err(1:zrmax));
    nstep=round((length(err)-zrmax)/(nbin-1));
    for i=2:nbin
        if i~=nbin
            x(i)=mean(pres(zrmax+(i-2)*nstep+1:zrmax+(i-1)*nstep));
            y(i)=mean(err(zrmax+(i-2)*nstep+1:zrmax+(i-1)*nstep));
        else
            x(i)=mean(pres(zrmax+(i-2)*nstep+1:end));
            y(i)=mean(err(zrmax+(i-2)*nstep+1:end));
        end
    end

    err=std(y,[],1);
    if iter==1
        img.title=['Reconstruction Error Vs Distance From Sensor#',num2str(subi),num2str(legi)];
        img.xlabel='Distance from Nearest Sensor';
        img.ylabel='RMS Error';
        img.xlim=[0 max(x)];
        img.xtick=0:1:max(x);
    end
    for i=1:length(img.xtick)
        img.xticklabel{1,i}=img.xtick(i);
    end
    
    img.legend{1,iter}=['K = ' num2str(4*kbest.K+kvalue(iter))];
    my_plot(x,y,err,img)
    img.hold=1;
end
img.hold=0;
end
%%%%%%%%%%%


%%%%%%%%%pressure vs error plot%%%%%%%%
if kbest.option==1
nbin=50;
for iter=1:length(kvalue)
    x=zeros(1,nbin);
    y=zeros(1,nbin);
    x_over=zeros(1,nbin);
    y_over=zeros(1,nbin);
    x_under=zeros(1,nbin);
    y_under=zeros(1,nbin);
    [pres,index]=sort(ReImg{1,iter}(:,1));
    err=ReImg{1,iter}(index,2)-pres;
    err_over=err(err>0);
    err_under=err(err<0);
    pres_over=pres(err>0);
    pres_under=pres(err<0);
    zrmax=find(pres==0,1,'last');
    zrmax_over=find(pres_over==0,1,'last');
    zrmax_under=0;
    x(1)=0;y(1)=(mean(err(1:zrmax).^2))^0.5;
    x_over(1)=0;y_over(1)=(mean(err_over(1:zrmax_over).^2))^0.5;
    x_under(1)=0;y_under(1)=0;
    nstep=round((length(err)-zrmax)/(nbin-1));
    nstep_over=round((length(err_over)-zrmax_over)/(nbin-1));
    nstep_under=round((length(err_under)-zrmax_under)/(nbin-1));
    for i=2:nbin
        if i~=nbin
            x(i)=mean(pres(zrmax+(i-2)*nstep+1:zrmax+(i-1)*nstep));
            y(i)=(mean(err(zrmax+(i-2)*nstep+1:zrmax+(i-1)*nstep).^2))^0.5;
            x_over(i)=mean(pres_over(zrmax_over+(i-2)*nstep_over+1:zrmax_over+(i-1)*nstep_over));
            y_over(i)=(mean(err_over(zrmax_over+(i-2)*nstep_over+1:zrmax_over+(i-1)*nstep_over).^2))^0.5;
            x_under(i)=mean(pres_under(zrmax_under+(i-2)*nstep_under+1:zrmax_under+(i-1)*nstep_under));
            y_under(i)=(mean(err_under(zrmax_under+(i-2)*nstep_under+1:zrmax_under+(i-1)*nstep_under).^2))^0.5;
        else
            x(i)=mean(pres(zrmax+(i-2)*nstep+1:end));
            y(i)=(mean(err(zrmax+(i-2)*nstep+1:end).^2))^0.5;
            x_over(i)=mean(pres_over(zrmax_over+(i-2)*nstep_over+1:end));
            y_over(i)=(mean(err_over(zrmax_over+(i-2)*nstep_over+1:end).^2))^0.5;
            x_under(i)=mean(pres_under(zrmax_under+(i-2)*nstep_under+1:end));
            y_under(i)=(mean(err_under(zrmax_under+(i-2)*nstep_under+1:end).^2))^0.5;
        end
    end

    err=std(y,[],1);
    err_over=std(y_over,[],1);
    err_under=std(y_under,[],1);
    
    img.title=['Reconstruction Error Vs Real Sensor Pressure #',num2str(subi),num2str(legi),' K = ' num2str(4*kbest.K+kvalue(iter))];
    img.xlabel='Real Pressure';
    img.ylabel='RMS Error';
    img.xlim=[0 max(x)];
    img.xtick=0:25:max(x);
    for i=1:length(img.xtick)
        img.xticklabel{1,i}=img.xtick(i);
    end
    
    img.legend{1,1}='both';
    my_plot(x,y,err,img)
    img.legend{1,2}='over-estimate';
    img.hold=1;
    my_plot(x_over,y_over,err_over,img)
    img.legend{1,3}='under-estimate';
    my_plot(x_under,y_under,err_under,img)
    img.hold=0;
end
img.hold=0;
end
%%%%%%%%%%%


%mutual coherence
% % img.title='Mutual Coherence';
% % img.ylabel='(1+mue^-1)/2';
% % 
% % maxcoh=(1+maxcoh.^-1)/2;
% % mue_coh=mean(maxcoh,1);
% % err=std(mue_coh,[],1);
% % my_plot(x,mue_coh,err,img)