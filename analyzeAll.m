clc
clear all
close all
warning('off')

load('Dataset\Data');

dict.dir='Dictionaries\';
dict.learn=0; %set to 1 to learn dictionary
dict.iter=25;
dict.S=2;
dict.errorFlag=0;
dict.errorGoal=1;
dict.totalerr=4;
dict.nn=1;
method.S=2;
kbest.sur=[3,3]; %surrounding points about each sensor
kbest.option=0; %set 1 to show pressure vs error

dist_Terr=[];
ReIMG=[];
for subi=1:5
    for legi=1:2
        [img,y((subi-1)*2+legi,:),dist_terr,ReImg]=testallfunc(data,subi,legi,dict,method,kbest);
        num((subi-1)*2+legi,1:size(y,2))=length(find(data.zeromask(:,(subi-1)*2+legi)==1));
        dist_Terr=[dist_Terr;dist_terr];
        ReIMG=[ReIMG;ReImg];
    end
    subi
end

LSET=(sum(y.^2.*num)./sum(num)).^0.5;
x=1:length(LSET);

img.title=['Reconstruction Using OMP for All Subjects'];

my_plot_all(x,LSET,img)


%%%%%%%%%distance plot%%%%%%%%
if kbest.option==1
kvalue=[0 ];%1 5 7 10 14 20];
kbest.K=1;
nbin=20;
for iter=1:size(dist_Terr,2)
    x=zeros(1,nbin);
    y=zeros(1,nbin);
    temp=vertcat(dist_Terr{:,iter});
    [dist,index]=sort(temp(:,1));
    err=temp(index,2);
    zrmax=find(dist==0,1,'last');
    x(1)=0;y(1)=mean(err(1:zrmax));
    nstep=round((length(err)-zrmax)/(nbin-1));
    for i=2:nbin
        if i~=nbin
            x(i)=mean(dist(zrmax+(i-2)*nstep+1:zrmax+(i-1)*nstep));
            y(i)=mean(err(zrmax+(i-2)*nstep+1:zrmax+(i-1)*nstep));
        else
            x(i)=mean(dist(zrmax+(i-2)*nstep+1:end));
            y(i)=mean(err(zrmax+(i-2)*nstep+1:end));
        end
    end

    err=std(y,[],1);
    if iter==1
        img.title=['Reconstruction Error Vs Distance From Sensor For All Subjects'];
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
    temp=vertcat(ReIMG{:,iter});
    [pres,index]=sort(temp(:,1));
    err=temp(index,2)-pres;
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
