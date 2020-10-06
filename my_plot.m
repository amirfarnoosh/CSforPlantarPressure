function my_plot(x,y,err,img)

% img.title
% img.xlabel
% img.ylabel
% img.xlim
% img.xtick
% img.xticklabel

if img.hold==0
     fig=figure;
     fig.Color='white';
     fig.WindowStyle='docked';
else
    hold on
end

plt=errorbar(x,y,err);
grid on

if img.hold==1
    legend(img.legend);
end

%plot properties
plt.LineStyle='-';
plt.LineWidth=1;
plt.Marker='hexagram';
plt.MarkerSize=3;

%title properties  
tlt=title(img.title);
tlt.LineWidth=3;
tlt.FontName='Times New Roman';
tlt.FontSize=10;

%x & y label properties
xl=xlabel(img.xlabel);
xl.LineWidth=3;
xl.FontName='Times New Roman';
xl.FontSize=8;

yl=ylabel(img.ylabel);
yl.LineWidth=3;
yl.FontName='Times New Roman';
yl.FontSize=8;

ax = gca;
ax.FontSize =6;

ax.XLim =img.xlim;
ax.XTick=img.xtick;
ax.XTickLabel=img.xticklabel;

%ax.YLim =img_ylim;
%ax.YTick=img_ytick;
%ax.YTickLabel=img_yticklabel;