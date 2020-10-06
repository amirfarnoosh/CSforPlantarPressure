function my_plot_all(x,y,img)

% img.title
% img.xlabel
% img.ylabel
% img.xlim
% img.xtick
% img.xticklabel

fig=figure;
fig.Color='white';
fig.WindowStyle='docked';

plt=plot(x,y);
grid on

%plot properties
plt.LineStyle='-';
plt.LineWidth=3;
plt.Marker='hexagram';
plt.MarkerSize=6;

%title properties  
tlt=title(img.title);
tlt.LineWidth=3;
tlt.FontName='Times New Roman';
tlt.FontSize=18;

%x & y label properties
xl=xlabel(img.xlabel);
xl.LineWidth=3;
xl.FontName='Times New Roman';
xl.FontSize=14;

yl=ylabel(img.ylabel);
yl.LineWidth=3;
yl.FontName='Times New Roman';
yl.FontSize=14;

ax = gca;
ax.FontSize =18;

ax.XLim =img.xlim;
ax.XTick=img.xtick;
ax.XTickLabel=img.xticklabel;

%ax.YLim =img_ylim;
%ax.YTick=img_ytick;
%ax.YTickLabel=img_yticklabel;