grid;
 set(gcf,'Position' , [0,220 ,640,480])
 set(gcf, 'Renderer', 'ZBuffer') 
xrange=get(gca,'XLim');
yrange=get(gca,'YLim');
yrange(1)=yrange(1)-max(yrange)*0.05;
yrange(2)=yrange(2)+max(yrange)*0.05;
xrange(1)=xrange(1)-max(xrange)*0.05;
xrange(2)=xrange(2)+max(xrange)*0.05;
set(gca,'XLim',xrange,'YLim',yrange);