function linkedPoints(o,thisTag)
% helps set up for linkedPoints
if nargin<2,
    thisTag='linkedPoints';
end
if nargin<1,
    o=gco;
end
set(o,'buttondownfcn',@linkedPoints_callback,'tag',thisTag);

function linkedPoints_callback(o,e)
% allows clicking on points in one figure to highlight the same points in another
% set up with setLinkedPoints()

p=e.IntersectionPoint;
xd=get(o,'xdata');yd=get(o,'ydata');
xdiff=abs(xd-p(1));
ydiff=abs(yd-p(2));
[~,index]=min(xdiff.^2+ydiff.^2);
index=index(1);
fprintf('linked point click: ( %1.3e , %1.3e )\n',xd(index),yd(index));
%index = intersect(find(abs(xd-p(1))<1e-7),find(abs(yd-p(2))<1e-7));
a = get(o,'parent');
f = get(a,'parent');
thisTag = get(o,'tag');
lp = findobj(f,'-regexp','tag',thisTag);
lp = [o; lp(:)];
h=[];
for ii = lp(:)'
    ax=get(ii,'parent');
    xd=get(ii,'xdata');
    yd=get(ii,'ydata');
    c=get(ii,'color');
    thish = get(ii,'userdata');
    if ishandle(thish)
        delete(thish);
    end
    h(end+1)=plot(ax,xd(index),yd(index),'o','markersize',20,'linewidth',2,'color',c);
end
set(o,'userdata',h);


