function matchXaxes(h)
if nargin<1,
    h=gca;
end
% helps set up for matchXaxes
addlistener(h, 'XLim', 'PostSet', @matchXaxes_callback); % so we can zoom in and compare p & Q data easily
set(h, 'userdata', 'matchXaxes');


function matchXaxes_callback(o,e)
% matches axes in the same figure; this should be set as a listener
% and then set the userdata of all axes to match to a string including
% 'matchXaxes': e.g.
% addlistener(gca, 'XLim', 'PostSet', @matchXaxes);
% set(gca, 'userdata', 'matchXaxes');
% OR use setMatchXaxes() ...
a=e.AffectedObject;
fig = get(a,'parent');
xl = get(a,'xlim');
cc=get(fig,'children');
for jj=1:numel(cc),
    if strcmp(get(cc(jj),'type'),'axes') && strcmp(get(cc(jj),'userdata'),'matchXaxes'),
        set(cc(jj),'xlim',xl);
    end
end