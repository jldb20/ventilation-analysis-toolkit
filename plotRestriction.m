function plotRestriction(s_corr,c_disp,maxme,s_corr_markers,s_corr_labs,cycles)
% plots the parameters related to the restriction between the two meters


% controls whether windows fill screen or not
if nargin<3,
    maxme=false;
end

% ugh. this code needs a redesign.
if nargin<4,
    s_corr_markers(2,1:2) = [numel(s_corr(1).t) numel(s_corr(2).t)];
    s_corr_labs{1} = 'single set';
end
if nargin<5,
    for ii=1:size(s_corr_markers,1),
        s_corr_labs{ii} = 'unlabelled';
    end
end
if nargin<6,
    cycles=[];
end
    
figure(301);clf;
if maxme && ~strcmpi(get(gcf,'units'),'normalized'), set(gcf,'units','normalized','position',[0 0 1 1]); end

% plot styles for multiple datasets
ptStyle={'.','s','x','^','o','<','>','*'};       % different styles

% now plot p and Q signals vs time and on same graph
subplot(2,2,1);hold on;
timePlot(s_corr,s_corr_markers,s_corr_labs,'p',cycles);

subplot(2,2,2);hold on;
timePlot(s_corr,s_corr_markers,s_corr_labs,'Q',cycles,true);

% plot datapoints for multiple datasets:
subplot(2,1,2);hold on; grid on
set(gca,'tag','RestrictorFig');
ptCol = [0    0.4470    0.7410]; % all the same color 
leg = {};
for ii=1:size(s_corr_markers,1)-1,
    s_corr_range = s_corr_markers(ii)+1:s_corr_markers(ii+1);
    h=plot(s_corr(3).Q(s_corr_range),s_corr(3).p(s_corr_range),ptStyle{ii},'color',ptCol);
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    leg{end+1}=sprintf('measured (%s)',s_corr_labs{ii});
end
xlabel('Q /Ls^{-1}'); ylabel('\Delta p /cmH2O'); grid on;

% get title from labels, or if the labels data file can be found then
% import that and use those titles.
tit = s_corr_labs{1};
for ii=2:numel(s_corr_labs),
    tit=[tit ', ' s_corr_labs{ii}];
end
tit = [tit ': ' s_corr(1).label];
title(tit);

% reconstruct best fit
s_fit.Q = linspace(0,max(s_corr(3).Q),30).';
s_fit.p = [s_fit.Q.*abs(s_fit.Q) s_fit.Q ones(size(s_fit.Q))] * c_disp(1,:)';
plot(s_fit.Q, s_fit.p, 'g-','linewidth',2);
leg{end+1}=sprintf('fit   c2: %0.1e   c1: %0.1e   c0: %0.1e',c_disp(1,:));
for ii=2:size(c_disp,1),
    s_fit.p = [s_fit.Q.*abs(s_fit.Q) s_fit.Q ones(size(s_fit.Q))] * c_disp(ii,:)';
    plot(s_fit.Q, s_fit.p,'linewidth',1);
    leg{end+1}=sprintf('fit   c2: %0.1e   c1: %0.1e   c0: %0.1e',c_disp(ii,:));
end
legend(leg,'location','northwest');

% set better limits
ylim(1) = min(0,min(s_corr(3).p(s_corr(3).Q~=0)));
ylim(2) = max(0,max(s_corr(3).p(s_corr(3).Q~=0)));
set(gca,'ylim',ylim);

sgtitle('restriction between meters');




function timePlot(s_corr,s_corr_markers,s_corr_labs,plotvar,cycles,plotCombined)
% plots flow measurements against time
% s_corr is the signal structure with fields t,p,V,Q
% s_corr_markers is the vector of start and end indices if s_corr contains
%                multiple datasets
% s_corr_labs is the dataset labels if s_corr contains multiple datasets
% plotvar is a string describing which field to plot (i.e. 'p', 'V', or 'Q')
% cycles is a cell array of the variable cycles that is returned from
%        splitCycle, as modified by identifyRestriction2 for plotting here.

% option for plotting difference or mean strored in s_corr(3)
if nargin<6 || isempty(plotCombined),
    plotCombined = false;
end

if nargin>4 && ~isempty(cycles), % plots the underlying data for the phase-averaged cases
    for ii=1:numel(cycles),
        for jj=1:numel(cycles{ii}{1})
            hx = plot(cycles{ii}{1}(jj).tau,cycles{ii}{1}(jj).(plotvar),'x','color','#bbbbff');
        end
        for jj=1:numel(cycles{ii}{2})
            o = plot(cycles{ii}{2}(jj).tau,cycles{ii}{2}(jj).(plotvar),'x','color','#ffbbbb');
        end
        hx=[hx;o]; % for legend below
    end
end
for ii=1:size(s_corr_markers,1)-1, % plots the phase-averaged cases or the full dataset if not phase-averaged
    s_corr_range = s_corr_markers(ii)+1:s_corr_markers(ii+1);
    h=plot(s_corr(1).t(s_corr_range),s_corr(1).(plotvar)(s_corr_range),'b-','linewidth',2);
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    h(2)=plot(s_corr(2).t(s_corr_range),s_corr(2).(plotvar)(s_corr_range),'r-','linewidth',2);
    linkedPoints(h(2),sprintf('linkedPts%1.0d',ii));
    if size(s_corr_markers,1)>2, % if there's more than one dataset, label datasets
        text(s_corr(2).t(s_corr_range(end)),s_corr(2).(plotvar)(s_corr_range(end)),s_corr_labs{ii});
    end
    if plotCombined,
        h(3)=plot(s_corr(2).t(s_corr_range),s_corr(3).(plotvar)(s_corr_range),'g--','linewidth',2); 
        linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    end
end
leg={'meter 1 (phase av.)','meter 2 (phase av.)'};
if plotCombined,
    switch plotvar,
        case 'p',
            leg{end+1}=('difference');
        case 'Q',
            leg{end+1}=('mean');
        otherwise
            warning('You''re plotting a field that was not anticipated when writing this function. Be sure you know what you''re doing.');
    end
end

% for displaying number of cycles
meter1_cycles = '';
meter2_cycles = '';
for ii=1:numel(cycles),
    meter1_cycles = sprintf('%s%2.0d,',meter1_cycles,numel(cycles{ii}{1}));
    meter2_cycles = sprintf('%s%2.0d,',meter2_cycles,numel(cycles{ii}{2}));
end
meter1_cycles(end)=[];
meter2_cycles(end)=[];
    
leg={leg{:},sprintf('meter 1 (%s cycles)',meter1_cycles),sprintf('meter 2 (%s cycles)',meter2_cycles)};
xlabel('time /s'); legend([h(:);hx(:)],leg); grid on;
switch plotvar,
    case 'p',
        ylabel('p /cmH2O');
    case 'Q',
        ylabel('Q /(L/s)');
    case 'V',
        ylabel('V /L');
    otherwise
        warning('You''re plotting a field that was not anticipated when writing this function. Be sure you know what you''re doing.');
end

%        figure(305);clf; hold on;
matchXaxes;

