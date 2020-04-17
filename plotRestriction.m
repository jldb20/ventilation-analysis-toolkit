function plotRestriction(s_corr,c_disp,maxme,s_corr_markers,s_corr_labs)
% plots the parameters related to the restriction between the two meters


% controls whether windows fill screen or not
if nargin<3,
    maxme=false;
end

if nargin<6,
    saveMyFig=false;
end

% ugh. this code needs a redesign.
if nargin<4,
    s_corr_markers(2,1:2) = [numel(s_corr(1).t) numel(s_corr(2).t)];
    s_corr_labs{1} = 'single set';
end
if nargin<5,
    for ii=1:size(s_corr_markers,1),
        s_corr_labs{1} = 'unlabelled';
    end
end
    
figure(301);clf;
if maxme && ~strcmpi(get(gcf,'units'),'normalized'), set(gcf,'units','normalized','position',[0 0 1 1]); end

% now plot p and Q signals on same graphs
subplot(2,2,1);hold on;
for ii=1:size(s_corr_markers,1)-1,
    s_corr_range = s_corr_markers(ii)+1:s_corr_markers(ii+1);
    h=plot(s_corr(1).t(s_corr_range),s_corr(1).p(s_corr_range),'b-');
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    h=plot(s_corr(2).t(s_corr_range),s_corr(2).p(s_corr_range),'r-');
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
end
xlabel('time /s');ylabel('p /cmH2O'); legend('meter 1','meter 2'); grid on;
%        figure(305);clf; hold on;
matchXaxes;
subplot(2,2,2);hold on;
for ii=1:size(s_corr_markers,1)-1,
    s_corr_range = s_corr_markers(ii)+1:s_corr_markers(ii+1);
    h=plot(s_corr(1).t(s_corr_range),s_corr(1).Q(s_corr_range),'b-');
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    h=plot(s_corr(2).t(s_corr_range),s_corr(2).Q(s_corr_range),'r-'); %legend('meter 1','meter 2'); grid on;
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    h=plot(s_corr(2).t(s_corr_range),s_corr(3).Q(s_corr_range),'g--'); 
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
end
legend('meter 1','meter 2','mean'); grid on;
xlabel('time /s');ylabel('Q /Ls^{-1}');
matchXaxes; % so we can zoom in and compare p & Q data easily

% plot datapoints for multiple datasets:
subplot(2,1,2);hold on; grid on
set(gca,'tag','RestrictorFig');
ptStyle={'.','s','x','^','o','<','>','*'};       % different styles
ptCol = [0    0.4470    0.7410]; % all the same color 
leg = {};
for ii=1:size(s_corr_markers,1)-1,
    s_corr_range = s_corr_markers(ii)+1:s_corr_markers(ii+1);
    h=plot(s_corr(3).Q(s_corr_range),s_corr(3).p(s_corr_range),ptStyle{ii},'color',ptCol);
    linkedPoints(h,sprintf('linkedPts%1.0d',ii));
    leg{end+1}=sprintf('measured (%s)',s_corr_labs{ii});
end
xlabel('Q /Ls^{-1}'); ylabel('\Delta p /cmH2O'); grid on;
tit = s_corr_labs{1};
for ii=2:numel(s_corr_labs),
    tit=[tit ', ' s_corr_labs{ii}];
end
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
