function s_corr = doCorrelate(s,cycle_length,debug)
% Returns correlated (aligned) versions of the signals in s
% s is a signal as returned by importWrapper.
% cycle_length is in samples and is used to determine how far to try
% shifting the two signals in order to find an alignment.
% debug is boolean to decide if debugging output is wanted

if nargin<3,
    debug=false;
end

% correlate using Q:
z1=s(1).Q; % variable to correlate in alignment
z2=s(2).Q; % variable to correlate in alignment
corr_length = min(length(z1),length(z2))-cycle_length;
for ii=1:cycle_length,
    corr(ii) = sum(z1(1:corr_length) .* z2(ii+(0:corr_length-1)));
    if debug,
        % debugging:
        Q2o = s(2).Q(ii+(0:corr_length-1));
        t2o = s(2).t(1:corr_length);
        figure(399); clf;
        subplot(2,1,2);plot(corr);
        title('press any key (repeatedly)');
        subplot(2,2,1); hold on;
        plot(s(1).t,s(1).Q,'b-'); plot(t2o,Q2o,'r-');
        setMatchXaxes;
        ax=subplot(2,2,2); hold on;
        plot(t2o,z1(1:corr_length),'b-'); plot(t2o,z2(ii+(0:corr_length-1)),'r-');
        setMatchXaxes;
        if exist('xl','var'), set(ax,'xlim',xl); end
        pause
        xl=get(ax,'xlim');
    end
end
[~,offset] = max(corr);
s_corr(2).Q = s(2).Q(offset(1)+(0:corr_length-1));
s_corr(2).p = s(2).p(offset(1)+(0:corr_length-1));
s_corr(2).t = s(2).t(1:corr_length);
s_corr(1).t = s(1).t(1:corr_length);
s_corr(1).p = s(1).p(1:corr_length);
s_corr(1).Q = s(1).Q(1:corr_length);

if debug,
    close(399);
end
