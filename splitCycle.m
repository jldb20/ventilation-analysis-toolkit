function [partitioned,stats,cycles]=splitCycle(s,threshold,margin)
%
% This is a bit rough and ready for now, and I can't help thinking there
% must be a simpler way to go about it. But it works, so it may never
% change...
%
% In simplest usage: 
%
% splitCycle(sig) produces some nice plots
%
% partitioned_sig = splitCycle(sig)
% sig is an input file from importWrapper() : an n-element structure where
%     n is the number of meters and the fields in the structure are p, V, Q
%     and t.
% partitioned_sig also contains p, Q, V, and t, but these are now divided
%                 into segments by NaN values, with each segment representing
%                 a complete inspiration or expiration part of a cycle.
%                 There is also a "tau" field, which starts the time at
%                 zero for each inspiration-expiration cycle (e.g. try
%                 plot(partitioned_sig.tau,partitioned_sig.p);).
%                 Finally, there are also _pos and _neg suffixes which
%                 separate the datasets into inspiration and expiration
%                 segments, respectively.
%
% [partitioned,stats,cycles]=splitCycle(s,threshold,margin)
% s as before
% threshold is the low flow cutout to use (flow below this threshold set to
%           zero)
% margin is the margin to remove at beginning and end of each inspiration
%        and expiration phase. It is a scalar or a two element vector
%        [start_margin end_margin] and defines the number of seconds to be
%        cropped. This is to avoid transient dynamic effects in the
%        identification of restriction parameters.
% stats contains: sample_rate, RR (respiration rate), Vin_frac (to identify
%                 volume going in as total volume passing through in either
%                 direction), and inspiration/expiration mode.
% cycles is a cell array, where each cell is a structure array for one of
%        the meters. The structure array has one element for each cycle,
%        and has fields containing, amongst other things, the timeseries
%        data for the positive (inspiration) and negative (expiration)
%        phases for that cycle.
%
%
%
% Some more somewhat cryptic documentation below. To be updated. Hopefully.
% =========================================================================
%
% splits sig data (pressure, flow, volume, time) into individual cycles
% inputs:
% -s is the signal from importWrapper
% -threshold is the low flow cutout to use (if it's been done in the data
% already, which it usually has, then feel free to omit this or give [] and
% the relatively low default of 0.01 L/s will be used.
% -margin is a margin to apply on either side of each positive region and
% each negative region: this amount of time will be shaved off both ends of
% the dataset. If margin is a 2-element vector then start and end are
% specified independently. NB the margin is applied to the indices vectors,
% covering all the points in the regions, but not to the individual start
% and end indices recorded for each region. From the vectors, it affects
% the datasets returned for p, Q, V and t. (see the plot it makes to get an
% idea what it's doing if this is unclear) Also NB, it does this for the
% pos and neg indices but not for the full cycle indices because this is
% not how the margin is meant to work (it's meant to be applied to the pos
% and neg regions separately).
% -plt is a boolean indicating whether or not to plot a figure (default no)
% outputs:
% -cycles is the data sorted into a structure array (kinda cluttered but
% lots of useful versions of the data and corresponding indices in there.
% Could do with documenting but that's not going to happen any time soon). 
% -cropped is a single set of p, Q, V, t and tau data with regions of
% unwanted data removed and replaced with a single nan value. (The NaN is
% so it can be plotted easily as individual disconnected segments with a
% single plot command. It can be removed with e.g.
% cropped.p(isnan(cropped.p))=[] ). This is probably the most useful
% output. That's why I've moved it first. By the way, tau is the cycle
% time, starting from zero for each cycle.
%
% this function is way more complicated than it needs to be. Needs to be
% simplified.
% I'd like to make it so it accepts only a single s input (not an array of
% structures). This would make it marginally less complicated and
% error-prone. Figure plotting would need to move outside this function,
% but that's not a bad thing.

% input validation & defualts
if nargin<2 || isempty(threshold), % threshold below which flow is considered zero
    threshold = 0.01;
end
if nargout==0, % plot figure(s) if no outputs are requested
    plt=true;
else           % otherwise, don't plot
    plt=false;
end
if nargin<3 || isempty(margin),
    margin = 0; % seconds
end
if numel(margin)<2,
    margin = margin * [1 1]; % start and end margins
end

%  number of meters saved in the dataset
n_meters = numel(s);

% so the figures for old data get closed and don't confuse things
fig1=findobj('-regexp','tag','ventilator-splitCycle1');
fig2=findobj('-regexp','tag','ventilator-splitCycle2');
if ~plt,
    close(fig1);
    close(fig2);
else
    close(fig1(2:end)); % just in case
    close(fig2(2:end)); % just in case
    if n_meters<2
        close(fig2);
    end
end

% now loop for all meters in dataset
for n = 1:n_meters,
    fprintf('METER %0.0d\n',n);
    
    % for brevity, and to be sure of vector shape
    p=s(n).p(:);
    Q=s(n).Q(:);
    t=s(n).t(:);
    V=s(n).V(:);
    
    % find the sample rate
    dt = diff(t);
    sample_time = mean(dt);
    if any(diff(dt)./mean(dt) > 1e-10),
        % This checks that the sampling rate is practically constant
        % (allowing for numerical rounding).
        error('Time series must be sampled at regular intervals.');
    end
    stats.sample_rate = 1/sample_time;
    
    % find respiration rate (to nearest sample time) using autocorrelation
    a=zeros(numel(Q),1);
    for ii=1:numel(Q),
        a(ii)=Q(ii:end)'*Q(1:end-ii+1);
    end
    % set first 0.5 s of correlation data to zero to avoid identifying
    % a period of zero. (Fine up to around RR = 2/s = 120/min)
    a(1:ceil(stats.sample_rate/2)) = 0; 
    [~,RRsamples] = max(a); % number of samples between cycles
    stats.RR = stats.sample_rate/RRsamples; % respiration rate in 1/s
    fprintf('Respiration rate is %0.3f breaths/sec or %0.1f breaths/min\n',stats.RR,stats.RR*60);
    if stats.RR<10/60 || stats.RR>25/60,
        warning('A respiratory rate of %1.1f breaths per minute seems unusual.',stats.RR*60);
    end

    % determine if this is a symmetric flow (i.e. if it is both inspiratory
    % and expiratory, or just one or the other, or something in the middle.
    crop = mod(numel(Q),RRsamples);            % leftover signal length after taking complete number of cycles out of signal
    Qcrop = Q(ceil(crop/2):end-floor(crop/2)); % ensure whole number of cycles to compare here
    Vin=sum(Qcrop(Qcrop>0))*sample_time;    % total volume in over whole number of cycles
    Vout=-sum(Qcrop(Qcrop<0))*sample_time;   % total volume out over whole number of cycles
    stats.Vin_frac = Vin/(Vin+Vout);
    fprintf('Inspiration volume fraction (*NOT* related to I:E ratio!) is %2.1f%%\n',stats.Vin_frac*100);
    if stats.Vin_frac>90/100,
        stats.mode = 'inspiration';
    elseif stats.Vin_frac>60/100,
        warning('It is not clear where in the circuit this flow meter is placed.');
        stats.mode = 'more_insp_than_exp';
    elseif stats.Vin_frac>40/100,
        stats.mode = 'bidirectional';
    elseif stats.Vin_frac>10/100,
        warning('It is not clear where in the circuit this flow meter is placed.');
        stats.mode = 'more_exp_than_insp';
    else
        stats.mode = 'expiration';
    end
    fprintf('Detected mode: %s\n',stats.mode);
        
    % partition the data according to flow direction
    Qpos = (Q>threshold);
    Qneg = (Q<-threshold);
    Qzero = ~(Qpos | Qneg);
    Q(Qzero)=0;  % ignore low flow, let's set to zero (helps identify sign changes too)
    
    % identify start and end indices of pos/neg sections
    sQ = sign(Q);     % i.e. sQ(Qpos)=1, sQ(Qzero)=0 and sQ(Qneg)=-1
    Qd = diff(sQ);    % identifies sign change locations
    Qdf = [0; Qd]; % this aligns sign changes with the data immediately after them
    Qdb = [Qd; 0]; % this aligns sign changes with the data immediately before them
    start_pos = (Qdf>0) & (sQ>0); % start of positive section
    start_neg = (Qdf<0) & (sQ<0); % start of negative section
    end_pos   = (Qdb<0) & (sQ>0); % end of positive section
    end_neg   = (Qdb>0) & (sQ<0); % end of negative section
    spi = find(start_pos); % start pos indices
    epi = find(end_pos); % etc
    sni = find(start_neg);
    eni = find(end_neg);
    
    % start on a positive cycle, and group indices into regions
    if numel(spi)<1 || numel(epi)<1,
        error('no complete positive flow sections detected in this data');
    end
    epi(epi<spi(1))=[];
    eni(eni<spi(1))=[];
    sni(sni<spi(1))=[];
    if numel(epi)<numel(spi)-1,
        error('something''s gone wrong, this error should not be possible.');
    end
    if numel(spi)>numel(epi), epi(end+1)=spi(end); end % we need the last spi marker to judge the end of the final cycle (not true of sni 7 lines below!)
    pos = [spi epi];
    
    % group negative indices into regions
    if numel(eni)<numel(sni)-1,
        error('something''s gone wrong, this error should not be possible.');
    end
    neg = [sni(1:numel(eni)) eni];
    
    % find start of first cycle
    % The start of a cycle is defined as start of largest contiguous region
    % of positive flow within 1/stats.RR seconds of previous cycle start.
    % (i.e. within one cycle time)
    indices = pos(:,1)-pos(1,1)<RRsamples;
    region_size = diff(pos(indices,:)');
    [~,beginning] = max(region_size);
    beginning = indices(beginning);
    % discard any regions before first cycle
    discard = neg(:,1)<pos(beginning,1);
    neg(discard,:)=[];
    discard = pos(:,1)<pos(beginning,1);
    pos(discard,:)=[];

    % now group into cycles
    cnt=0;
    unusedSegmentFlag=false;
    while ~isempty(pos),
        cnt = cnt+1;
        
        % find start of next cycle (largest positive region starting
        % between 95-150% of current cycle)
        interval = pos(:,1)-pos(1,1);
        indices = find(interval<RRsamples*1.5 & interval>RRsamples*0.95);
        region_size = diff(pos(indices,:)');
        [~,next_pos] = max(region_size);
        next_pos = indices(next_pos);

        % if we can't find the start of the next cycle then we haven't got
        % a complete cycle here so finish:
        if isempty(next_pos),
            break;
        end
        
        % otherwise, next contiguous postive section is to be recorded as an
        % inspiration cycle:
        cycles{n}(cnt).pos = pos(1,:);
        
        % record the largest negative region starting within one cycle as
        % the expiration cycle
        indices = neg(:,1)-pos(1,1)<RRsamples;
        region_size = diff(neg(indices,:)');
        [~,next_neg] = max(region_size);
        next_neg = indices(next_neg);
        cycles{n}(cnt).neg = neg(next_neg,:);
        
        % check for anything we didn't use in this cycle and remove
        % everything before next cycle from the list
        remove_pos = pos(:,1)<pos(next_pos,1);
        remove_neg = neg(:,1)<pos(next_pos,1);
        if sum(remove_pos)>1,
            unusedSegmentFlag=true;
        end
        if sum(remove_neg)>1,
            unusedSegmentFlag=true;
        end
        pos(remove_pos,:)=[];
        neg(remove_neg,:)=[];
    end
    
    % warn if there were discarded data segments (orphaned by zero flow
    % regions)
    if unusedSegmentFlag,
        warning('There were sections of nonzero flow orphaned from the contiguous flow regions. These have been discarded.');
    end
    
    % collect actual timeseries data in cycles too, applying margin at the same time.
    margin_samples = ceil(stats.sample_rate * margin); % this is hardcoded at 50 Hz because I'm lazy and this is probably all that's needed
    partitioned(n).p=[];partitioned(n).Q=[];partitioned(n).V=[];partitioned(n).t=[];partitioned(n).tau=[];
    partitioned(n).p_pos = []; partitioned(n).p_neg = []; partitioned(n).Q_pos = []; partitioned(n).Q_neg = [];
    partitioned(n).V_pos = []; partitioned(n).V_neg = []; partitioned(n).t_pos = []; partitioned(n).t_neg = [];
    partitioned(n).tau_pos = []; partitioned(n).tau_neg = []; % tau is the cycle time
    for ii = 1:numel(cycles{n})
        cycles{n}(ii).start_pos = cycles{n}(ii).pos(1,1);
        cycles{n}(ii).end_pos = cycles{n}(ii).pos(end,end);
        cycles{n}(ii).start = cycles{n}(ii).pos(1,1);
        cycles{n}(ii).end = cycles{n}(ii).pos(end,end); % will get overwritten below if there are any neg regions
        if numel(cycles{n}(ii).neg)>0,
            cycles{n}(ii).start_neg = cycles{n}(ii).neg(1,1);
            cycles{n}(ii).end_neg = cycles{n}(ii).neg(end,end);
            cycles{n}(ii).end = cycles{n}(ii).neg(end,end); % overwrites pos end
        else
            cycles{n}(ii).start_neg = [];
            cycles{n}(ii).end_neg = [];
        end
        if cycles{n}(ii).end_pos - cycles{n}(ii).start_pos > sum(margin_samples),
            cycles{n}(ii).start_pos_margin = cycles{n}(ii).start_pos + margin_samples(1);
            cycles{n}(ii).end_pos_margin = cycles{n}(ii).end_pos - margin_samples(2);
            cycles{n}(ii).start_neg_margin = cycles{n}(ii).start_neg + margin_samples(1);
            cycles{n}(ii).end_neg_margin = cycles{n}(ii).end_neg - margin_samples(2);
        else
            cycles{n}(ii).start_pos_margin = [];
            cycles{n}(ii).end_pos_margin = [];
            cycles{n}(ii).start_neg_margin = [];
            cycles{n}(ii).end_neg_margin = [];
        end
        cycles{n}(ii).indices = cycles{n}(ii).start:cycles{n}(ii).end;
        cycles{n}(ii).indices_pos = cycles{n}(ii).start_pos_margin:cycles{n}(ii).end_pos_margin;
        cycles{n}(ii).indices_neg = cycles{n}(ii).start_neg_margin:cycles{n}(ii).end_neg_margin;
        cycles{n}(ii).p = p(cycles{n}(ii).indices);
        cycles{n}(ii).Q = Q(cycles{n}(ii).indices);
        cycles{n}(ii).t = t(cycles{n}(ii).indices);
        cycles{n}(ii).V = V(cycles{n}(ii).indices);
        cycles{n}(ii).p_pos = p(cycles{n}(ii).indices_pos);
        cycles{n}(ii).Q_pos = Q(cycles{n}(ii).indices_pos);
        cycles{n}(ii).t_pos = t(cycles{n}(ii).indices_pos);
        cycles{n}(ii).V_pos = V(cycles{n}(ii).indices_pos);
        cycles{n}(ii).p_neg = p(cycles{n}(ii).indices_neg);
        cycles{n}(ii).Q_neg = Q(cycles{n}(ii).indices_neg);
        cycles{n}(ii).t_neg = t(cycles{n}(ii).indices_neg);
        cycles{n}(ii).V_neg = V(cycles{n}(ii).indices_neg);
        
        % amalgamated data:
        if ~isempty(cycles{n}(ii).t_pos),
            partitioned(n).p = [partitioned(n).p ; cycles{n}(ii).p_pos ; nan];
            partitioned(n).Q = [partitioned(n).Q ; cycles{n}(ii).Q_pos ; nan];
            partitioned(n).V = [partitioned(n).V ; cycles{n}(ii).V_pos ; nan];
            partitioned(n).t = [partitioned(n).t ; cycles{n}(ii).t_pos ; nan];
            partitioned(n).tau = [partitioned(n).tau ; cycles{n}(ii).t_pos-cycles{n}(ii).t_pos(1) ; nan];
            partitioned(n).tau_pos = [partitioned(n).tau_pos ; cycles{n}(ii).t_pos-cycles{n}(ii).t_pos(1) ; nan];
            partitioned(n).p_pos = [partitioned(n).p_pos ; cycles{n}(ii).p_pos ; nan];
            partitioned(n).Q_pos = [partitioned(n).Q_pos ; cycles{n}(ii).Q_pos ; nan];
            partitioned(n).V_pos = [partitioned(n).V_pos ; cycles{n}(ii).V_pos ; nan];
            partitioned(n).t_pos = [partitioned(n).t_pos ; cycles{n}(ii).t_pos ; nan];
        end
        if ~isempty(cycles{n}(ii).t_neg),
            partitioned(n).p = [partitioned(n).p ; cycles{n}(ii).p_neg ; nan];
            partitioned(n).Q = [partitioned(n).Q ; cycles{n}(ii).Q_neg ; nan];
            partitioned(n).V = [partitioned(n).V ; cycles{n}(ii).V_neg ; nan];
            partitioned(n).t = [partitioned(n).t ; cycles{n}(ii).t_neg ; nan];
            partitioned(n).tau = [partitioned(n).tau ; cycles{n}(ii).t_neg-cycles{n}(ii).t_pos(1) ; nan];
            partitioned(n).tau_neg = [partitioned(n).tau_neg ; cycles{n}(ii).t_neg-cycles{n}(ii).t_pos(1) ; nan];
            partitioned(n).p_neg = [partitioned(n).p_neg ; cycles{n}(ii).p_neg ; nan];
            partitioned(n).Q_neg = [partitioned(n).Q_neg ; cycles{n}(ii).Q_neg ; nan];
            partitioned(n).V_neg = [partitioned(n).V_neg ; cycles{n}(ii).V_neg ; nan];
            partitioned(n).t_neg = [partitioned(n).t_neg ; cycles{n}(ii).t_neg ; nan];
        end
    end
    
    % plotting
    if plt,
        QplotP=Q; QplotP(~Qpos)=nan;
        QplotN=Q; QplotN(~Qneg)=nan;
        pplotP=p; pplotP(~Qpos)=nan;
        pplotN=p; pplotN(~Qneg)=nan;
        fig=findobj('-regexp','tag',sprintf('ventilator-splitCycle%1.0d',n));
        if isempty(fig),
            fig=figure('tag',sprintf('ventilator-splitCycle%1.0d',n));
        else
            figure(fig);clf;
        end
        subplot(3,1,1); grid on; hold on;
        xlabel('time /s'); ylabel('pressure /cmH2O');
        plot(t,p,'k-');
        matchXaxes;
        subplot(3,1,2); grid on; hold on;
        xlabel('time /s'); ylabel('flow /(L/s)');
        plot(t,Q,'k-');
        matchXaxes;
        subplot(3,1,3); grid on; hold on;
        xlabel('time /s'); ylabel('volume /L)');
        plot(t,V,'k-');
        matchXaxes;

        for ii=1:numel(cycles{n}),
            subplot(3,1,1);
            plot(cycles{n}(ii).t_pos,cycles{n}(ii).p_pos,'r-','linewidth',1.5);
            plot(cycles{n}(ii).t_neg,cycles{n}(ii).p_neg,'b-','linewidth',1.5);
            plot(t([cycles{n}(ii).start_pos_margin cycles{n}(ii).end_pos_margin]),p([cycles{n}(ii).start_pos_margin cycles{n}(ii).end_pos_margin]),'ro','markerfacecolor','red');
            plot(t([cycles{n}(ii).start_neg_margin cycles{n}(ii).end_neg_margin]),p([cycles{n}(ii).start_neg_margin cycles{n}(ii).end_neg_margin]),'bo','markerfacecolor','blue');
            plot(t([cycles{n}(ii).start_pos cycles{n}(ii).end_pos]),p([cycles{n}(ii).start_pos cycles{n}(ii).end_pos]),'ro');
            plot(t([cycles{n}(ii).start_neg cycles{n}(ii).end_neg]),p([cycles{n}(ii).start_neg cycles{n}(ii).end_neg]),'bo');

            subplot(3,1,2);
            plot(cycles{n}(ii).t_pos,cycles{n}(ii).Q_pos,'r-','linewidth',1.5);
            plot(cycles{n}(ii).t_neg,cycles{n}(ii).Q_neg,'b-','linewidth',1.5);
            plot(t([cycles{n}(ii).start_pos_margin cycles{n}(ii).end_pos_margin]),Q([cycles{n}(ii).start_pos_margin cycles{n}(ii).end_pos_margin]),'ro','markerfacecolor','red');
            plot(t([cycles{n}(ii).start_neg_margin cycles{n}(ii).end_neg_margin]),Q([cycles{n}(ii).start_neg_margin cycles{n}(ii).end_neg_margin]),'bo','markerfacecolor','blue');
            plot(t([cycles{n}(ii).start_pos cycles{n}(ii).end_pos]),Q([cycles{n}(ii).start_pos cycles{n}(ii).end_pos]),'ro');
            plot(t([cycles{n}(ii).start_neg cycles{n}(ii).end_neg]),Q([cycles{n}(ii).start_neg cycles{n}(ii).end_neg]),'bo');

            subplot(3,1,3);
            plot(cycles{n}(ii).t_pos,cycles{n}(ii).V_pos,'r-','linewidth',1.5);
            plot(cycles{n}(ii).t_neg,cycles{n}(ii).V_neg,'b-','linewidth',1.5);
            plot(t([cycles{n}(ii).start_pos_margin cycles{n}(ii).end_pos_margin]),V([cycles{n}(ii).start_pos_margin cycles{n}(ii).end_pos_margin]),'ro','markerfacecolor','red');
            plot(t([cycles{n}(ii).start_neg_margin cycles{n}(ii).end_neg_margin]),V([cycles{n}(ii).start_neg_margin cycles{n}(ii).end_neg_margin]),'bo','markerfacecolor','blue');
            plot(t([cycles{n}(ii).start_pos cycles{n}(ii).end_pos]),V([cycles{n}(ii).start_pos cycles{n}(ii).end_pos]),'ro');
            plot(t([cycles{n}(ii).start_neg cycles{n}(ii).end_neg]),V([cycles{n}(ii).start_neg cycles{n}(ii).end_neg]),'bo');
        end
        sgtitle(sprintf('meter %1.0d, margin = %0.2f s (start), %0.2f s (end)',n,margin(1),margin(2)));
    end
end
