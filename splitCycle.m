function [cropped,cycles]=splitCycle(s,threshold,margin,plt)
%
% This is a bit rough and ready for now. In simplest usage:
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
% Don't get too attached to the format of the data output from this
% function; I plan to revamp it in the next day or so.
%
%
% Some more somewhat cryptic documentation below. To be updated.
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
if nargin<2 || isempty(threshold),
    threshold = 0.01;
end
if nargin<4,
    if nargout==0, % by default, plot if no outputs are requested
        plt=true;
    else           % otherwise, don't plot
        plt=false;
    end
end
if nargin<3 || isempty(margin),
    margin = 0; % seconds
end
if numel(margin)<2,
    margin = margin * [1 1];
end

%  number of meters saved in the dataset
n_meters = numel(s);

% so the figures for unavailable meters get closed
for n = 1:2,
    if ishandle(100*n+69) && n>n_meters, close(100*n+69); end
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
        % (allowing for numerical rounding). This could potentially be
        % relaxed to a smaller percentage, then subtly irregular sampling
        % could just be corrected by resampling here.
        error('Time series must be sampled at regular intervals. Resample or interpolate to get a useable signal.');
    end
    sample_rate = 1/sample_time;
    
    % find respiration rate using autocorrelation
    % NB it would be nice to interpolate to get an even closer fit - but
    % not until I've got the cycles plotting nicely over the top of each
    % other to get an indication of whether it improves anything!
    a = autocorr(Q,'numlags',numel(Q)-1);
    % now set first 0.5 s of correlation data to zero to avoid identifying
    % a period of zero. This will only start to cause problems if the RR is
    % well over 60 - probably unheard of. (We could also consider using the
    % gradient, to exclude anything before the first trough but that would
    % be susceptible to errors in the presence of noisy signals. This shold
    % be fine for what we're doing.)
    a(1:ceil(sample_rate/2)) = 0; 
    [~,RRsamples] = max(a); % number of samples between cycles
    RR = sample_rate/RRsamples; % respiration rate in 1/s
    if RR<10/60 || RR>25/60,
        warning('A respiratory rate of %1.1f breaths per minute seems unusual.',RR*60);
    end

    % determine if this is a symmetric flow (i.e. if it is both inspiratory
    % and expiratory, or just one or the other, or something in the middle.
    crop = mod(numel(Q),RRsamples);            % leftover signal length after taking complete number of cycles out of signal
    Qcrop = Q(ceil(crop/2):end-floor(crop/2)); % ensure whole number of cycles to compare here
    Vin=sum(Qcrop(Qcrop>0))*sample_time;    % total volume in over whole number of cycles
    Vout=-sum(Qcrop(Qcrop<0))*sample_time;   % total volume out over whole number of cycles
    Vfrac = Vin/(Vin+Vout); 
    fprintf('Inspiration volume fraction (*NOT* related to I:E ratio!) is %2.1f%%\n',Vfrac*100);
    if Vfrac>90/100,
        leg = 'inspiration';
    elseif Vfrac>60/100,
        warning('It is not clear where in the circuit this flow meter is placed.');
        leg = 'more_insp_than_exp';
    elseif Vfrac>40/100,
        leg = 'bidirectional';
    elseif Vfrac>10/100,
        warning('It is not clear where in the circuit this flow meter is placed.');
        leg = 'more_exp_than_insp';
    else
        leg = 'expiration';
    end
    fprintf('Detected mode: %s\n',leg);
    
    % find the I:E ratio (this is not working as I'd hoped - shelved for now) 
    % First, low pass filter at 8 times the respiration rate, to smooth any
    % noise near the flow reversal; I'm really hoping for only one zero
    % crossing in the smoothed curve.
%     Qsmooth = lpfilt(Q,1/sample_rate,1/RR/32); 
%     dQ = diff(Q); ddQ = diff(dQ); ddQ=[0;ddQ(:)];
%     figure; subplot(3,1,1);
%     plot([Q(:) Qsmooth(:)]); grid on;
%     setMatchXaxes;
%     subplot(3,1,2);
%     plot(dQ);grid on;
%     setMatchXaxes;
%     subplot(3,1,3);
%     plot(ddQ);grid on;
%     setMatchXaxes;
%     keyboard
%     return
    
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
    
    % find start of cycle (using intervals between postive
    % regions) and discard any regions before this
    interval = pos(2:end,1)-pos(1:end-1,2);
    beginning = find(interval>0.95*max(interval));
    beginning = beginning(1)+1;
    if pos(beginning,1)>1.1*RRsamples, % this means we started in the middle of an interval, with a 10% margin
        [~,beginning] = min((pos(:,1)-pos(beginning,1)+RRsamples).^2);
    end
    discard = neg(:,1)<pos(beginning,1);
    neg(discard,:)=[];
    discard = pos(:,1)<pos(beginning,1);
    pos(discard,:)=[];
    
    % now group into cycles
    cnt=0;
    while ~isempty(pos),
        cnt = cnt+1;
        
        % add positive sections
        indices = pos(:,1)-pos(1,1) < RRsamples*0.95;
        cycles{n}(cnt).pos = pos(indices,:);
        pos(indices,:)=[];
        
        if isempty(pos), % no full cycles left
            cycles{n}(cnt)=[];
            break;
        else
            % check there are no negative sections overlapping (I don't
            % know what to do with those)
            if any(any(neg<=cycles{n}(cnt).pos(end,end))),
                error('The positive and negative regions in this cycle overlap. Not sure how to treat these.');
                % could just stick the overlapping bits in the positive side???
            end
            
            % now add the negative sections
            indices = neg(:,1) < pos(1,1);
            cycles{n}(cnt).neg = neg(indices,:);
            neg(indices,:)=[];
        end
    end
    
    % collect actual timeseries in cycles too, applying margin at the same time.
    margin_samples = ceil(50 * margin); % this is hardcoded at 50 Hz because I'm lazy and this is probably all that's needed
    cropped(n).p=[];cropped(n).Q=[];cropped(n).V=[];cropped(n).t=[];cropped(n).tau=[];
    cropped(n).p_pos = []; cropped(n).p_neg = []; cropped(n).Q_pos = []; cropped(n).Q_neg = [];
    cropped(n).V_pos = []; cropped(n).V_neg = []; cropped(n).t_pos = []; cropped(n).t_neg = [];
    cropped(n).tau_pos = []; cropped(n).tau_neg = []; % tau is the cycle time
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
            cropped(n).p = [cropped(n).p ; cycles{n}(ii).p_pos ; nan];
            cropped(n).Q = [cropped(n).Q ; cycles{n}(ii).Q_pos ; nan];
            cropped(n).V = [cropped(n).V ; cycles{n}(ii).V_pos ; nan];
            cropped(n).t = [cropped(n).t ; cycles{n}(ii).t_pos ; nan];
            cropped(n).tau = [cropped(n).tau ; cycles{n}(ii).t_pos-cycles{n}(ii).t_pos(1) ; nan];
            cropped(n).tau_pos = [cropped(n).tau_pos ; cycles{n}(ii).t_pos-cycles{n}(ii).t_pos(1) ; nan];
            cropped(n).p_pos = [cropped(n).p_pos ; cycles{n}(ii).p_pos ; nan];
            cropped(n).Q_pos = [cropped(n).Q_pos ; cycles{n}(ii).Q_pos ; nan];
            cropped(n).V_pos = [cropped(n).V_pos ; cycles{n}(ii).V_pos ; nan];
            cropped(n).t_pos = [cropped(n).t_pos ; cycles{n}(ii).t_pos ; nan];
        end
        if ~isempty(cycles{n}(ii).t_neg),
            cropped(n).p = [cropped(n).p ; cycles{n}(ii).p_neg ; nan];
            cropped(n).Q = [cropped(n).Q ; cycles{n}(ii).Q_neg ; nan];
            cropped(n).V = [cropped(n).V ; cycles{n}(ii).V_neg ; nan];
            cropped(n).t = [cropped(n).t ; cycles{n}(ii).t_neg ; nan];
            cropped(n).tau = [cropped(n).tau ; cycles{n}(ii).t_neg-cycles{n}(ii).t_pos(1) ; nan];
            cropped(n).tau_neg = [cropped(n).tau_neg ; cycles{n}(ii).t_neg-cycles{n}(ii).t_pos(1) ; nan];
            cropped(n).p_neg = [cropped(n).p_neg ; cycles{n}(ii).p_neg ; nan];
            cropped(n).Q_neg = [cropped(n).Q_neg ; cycles{n}(ii).Q_neg ; nan];
            cropped(n).V_neg = [cropped(n).V_neg ; cycles{n}(ii).V_neg ; nan];
            cropped(n).t_neg = [cropped(n).t_neg ; cycles{n}(ii).t_neg ; nan];
        end
    end
    
    % plotting
    if plt,
        QplotP=Q; QplotP(~Qpos)=nan;
        QplotN=Q; QplotN(~Qneg)=nan;
        pplotP=p; pplotP(~Qpos)=nan;
        pplotN=p; pplotN(~Qneg)=nan;
        figure(100*n+69); clf;
        subplot(3,1,1); grid on; hold on;
        xlabel('time /s'); ylabel('pressure /cmH2O');
        plot(t,p,'k-');
        matchXaxes;
 %       plot(cropped(n).t,cropped(n).p,'r-');
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
