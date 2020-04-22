function models = identifyRestriction2(fn,margin,withFig,c_selection)
% function to estimate restriction from two .sig files, [fn '.sig'],
% located in 'FlowMeter1/' and 'FlowMeter2/' respectively.
%
% If fn is a cell array of strings, it treats all of the datasets together
% for the purposes of fitting the parameters.
%
% margin is either a scalar defining the time interval to exclude from both
% the start and end of each inspiration and each expiration phase, or a
% 2-element vector defining the time to remove from start and end,
% respectively.
%
% c_selection is a cell array of models to include (see output, below).
%             Default is all models. (only relevant to plotting).
%             e.g. {'mixed','turbulent','linear'}. If a string is given it
%             is automatically converted to a 1x1 cell array. Combinations
%             can also be supported, so far:
%             'standard' = {'mixed','turbulent','linear'}
%
% If no output arguments are given then it plots results. If arguments are
% given and you still want figures, set the third input argument as the
% string 'withFig'.
%
% Examples:
% identifyRestriction('v01',0.2)
% identifyRestriction({'v01','v02','v03'},[0.2 0.2])
%
% output is a structure with different models fitted:
% models(n).label = 'mixed'     quadratic fit (mixed flow model, 2000<Re<4000), constrained to origin
% models(n).label = 'turbulent' squared fit (turbulent flow model, Re>4000), constrained to origin
% models(n).label = 'laminar'   linear fit (laminar flow model, Re<2000), constrained to origin
% models(n).label = 'valve_lin' check valve model; linear fit, zero-flow pressure offset
% models(n).c2 quadratic loss coefficient
% models(n).c1 linear loss coefficient
% models(n).c0 constant loss coefficient
% models(n).c  vector of all coefficients (lowest to highest order)
%
% Output units are L for volume, cmH2O for pressure, s for time, and
% permutations thereof. No exceptions. No mL or lpm or 1/min or anything
% like that.
%
% This uses Method 2, based on phase-averaged signals. Method 1 is also
% available, which uses alignment of full timeseries.
%

% setup
if nargin<3, withFig = 'noThanks'; end
if nargout>0 && ~strcmp(withFig,'withFig'), % argument to disable plotting if numerical output is requested, unless figs are specifically requested in the input args
    plt = false;
else
    plt = true;
end
maxme=true; % maximise all multi-plot figures?
if nargin<2 || isempty(margin),
    warning('No margin specified; defaulting to 0.2 seconds at start and end.');
    margin = 0.2;
end
if nargin<4,
    c_selection = [];
end
if ischar(c_selection),
    c_selection = {c_selection};
end
if numel(margin)<2,
    margin = [1 1] * margin;
end
start_crop = margin(1);
end_crop = margin(2);
Qthreshold = 0;

% added functionality to allow several datasets to be included in fitting
if ~iscell(fn),
    fns={fn};
else
    fns=fn;
end

% fudge because flow rates too low for v111 for this to work
% (do something about it in importWrapper??)
% We need a nicer way of dealing with this kind of problematic input data in general.
for ii=1:numel(fns),
    if strcmp(fns{ii},'v111'),
        warning('v111 is not going to be treated well by splitCycle, due to low-flow cutout. It is being excluded from results here.');
        fns(ii)=[];
        if numel(fns)<1,
            models = [];
            return 
        end
    end
end

% load data
[s,par] = importWrapper(fns{1});
n_meters = numel(s);
if n_meters<2,
    error('This function needs two meters to estimate restriction between them');
end

% ************************************
% **** METHOD 2 **********************
% ************************************

% estimate restriction coefficients between two meters
% First, split signals into phases, using splitCycle, then do phase
% averaging, and compare mean cycles.
[s_part,stats,cycles] = splitCycle(s,[],margin);
fprintf('Number of cycles for comparison: %3.0d (meter 1), %3.0d (meter 2)\n',numel(cycles{1}),numel(cycles{2}));
interp_method='spline';
% determine common start and end times (cropping signals where necessary)
starttime=-9e9; endtime=9e9;
starttime_pos=-9e9; endtime_pos=9e9;
starttime_neg=-9e9; endtime_neg=9e9;
for n = 1:n_meters,
    for ii=1:numel(cycles{n}),
        starttime = max(starttime,min(cycles{n}(ii).tau));
        endtime   = min(endtime  ,max(cycles{n}(ii).tau));
        starttime_pos = max(starttime_pos,min(cycles{n}(ii).tau_pos)); % not sure which I'll need at time of writing (I'd like to say I plan more...)
        endtime_pos   = min(endtime_pos  ,max(cycles{n}(ii).tau_pos));
        starttime_neg = max(starttime_neg,min(cycles{n}(ii).tau_neg));
        endtime_neg   = min(endtime_neg  ,max(cycles{n}(ii).tau_neg));
    end
end
% resample just in case (this is with half a mind to provide
% inter-timestep alignment within splitCycle, which would adjust the
% tau values for each cycle while leaving the datapoints unchanged.)
% And then phase-average.
numsamples = floor((endtime-starttime)*stats.sample_rate)+1;
newtau = linspace(starttime,endtime,numsamples);
numsamples_pos = floor((endtime_pos-starttime_pos)*stats.sample_rate)+1;
newtau_pos = linspace(starttime_pos,endtime_pos,numsamples_pos);
if ~isempty(starttime_neg), % not sure why I'm bothering; pretty sure I'm only going to use positive flow in this function
    numsamples_neg = floor((endtime_neg-starttime_neg)*stats.sample_rate)+1;
    newtau_neg = linspace(starttime_neg,endtime_neg,numsamples_neg);
end
for n = 1:n_meters,
    phase_av(n).p=zeros(numel(newtau),1);phase_av(n).Q=zeros(numel(newtau),1);phase_av(n).V=zeros(numel(newtau),1);
    phase_av(n).p_pos=zeros(numel(newtau_pos),1);phase_av(n).Q_pos=zeros(numel(newtau_pos),1);phase_av(n).V_pos=zeros(numel(newtau_pos),1);
    if ~isempty(starttime_neg), % not sure why I'm bothering; pretty sure I'm only going to use positive flow in this function
        phase_av(n).p_neg=zeros(numel(newtau_neg),1);phase_av(n).Q_neg=zeros(numel(newtau_neg),1);phase_av(n).V_neg=zeros(numel(newtau_neg),1);
    end
    for ii=1:numel(cycles{n}),
        % resample
        cycles{n}(ii).p=interp1(cycles{n}(ii).tau,cycles{n}(ii).p,newtau,interp_method);
        cycles{n}(ii).Q=interp1(cycles{n}(ii).tau,cycles{n}(ii).Q,newtau,interp_method);
        cycles{n}(ii).V=interp1(cycles{n}(ii).tau,cycles{n}(ii).V,newtau,interp_method);
        cycles{n}(ii).tau = newtau;
        cycles{n}(ii).p_pos=interp1(cycles{n}(ii).tau,cycles{n}(ii).p,newtau_pos,interp_method);
        cycles{n}(ii).Q_pos=interp1(cycles{n}(ii).tau,cycles{n}(ii).Q,newtau_pos,interp_method);
        cycles{n}(ii).V_pos=interp1(cycles{n}(ii).tau,cycles{n}(ii).V,newtau_pos,interp_method);
        cycles{n}(ii).tau_pos = newtau_pos;
        % and phase average
        phase_av(n).p = phase_av(n).p + cycles{n}(ii).p(:);
        phase_av(n).Q = phase_av(n).Q + cycles{n}(ii).Q(:);
        phase_av(n).V = phase_av(n).V + cycles{n}(ii).V(:);
        phase_av(n).p_pos = phase_av(n).p_pos + cycles{n}(ii).p_pos(:);
        phase_av(n).Q_pos = phase_av(n).Q_pos + cycles{n}(ii).Q_pos(:);
        phase_av(n).V_pos = phase_av(n).V_pos + cycles{n}(ii).V_pos(:);
        if ~isempty(starttime_neg), % not sure why I'm bothering; pretty sure I'm only going to use positive flow in this function
            % resample
            cycles{n}(ii).p_neg=interp1(cycles{n}(ii).tau,cycles{n}(ii).p,newtau_neg,interp_method);
            cycles{n}(ii).Q_neg=interp1(cycles{n}(ii).tau,cycles{n}(ii).Q,newtau_neg,interp_method);
            cycles{n}(ii).V_neg=interp1(cycles{n}(ii).tau,cycles{n}(ii).V,newtau_neg,interp_method);
            cycles{n}(ii).tau_neg = newtau_neg;
            % and phase average
            phase_av(n).p_neg = phase_av(n).p_neg + cycles{n}(ii).p_neg(:);
            phase_av(n).Q_neg = phase_av(n).Q_neg + cycles{n}(ii).Q_neg(:);
            phase_av(n).V_neg = phase_av(n).V_neg + cycles{n}(ii).V_neg(:);
        end
    end
    phase_av(n).t = newtau;
    phase_av(n).t_pos = newtau_pos;
    phase_av(n).p = phase_av(n).p/numel(cycles{n});
    phase_av(n).Q = phase_av(n).Q/numel(cycles{n});
    phase_av(n).V = phase_av(n).V/numel(cycles{n});
    phase_av(n).p_pos = phase_av(n).p_pos/numel(cycles{n});
    phase_av(n).Q_pos = phase_av(n).Q_pos/numel(cycles{n});
    phase_av(n).V_pos = phase_av(n).V_pos/numel(cycles{n});
    if ~isempty(starttime_neg), % not sure why I'm bothering; pretty sure I'm only going to use positive flow in this function
        phase_av(n).t_neg = newtau_neg;
        phase_av(n).p_neg = phase_av(n).p_neg/numel(cycles{n});
        phase_av(n).Q_neg = phase_av(n).Q_neg/numel(cycles{n});
        phase_av(n).V_neg = phase_av(n).V_neg/numel(cycles{n});
    end
end
allcycles{1} = cycles; % for plotting

% this can add more than one dataset into the fit:
s_corr_markers = [0 0];
if numel(fns)>1,
    for jj=2:numel(fns),
        s2 = importWrapper(fns{jj});
        % First, split signals into phases, using splitCycle, then do phase
        % averaging, and compare mean cycles.
        [~,stats,cycles] = splitCycle(s2,[],margin);
        fprintf('Number of cycles for comparison: %3.0d (meter 1), %3.0d (meter 2) [case: %s]\n',numel(cycles{1}),numel(cycles{2}),fns{jj});
        % determine common start and end times (cropping signals where necessary)
        starttime_pos=-9e9; endtime_pos=9e9;
        for n = 1:n_meters,
            for ii=1:numel(cycles{n}),
                starttime_pos = max(starttime_pos,min(cycles{n}(ii).tau_pos));
                endtime_pos   = min(endtime_pos  ,max(cycles{n}(ii).tau_pos));
            end
        end
        % resample just in case (this is with half a mind to provide
        % inter-timestep alignment within splitCycle, which would adjust the
        % tau values for each cycle while leaving the datapoints unchanged.)
        % And then phase-average.
        numsamples_pos = floor((endtime_pos-starttime_pos)*stats.sample_rate)+1;
        newtau_pos = linspace(starttime_pos,endtime_pos,numsamples_pos);
        s_corr_markers(end+1,1:2) = [numel(phase_av(1).t_pos) numel(phase_av(2).t_pos)]; % already resigned to just using positive flow here
        for n = 1:n_meters,
            this_phase_av(n).p_pos=zeros(numel(newtau_pos),1);this_phase_av(n).Q_pos=zeros(numel(newtau_pos),1);this_phase_av(n).V_pos=zeros(numel(newtau_pos),1);
            for ii=1:numel(cycles{n}),
                % resample
                cycles{n}(ii).p_pos=interp1(cycles{n}(ii).tau,cycles{n}(ii).p,newtau_pos,interp_method);
                cycles{n}(ii).Q_pos=interp1(cycles{n}(ii).tau,cycles{n}(ii).Q,newtau_pos,interp_method);
                cycles{n}(ii).V_pos=interp1(cycles{n}(ii).tau,cycles{n}(ii).V,newtau_pos,interp_method);
                cycles{n}(ii).tau_pos = newtau_pos;
                % and phase average
                this_phase_av(n).p_pos = this_phase_av(n).p_pos + cycles{n}(ii).p_pos(:);
                this_phase_av(n).Q_pos = this_phase_av(n).Q_pos + cycles{n}(ii).Q_pos(:);
                this_phase_av(n).V_pos = this_phase_av(n).V_pos + cycles{n}(ii).V_pos(:);
                % this is to prepare for plotting all cases alongside phase
                % average:
                cycles{n}(ii).tau = cycles{n}(ii).tau + phase_av(n).t_pos(end) + phase_av(n).t_pos(1)*2; % use start time as well to offset end time - this is a fudge assuming start and end are similar
            end
            phase_av(n).t_pos = [phase_av(n).t_pos(:); newtau_pos(:) + phase_av(n).t_pos(end) + phase_av(n).t_pos(1)*2]; % use start time as well to offset end time - this is a fudge assuming start and end are similar
            phase_av(n).p_pos = [phase_av(n).p_pos; this_phase_av(n).p_pos/numel(cycles{n})];
            phase_av(n).Q_pos = [phase_av(n).Q_pos; this_phase_av(n).Q_pos/numel(cycles{n})];
            phase_av(n).V_pos = [phase_av(n).V_pos; phase_av(n).V_pos/numel(cycles{n})];
        end
        allcycles{jj} = cycles; % for plotting
    end
end
s_corr_markers(end+1,1:2) = [numel(phase_av(1).t_pos) numel(phase_av(2).t_pos)]; % already resigned to just using positive flow here

% and compare deltaP-Qbar (pressure drop - mean flow)
phase_av(3).p_pos = phase_av(1).p_pos-phase_av(2).p_pos;
phase_av(3).Q_pos = (phase_av(1).Q_pos+phase_av(2).Q_pos)/2;
if ~isempty(starttime_neg), % not sure why I'm bothering; pretty sure I'm only going to use positive flow in this function
    phase_av(3).p_neg = phase_av(1).p_neg-phase_av(2).p_neg;
    phase_av(3).Q_neg = (phase_av(1).Q_neg+phase_av(2).Q_neg)/2;
end

% fitting: default is quadratic with origin constraint.
% (only doing positive flow)
x = phase_av(3).Q_pos;  % the abscissa for the curve fit (we're using positive flow only)
y = phase_av(3).p_pos;  % the ordinate for the curve fit

c(1,1:3) = fitFunc(x,y,[2 1]); % quad_origin ('mixed')
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(1,1:3));
c(2,1:3) = fitFunc(x,y,[2]); % squared ('turbulent')
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(1,1:3));
c(3,1:3) = fitFunc(x,y,[1]);   % lin_origin ('laminar')
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(2,1:3));
c(4,1:3) = fitFunc(x,y,[1 0]);   % linear (for valve cracking pressure) ('valve')
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(2,1:3));


% plot p-Q relationship
if plt
    % selecting which fits to plot
    if isempty(c_selection), 
        c_sel = [1:4];
    else
        c_sel =[];
        for ii=1:numel(c_selection),
            switch c_selection{ii},
                case 'mixed',
                    c_sel(end+1) = 1;
                case 'turbulent',
                    c_sel(end+1) = 2;
                case 'laminar',
                    c_sel(end+1) = 3;
                case 'valve_lin',
                    c_sel(end+1) = 4;
                case 'standard',
                    c_sel(end+1) = 1;
                    c_sel(end+1) = 2;
                    c_sel(end+1) = 3;
                otherwise
                    warning('unrecognised model type in c_selection input');
            end
        end
    end
        
    % this is a bit of a fudge (leftover from identifyRestriction1...)
    for ii=1:3,
        s_corr(ii).p = phase_av(ii).p_pos;
        s_corr(ii).Q = phase_av(ii).Q_pos;
        s_corr(ii).V = phase_av(ii).V_pos;
        s_corr(ii).t = phase_av(ii).t_pos;
        s_corr(ii).label = s(1).label;
    end
    % end fudge
    plotRestriction(s_corr,c(c_sel,:),maxme,s_corr_markers,fns,allcycles);
end

% assign outputs
models(1).label = 'mixed';
models(1).c2 = c(1,1);
models(1).c1 = c(1,2);
models(1).c0 = c(1,3);
models(1).c = c(1,end:-1:1);
models(2).label = 'turbulent';
models(2).c2 = c(2,1);
models(2).c1 = c(2,2);
models(2).c0 = c(2,3);
models(2).c = c(2,end:-1:1);
models(3).label = 'laminar';
models(3).c2 = c(3,1);
models(3).c1 = c(3,2);
models(3).c0 = c(3,3);
models(3).c = c(3,end:-1:1);
models(4).label = 'valve_lin';
models(4).c2 = c(4,1);
models(4).c1 = c(4,2);
models(4).c0 = c(4,3);
models(4).c = c(4,end:-1:1);



