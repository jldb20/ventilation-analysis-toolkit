function models = identifyRestriction(fn,margin,withFig,c_selection)
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
% This uses Method 1, based on aligning full timeseries. A Method 2 is in
% development at the time of writing, which uses phase-averaged signals
% from the datasets.
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

% load data
[s,par] = importWrapper(fns{1});
n_meters = numel(s);
if n_meters<2,
    error('This function needs two meters to estimate restriction between them');
end

% ************************************
% **** METHOD 1 **********************
% ************************************

% estimate restriction coefficients between two meters
% Use correlation to align signals from two meters
% Previous versions have tried to do something clever with identifying
% which way around the meters are connected, but this is a bad idea. This
% version will assume the data has already been processed in
% importWrapper() so the sign convention is consistent between the two
% datasets.
cycle_length = 50 * 60/10 ; % 50Hz, 60 seconds, min. 10 breaths/min (=>max. cycle length)
s_corr = doCorrelate(s,cycle_length,false);

% fudge because flow rates too low for v111 for this one (see fudge in
% loop below too. mmmm fudge.)
% (just do this in importWrapper?)
if strcmp(fns{1},'v111'),
    s_corr(1).Q = s_corr(2).Q;
end

% this can add more than one dataset into the fit:
s_corr_markers = [0 0];
if numel(fns)>1,
    for ii=2:numel(fns),
        s2{ii} = importWrapper(fns{ii});
        s2_corr{ii} = doCorrelate(s2{ii},cycle_length,false);
        if strcmp(fns{ii},'v111'), s2_corr{ii}(1).Q = s2_corr{ii}(2).Q; end % fudge
        for jj=1:2
            s_corr_markers(ii,jj) = numel(s_corr(jj).t);
            s_corr(jj).p = [s_corr(jj).p; s2_corr{ii}(jj).p];
            s_corr(jj).Q = [s_corr(jj).Q; s2_corr{ii}(jj).Q];
            s_corr(jj).t = (0:numel(s_corr(jj).p)-1)/50; % 50Hz, easy route.
        end
    end
end
s_corr_markers(end+1,1:2) = [numel(s_corr(1).t) numel(s_corr(2).t)];

% and compare deltaP-Qbar (pressure drop - mean flow)
s_corr(3).p = s_corr(1).p-s_corr(2).p;
s_corr(3).Q = (s_corr(1).Q+s_corr(2).Q)/2;

% We need to omit the data at the very start and end of the inspiratory
% phase, where the pressure is rising/falling rapidly, because this
% data is (a) ill-conditioned and (b) contains unwanted dynamic effects.
Q1zero = abs(s_corr(1).Q)<=Qthreshold;
Q2zero = abs(s_corr(2).Q)<=Qthreshold;
Qzero = Q1zero | Q2zero; % indicates where either flow rate is below threshold
s_corr(3).Q(Qzero) = 0; % only include data with nonzero meter flow readings
start_flow = find(diff(Qzero)<0); % NB this identifies the flow magnitude rising above the threshold in *either* direction (pos or neg)
end_flow = find(diff(Qzero)>0); % NB this identifies the flow magnitude dropping below the threshold from *either* direction (pos or neg)
data_length = numel(s_corr(3).Q);
num_samples = ceil(start_crop * 50); % 50Hz
for jj=1:numel(start_flow),
    window = start_flow(jj)+(1:num_samples);
    window(window<1 | window>data_length) = [];
    s_corr(3).Q(window)=0;
end
num_samples = ceil(end_crop * 50); % 50Hz
end_flow = end_flow - num_samples;
for jj=1:numel(end_flow),
    window = end_flow(jj)+(1:num_samples);
    window(window<1 | window>data_length) = [];
    s_corr(3).Q(window)=0;
end

% fitting: default is quadratic with origin constraint.
pos_d = find(s_corr(3).Q>Qthreshold); % (we're fitting positive flow only)
x = s_corr(3).Q(pos_d);  % the abscissa for the curve fit (we're using positive flow only)
y = s_corr(3).p(pos_d);  % the ordinate for the curve fit
c(1,1:3) = fitFunc(x,y,[2 1]); % quad_origin
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(1,1:3));
c(2,1:3) = fitFunc(x,y,[2]); % squared
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(1,1:3));
c(3,1:3) = fitFunc(x,y,[1]);   % lin_origin
fprintf('loss coeffs between meters: c2 = %0.3e cmH2O/(L/s)^2   c1 = %0.3e cmH2O/(L/s)   c0 = %0.3e cmH2O\n',c(2,1:3));
c(4,1:3) = fitFunc(x,y,[1 0]);   % linear (for valve cracking pressure)
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
        
    plotRestriction(s_corr,c(c_sel,:),maxme,s_corr_markers,fns);
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



