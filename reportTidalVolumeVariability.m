function [rep,TV]=reportTidalVolumeVariability(prefix,range,threshold,margin,numcycles)
% reports on tidal volume stats for range of cases

% input validation & defualts
if nargin<5,
    numcycles=[];
end
if nargin<3 || isempty(threshold), % threshold below which flow is considered zero
    threshold = 0.01;
end
if nargout==0, % plot figure(s) if no outputs are requested
    plt=true;
else           % otherwise, don't plot
    plt=false;
end
if nargin<4 || isempty(margin),
    margin = 0; % seconds
    warning('No margin specified; defaulting to 0 seconds at start and end.');
end
if numel(margin)<2,
    margin = margin * [1 1]; % start and end margins
end

% assemble case names from given range
for ii=1:numel(range)
    fns{ii}=sprintf('%s%02.0f',prefix,range(ii));
end

rep=[]; % save output report for function output, or copying to clipboard etc.
rep = [rep sprintf('Case\tMeter\tNo. of cycles\tMean Tidal Vol. (L)\tSt. Dev.\tInterquart. range\n')];

for ii=1:numel(fns)
    [sig,par] = importWrapper(fns{ii});
    [p,st,cy]=splitCycle(sig,threshold,margin);
    if plt,
        splitCycle(sig,threshold,margin);
    end
    
    % get tidal volumes
    for n=1:numel(cy), 
        TVs = [];
        for jj=1:numel(cy{n}),
            TVs(jj) = max(cy{n}(jj).V)-min(cy{n}(jj).V);
        end
        if plt,
            figure(13+n);
            bar(TVs);
        end

        % aggregate stats
        if isempty(numcycles), % this decides number of cycles to include
            rng = numel(TVs);
        else
            rng = min(numcycles,numel(TVs));
        end
        TV(ii,n).num  = rng;
        TV(ii,n).mean = mean(TVs(1:rng));
        TV(ii,n).std  = std(TVs(1:rng));
        TV(ii,n).iqr  = iqr(TVs(1:rng));
    
        % output in tab-delimited format:
        rep = [rep sprintf('%s\t%1.0d\t%1.0f\t%1.3f\t%3.3e\t%3.3e\n',fns{ii}, n, TV(ii,n).num, TV(ii,n).mean, TV(ii,n).std, TV(ii,n).iqr)];
    end
        
    if plt,
        pause;
    end
end

% output in tab-delimited format:
fprintf(rep);
clipboard('copy',rep);
fprintf('\n Table copied to clipboard.\n');
