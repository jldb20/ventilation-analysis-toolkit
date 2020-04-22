function models=restrictorSweep(prefix,range,doPics,method,margin)
% Sweeps a range of test cases to produce estimates of the restriction
% between two meters. Can also amalagamate test sets into groups for the
% purpose of restriction parameter estimation.
%
% prefix is the test case prefix (e.g. 'M', 'v')
% range is either a vector of test case numbers, of a cell array of vectors.
%       If it's a vector then each case will be considered separately. If
%       it's a cell array then each cell should contain a vector and each
%       set of cases described by a vector will be treated as a group for
%       the purposes of parameter identification. e.g. [1:3], or
%       {1:3,2:5,6:9}
% doPics is boolean to decide whether to produce graphs in ./graphs/ folder
% method, if given, is a function handle for the identification algorithm
%         to use. Defaults to @identifyRestriction2. There's also
%         @identifyRestriction available.
% margin is either a scalar defining the time interval to exclude from both
%        the start and end of each inspiration and each expiration phase,
%        or a 2-element vector defining the time to remove from start and
%        end, respectively.
%
% models is a cell array of model outputs from identifyRestriction
%

if nargin<5 || isempty(margin),
    warning('No margin specified; defaulting to 0.2 seconds at start and end.');
    margin = 0.2;
end
if nargin<3 || isempty(doPics),
    doPics=false;
end
if doPics, % controls whether identifyRestriction will plot or not.
    withFig='withFig';
else
    withFig='noThanksVeryMuch';
end
if nargin<4 || isempty(method),
    method = @identifyRestriction2;
end

% if only a vector is given, then rearrange into cell array.
if ~iscell(range),
    for ii=1:numel(range),
        range2{ii} = range(ii);
    end
    range=range2;
end

% loop for each identifcation needed
for jj=1:numel(range)
    % convert ranges in cell array into cell arrays of strings for identifyRestriction
    for ii=1:numel(range{jj})
        fns{jj}{ii}=sprintf('%s%02.0f',prefix,range{jj}(ii));
    end
    
    % do identification
    models{jj} = method(fns{jj},margin,withFig); % use default margin
    
    % save graphs if requested
    if doPics,
        saveRestrictionFig;
        close;
    end
end

% these *should* all be the same for each test case. I haven't checked, but
% hard to envisage them changing.
fn = fieldnames(models{1}(1));
discard=[]; % we need to ignore the fields 'label' and 'c'. I should probably stop trying to be clever and just hard code fn. Or try to be cleverer about how I store the data.
for ii=1:numel(fn),
        if strcmp(fn{ii},'label') || strcmp(fn{ii},'c'),
            discard(end+1)=ii;
        end
end
fn(discard)=[];

% initalise string storage
str=[];

str=[str sprintf('Name')];
for jj = 1:numel(models{1})
    for ii = 1:numel(fn),
        str=[str sprintf('\t%s:%s',models{1}(jj).label,fn{ii})];
    end
end
str=[str sprintf('\n')];

headerstr = str;
str=[];

for jj=1:numel(models),
    % data rows
    for ii=1:numel(fns{jj}),
        str=[str sprintf('%s,',fns{jj}{ii})];
    end
    str(end)=[];
    for kk = 1:numel(models{jj}),
        for ii = 1:numel(fn),
            str=[str sprintf('\t%e',models{jj}(kk).(fn{ii}))];
        end
    end
    str=[str sprintf('\n')];
end

fprintf('\n\n');
fprintf(headerstr);
fprintf(str);
clipboard('copy', [headerstr str]);
fprintf('\nData rows copied to clipboard!\n');


function saveRestrictionFig
% saveRestrictionFig saves a nice version of the main
% p-Q plot in ./graphs/%s (%s = figure title)
ax = findobj(gcf,'tag','RestrictorFig');
fig=figure;
copyobj([ax get(ax,'Legend')],fig);
set(gca,'units','normalized','position',[0.0903    0.1769    0.8708    0.7606]);
set(gcf,'units','pixels','position',[161   296   702   281]);
tit=get(gca,'title');
name = tit.String;
name = split(name,':'); % remove description from filename
name = name{1};
%delete(tit);
saveas(fig,['./reference_material/restrictor_graphs/' name '.fig']);
saveas(fig,['./reference_material/restrictor_graphs/' name '.png']);
close(fig);

