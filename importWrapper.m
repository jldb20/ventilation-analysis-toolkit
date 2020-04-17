function [sig,par] = importWrapper(fn)
% does some pre-processing of VT+ input files
% - imports both sig and par files
% - removes invalid data (checked & identified manually)
% - organises data into convenient structures
% - checks for files from second meter and imports both datasets if they exist 
%
% sig(n) has fields: p, Q, V, t (cmH2O, L/s, L, s)
%                    pressure, flow rate, volume: raw data from Fluke VT+
%                    time added afterward using known sample rate
% par(n) has fields estimated by Fluke (taken as mean of rows in .par file):
% par(n).C compliance (L/cmH2O)
% par(n).PIP = Peak Inspiration Pressure (absolute value) (cmH2O)
% par(n).PEEP = Positive End Expiratory Pressure (cmH2O)
% par(n).TV = tidal volume (L)
% par(n).RR = respiration rate (1/s)
% par(n).inhold = inspiration hold time (s)
% par(n).exhold = expiration hold time (s)
% par(n).intime = inspiration time (s)
% par(n).extime = expiration time (s)
% par(n).IE = I:E ratio (as decimal: intime/extime)
% where n is the meter number
%

sig_in(1) = importVT(['FlowMeter1/' fn '.sig']);
par_in(1) = importVT(['FlowMeter1/' fn '.par']);

% check for second flow meter file
fn2 = ['FlowMeter2/' fn];
if exist([fn2 '.sig'],'file'),
    sig_in(2) = importVT([fn2 '.sig']);
    par_in(2) = importVT([fn2 '.par']);
end
n_meters = numel(sig_in);

% fix known problems with some of the datasets (cut data from start and end):
% Define default start and end indices as full range for each dataset, as
% well as rows in .par file to get rid of (none by default), then change
% below for specific datasets as needed.
for ii=1:n_meters      
    startindex(ii)=1;  
    endindex(ii)=numel(sig_in(ii).time);
    parcull{ii}=[];
end
switch fn
    case 'v08'
        startindex(2) = find(sig_in(2).time>13.66,1);
        parcull{2} = 1:4;
    case 'v18'
        endindex(2) = find(sig_in(2).time>49.86,1);
    case 'v28'
        startindex(2) = find(sig_in(2).time>54.6,1);
    case 'v29'
        endindex(1) = find(sig_in(1).time>24.6,1);
    case 'v32'
        startindex(1) = find(sig_in(1).time>33.7,1);
        startindex(2) = find(sig_in(2).time>33.7,1);
    case 'v116'
        endindex(2) = find(sig_in(2).time>40.2,1);
    case 'M06'
        startindex(1) = find(sig_in(1).time>75.6,1);
        parcull{1} = 1:18;
    case 'M09'
        endindex(1) = find(sig_in(1).time>31.8,1);
        parcull{1}=8:11;
    case 'M10'
        endindex(1) = find(sig_in(1).time>73.2,1);
        parcull{1} = 20:21;
    case 'M15'
        startindex(1) = find(sig_in(1).time>13.9,1);
        parcull{1} = 1:4;
    case 'M17'
        startindex(1) = find(sig_in(1).time>32.72,1);
    case 'M26'
        startindex(1) = find(sig_in(1).time>19.74,1);
        parcull{1}=1:4;
    case 'M29'
        endindex(1) = find(sig_in(1).time>39.54,1);
    case {'M34','M35','M36','M37','M38','M39','M40','M41','M42'}
        % flow data back to front in meter 2, needs inverting to match
        % meter 1.
        sig_in(2).data(:,1)=-sig_in(2).data(:,1);
        sig_in(2).data(:,2)=-sig_in(2).data(:,2);
    case {'S01','S02','S03','M46','M47','M48','M49','M50','M51','M52','M53','M54','M55','M56','M57','M58','M59','M60','M61','M62','M63','M64','M65','M66','M67'}
        error('results past M45 have not yet been checked for valid raw input data');
    otherwise
        % no action
end
for ii=1:n_meters % now crop datasets
    sig_in(ii).data(1:startindex(ii)-1,:) = [];
    sig_in(ii).time(1:startindex(ii)-1)   = [];
    sig_in(ii).data(endindex(ii)+1:end,:) = [];
    sig_in(ii).time(endindex(ii)+1:end)   = [];
    par_in(ii).data(parcull{ii},:)        = [];
    sig_in(ii).time = sig_in(ii).time - sig_in(ii).time(1);
end


% extract relevant data from relevant columns
for ii = 1:n_meters,
    d{ii} = sig_in(ii).data;
    t{ii} = sig_in(ii).time; sig(ii).t = t{ii};
    h{ii} = sig_in(ii).header;
    for jj=1:numel(h{ii}),  % sig file
        if strfind(lower(h{ii}{jj}),'flow'),
            sig(ii).Q = d{ii}(:,jj);
        elseif strfind(lower(h{ii}{jj}),'airway pressure'),
            sig(ii).p = d{ii}(:,jj); % Paw
        elseif strfind(lower(h{ii}{jj}),'volume'),
            sig(ii).V = d{ii}(:,jj);
        end
    end
    
    if size(par_in(ii).data,1)>0,
        for jj=1:numel(par_in(ii).header),  % par file
            if strfind(lower(par_in(ii).header{jj}),'compliance'),
                par(ii).C = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'pip'),
                par(ii).PIP = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'peep'),
                par(ii).PEEP = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'tidal volume'),
                par(ii).TV = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'breath rate'),
                par(ii).RR = mean(par_in(ii).data(:,jj)); % "respiration rate"
            elseif strfind(lower(par_in(ii).header{jj}),'in. hold'),
                par(ii).inhold = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'ex. hold'),
                par(ii).exhold = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'in. time'),
                par(ii).intime = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'ex. time'),
                par(ii).extime = mean(par_in(ii).data(:,jj));
            elseif strfind(lower(par_in(ii).header{jj}),'i:e ratio'),
                par(ii).IE = mean(par_in(ii).data(:,jj));
            end
        end
    else % in case .par file has no data rows
        par(ii).C = 0;
        par(ii).PIP = 0;
        par(ii).PEEP = 0;
        par(ii).TV = 0;
        par(ii).RR = 0;
        par(ii).inhold = 0;
        par(ii).exhold = 0;
        par(ii).intime = 0;
        par(ii).extime = 0;
        par(ii).IE = 0;
    end
    
    % Use consistent units
    % change lpm to L/s and mL to L and min to sec
    sig(ii).Q = sig(ii).Q/60;
    sig(ii).V = sig(ii).V/1000;
    par(ii).C = par(ii).C/1000;
    par(ii).TV = par(ii).TV/1000;
    par(ii).RR = par(ii).RR/60;
end
