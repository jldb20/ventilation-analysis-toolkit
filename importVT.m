function sig = importVT(fn)
% imports data from VT+ flow meter .par or .sig file with filename fn

% open file
fh = fopen(fn);
if fh==-1, error('Couldn''t find file %s\n',fn); end

try
    % get header row
    s = fgetl(fh);
    sig.header = split(s,char(9)); % splits at tab character
    sig.header(end) = []; % ignore time column (only timestamps start of file, or when specifically actioned on the meter - we're not using this function)
    
    % read numeric lines
    sig.data=[]; cnt = 0;
    while ~feof(fh),
        cnt=cnt+1;
        for ii=1:numel(sig.header)-1, % last column is done separately (finishes with '\n' not '\t')
            s = fscanf(fh,'%s\t',1);
            if numel(s)<2 || strcmp(s(1:2),'--'), % check for no data ('-----')
                n=nan;
            else
                if sig.header{ii}(1:3)=='I:E',
                    a=sscanf(s,'%f:%f',[2,1]);
                    n = a(1)/a(2); % I:E ratio needs converting to decimal
                else
                    n=sscanf(s,'%f'); % anything except I:E ratio
                end
            end
            sig.data(ii,cnt)=n;
        end
        s = fgetl(fh); % get the last column
        if numel(s)<2 || strcmp(s(1:2),'--'), % check for no data ('-----')
            n=nan;
        else
            n=sscanf(s,'%f',1);
        end
        sig.data(ii+1,cnt)=n;
    end
    
    % transpose into columns as in file
    sig.data = sig.data';
    
    % guess timestamps assuming 50Hz logging (from VT+ manual)
    sig.time = linspace(0,(size(sig.data,1)-1)/50,size(sig.data,1));
    sig.time = sig.time(:); % to match orientation of data matrix.
    
    % close file
    fclose(fh);
catch err
    fclose(fh);
    rethrow(err);
end
