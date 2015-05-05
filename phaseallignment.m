function phaseallignment(mainfold)

if ~exist('mainfold','var')
    mainfold=pwd;
end

display('Select the movies.')
[datalist,dataloc,findex]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
    'Matlab Files','multiselect','on');

if findex==0
    fprintf('no data selected\n')
    return
end

% make sure the variable class is correct
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);

% Write bin files for all movies
for ii=1:numel(dnames);
    
    % if selected files are not bin files
    if strcmp(dexts{ii},'.nd2')||strcmp(dexts{ii},'.tif')||...
            strcmp(dexts{ii},'.tiff')
        writebin(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]));
    end
end

display('Select the phase images.')
[phaselist,phaselistloc,findex]=uigetfile([mainfold filesep...
    '*.nd2;*.tif*;*.bin'],'Matlab Files','multiselect','on');
if findex==0
    fprintf('no phase images selected. try again.\n')
    return
end

if ~iscell(phaselist); phaselist={phaselist}; end
for ii=1:numel(phaselist); phaselist{ii}=[phaselistloc phaselist{ii}]; end
[plocs,pnames,pexts]=cellfun(@fileparts,phaselist,'uniformoutput',false);

% loop through data/phase image pairs
for ii=1:numel(dnames)
    
    if strcmp(pexts{ii},'.nd2')||strcmp(pexts{ii},'.tif')||...
            strcmp(pexts{ii},'.tiff')
        % convert to binary if selected version is not binary
        writebin(phaselist{ii});
        img=bingetframes(fullfile(plocs{ii},pnames{ii},'.bin'));
    else
        % otherwise, just load the binary version
        img=bingetframes(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
    end
    
    % get the length and size of the movie
    [~,n,s]=bingetframes(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]),1);
    
    % show the phase image
    subplot(1,2,1)
    imshow(img,[])

    % continuous average of the movie
    datimg=zeros(s);
    subplot(1,2,2)
    for jj=1:n
        datimg=datimg+double(bingetframes(fullfile(dlocs{ii},...
            [dnames{ii},dexts{ii}]),jj));
        
        if rem(jj,50)==0
            % show the continuously averaging movie
            imshow(datimg,[])
            title(sprintf(['average of frames\n 1:' num2str(jj)...
                ' of ' num2str(n) '.']))
            drawnow
        end
    end
end
end