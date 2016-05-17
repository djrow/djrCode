function paramsAn=masterfit3(mainfold,varargin)
% mainfold should have nd2 files or tiff stacks. if there are phasemask
% tiff files, they should correspond 1:1 with the nd2 files. the outputs
% written to disk are a 2Dtracks.fig, fits.fig, _analysis.mat, and a binary
% version of the nd2 (or tiff stack) movie.
%
% inside the _analysis.mat file there are the good fits, the tracked
% particles, the guessed particle positions, and a human-readable output
%
% there are various usage settings below.

% pay no attention to the man behind the curtain
%#ok<*PFBNS>
c=onCleanup(@()eval('fclose all'));
opts.WindowStyle='normal';
numDlgLines=1;

%% parameters
% CAMERA PARAMETERS
cParams.pixelSize = 49;
cParams.frameRate = 1/.04;

% ANALYSIS PARAMETERS
paramsAn.checkVals = 0;
paramsAn.phase = 1;
paramsAn.oldPhaseMasks = 1;
paramsAn.fitting = 1;
paramsAn.tracking = 0;

paramsAn.makeXls = 0;
paramsAn.bgsub = 0;
paramsAn.bgsubWidth = 50;
paramsAn.offset = 1000;
paramsAn.viewFits = 0;
paramsAn.minSep = 5;
paramsAn.widthLB = 1;
paramsAn.widthUB = 15;
paramsAn.loadOldParams = 0;
% fieldsAn = fieldnames(paramsAn);

% PHASEMASK PARAMETERS.
paramsPhase.nDilation = 2;
paramsPhase.lowThresh = 0;
paramsPhase.highThresh = 1;
paramsPhase.autofill = 1;
paramsPhase.minArea = 100;
paramsPhase.maxArea = 1e4;
fieldsPhase = fieldnames(paramsPhase);

% PEAK GUESSING PARAMETERS,
paramsPeaks.spotSizeLB = 1;
paramsPeaks.spotSizeUB = 10;
paramsPeaks.intThresh = 300;
paramsPeaks.hMax = 200;
paramsPeaks.lZero = 10;
fieldsPeaks = fieldnames(paramsPeaks);

% TRACKING PARAMETERS
paramsTr.minMerit = .08;
paramsTr.intTime = 1/cParams.frameRate;
paramsTr.gamma = 3;
paramsTr.minTrLength = 3;
paramsTr.maxStepSize = 10;
paramsTr.swh = 1;
paramsTr.delay = 0;
fieldsTr = fieldnames(paramsTr);

%% find all the movie files
[datalist,dataloc,~]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
    'Select the movies','multiselect','on');

if ~dataloc
    display('no data selected')
    return
end

if ~iscell(datalist); datalist={datalist}; end
% for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
datalist = cellfun(@(x)[dataloc, x],datalist,'uniformoutput',false);
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% initialize analysis *.mat files
if paramsAn.phase
    [phaselist,phaselistloc,~]=uigetfile([mainfold filesep...
        '*.nd2;*.tif*;*.bin'],'Select the phase images','multiselect','on');
    
    if ~phaselistloc
        display('no data selected')
        return
    end
    
    if ~iscell(phaselist); phaselist={phaselist}; end
%     for ii=1:numel(phaselist); phaselist{ii}=[phaselistloc phaselist{ii}]; end
    phaselist = cellfun(@(x)[phaselistloc, x],phaselist,'uniformoutput',false);
    [plocs,pnames,pexts]=cellfun(@fileparts,phaselist,'uniformoutput',false);
    
    for ii=1:numel(pnames)
        m=matfile([fullfile(plocs{ii},pnames{ii}),'_analysis.mat'],'Writable',true);
        
        if strcmp(pexts{ii},'.nd2')
            vidid = bfGetReader(phaselist{ii});
            phaseImg = bfGetPlane(vidid,ii);
        else
            phaseImg=imread(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
        end
        
        m.phaseMask=phaseImg;
    end    
else
    for ii=1:numel(dnames)
        [~,~,sz]=binGetFrames2([fullfile(dlocs{ii},dnames{ii}),'.bin'],1);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phaseMask=ones(sz);
        m.phaseMask=phaseMask;
    end
end

%% WRITE BINARY FILES
if paramsAn.bgsub
    bFiles=dir([dataloc,'*_bgsub.bin']);
    for ii=1:numel(dlocs);
        if ~any(ismember({bFiles.name},[dnames{ii},'_bgsub.bin']))
            
            fid=fopen([fullfile(dlocs{ii},dnames{ii}) '_bgsub.bin'],'W');
            [~,nframes,sz]=binGetFrames2([fullfile(dlocs{ii},dnames{ii}),'.bin'],1);
            
            for jj=1:floor(nframes/subwidth)
                waitbar(jj/floor(nframes/subwidth),h1)
                if jj~=floor(nframes/subwidth)
                    v=binGetFrames2([fullfile(dlocs{ii},dnames{ii}) '.bin'],...
                        1+subwidth*(jj-1):subwidth*jj);
                else
                    v=binGetFrames2([fullfile(dlocs{ii},dnames{ii}) '.bin'],...
                        1+subwidth*(jj-1):nframes);
                end
                
                fwrite(fid,uint16(bsxfun(@plus,double(v),-mean(v,3))+offset),'uint16');
            end
            fwrite(fid,sz,'double');
            fwrite(fid,nframes,'double');
            
            fclose(fid);
            
            % rename the data file names
            dnames{ii}=[dnames{ii},'_bgsub'];
        end
        
        % rename the data names variable
        dnames{ii}=[dnames{ii},'_bgsub'];
    end
elseif paramsAn.bgsub == 0
    bFiles=dir([dataloc,'*.bin']);
    for ii=1:numel(dnames);
        if ~any(ismember({bFiles.name},[dnames{ii},'.bin']))
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]));
        end
    end
end

%% compute phase mask
if paramsAn.phase&&~paramsAn.oldPhaseMasks
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        % phase image was first saved as the phase mask
        phaseImg=m.phaseMask;
        
        going = 1;
        while going
            oldParams = [paramsPhase.nDilation,paramsPhase.lowThresh,paramsPhase.highThresh,...
                paramsPhase.autofill,paramsPhase.minArea,paramsPhase.maxArea];
            
            phaseMask=valley2(phaseImg,paramsPhase);
            
            if paramsAn.checkVals
                subplot(121)
                imshow(phaseMask,[])
                subplot(122)
                imshow(phaseImg,[])
                
                % prompt user for parameter changes
                dlgPrompt=fieldsPhase;
                dlgTitle='Phasemask Parameters';
                dValues = cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
                                
                newParams=cellfun(@str2double,...
                    inputdlg(dlgPrompt,dlgTitle,numDlgLines,dValues,opts))';
                
                if isempty(newParams)
                    going = 0;
                elseif any(newParams~=oldParams)
                    paramsPhase.nDilation = newParams(1);
                    paramsPhase.lowThresh = newParams(2);
                    paramsPhase.highThresh = newParams(3);
                    paramsPhase.autofill = newParams(4);
                    paramsPhase.minArea = newParams(5);
                    paramsPhase.maxArea = newParams(6);
                else
                    going = 0;      % advance to next phase image if no parameters changed
                end
            else
                going = 0;
            end
        end
        
        phaseMask = selectCells(phaseMask,phaseImg);
        
        m.phaseMask = phaseMask;
        m.paramsPhase = paramsPhase;
    end
end

%% analyze movie files
for currMovie=1:numel(dnames)
    [~,nFrames,vSize]=binGetFrames2([fullfile(...
        dlocs{currMovie},dnames{currMovie}),'.bin'],1);
    
    m=matfile([fullfile(dlocs{currMovie},dnames{currMovie}),...
        '_analysis.mat'],'Writable',true);
    
    try
        phaseMask = m.phaseMask;
    catch
        warning('analysis has not been properly initialized. missing phase masks.')
        phaseMask = ones(vSize);
    end
    
    phaseImg=imread(fullfile(plocs{currMovie},[pnames{currMovie},pexts{currMovie}]));
    
    goodFits = zeros(0,23);
    if  ~paramsAn.fitting
        try
            goodFits = m.goodFits;
            guesses = m.guesses;
        catch
            warning(['missing data in analysis file number ' num2str(currMovie)])
            continue
        end
    else
        %% fitting
        fprintf(['Fitting this movie: ',dnames{currMovie},'\n'])
        
        p = cell(1,nFrames);
        ers = cell(1,nFrames);
        guesses = cell(1,nFrames);
        outPut = cell(1,nFrames);
        if paramsAn.checkVals
            % parameter checking version
            currFrame = 1;
            while currFrame<=nFrames
                thisframe=binGetFrames2(...
                    [fullfile(dlocs{currMovie},dnames{currMovie}) '.bin'],currFrame);
                
                [p{currFrame},ers{currFrame},guesses{currFrame},outPut{currFrame}] = ...
                    gaussFit(thisframe,phaseMask,...
                    'spotSizeLB',paramsPeaks.spotSizeLB,...
                    'spotSizeUB',paramsPeaks.spotSizeUB,...
                    'intThresh',paramsPeaks.intThresh,...
                    'hMax',paramsPeaks.hMax,...
                    'lZero',paramsPeaks.lZero,...
                    'showGuessing',paramsAn.checkVals,...
                    'frameNumber',currFrame);
                
                % prompt user for parameter changes (and advance to next frame)
                oldParams = [paramsPeaks.spotSizeLB,paramsPeaks.spotSizeUB,...
                    paramsPeaks.intThresh,paramsPeaks.hMax,paramsPeaks.lZero,currFrame+1];
                dlgPrompt=[fieldsPeaks;'next frame number'];
                dlgTitle='Peak Guessing. Close to advance to next movie.';
                dValues=cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
                
                newParams = cellfun(@str2double,...
                    inputdlg(dlgPrompt,dlgTitle,numDlgLines,dValues,opts))';
                
                if ~isempty(newParams)&&any(newParams ~= oldParams)
                    paramsPeaks.spotSizeLB = newParams(1);
                    paramsPeaks.spotSizeUB = newParams(2);
                    paramsPeaks.intThresh = newParams(3);
                    paramsPeaks.hMax = newParams(4);
                    paramsPeaks.lZero = newParams(5);
                    newParams(6) = newParams(6) - 1;    % don't advance to next frame
                end
                
                if isempty(newParams)
                    currFrame = nFrames+1;              % advance to next movie
                else
                    currFrame = newParams(6);           % advance to user input frame.
                end
            end
            
        else
            % parfor loop version
            parfor currFrame = 1:nFrames
                thisframe=binGetFrames2(...
                    [fullfile(dlocs{currMovie},dnames{currMovie}) '.bin'],currFrame);
                
                [p{currFrame},ers{currFrame},guesses{currFrame},outPut{currFrame}] = ...
                    gaussFit(thisframe,phaseMask,...
                    'spotSizeLB',paramsPeaks.spotSizeLB,...
                    'spotSizeUB',paramsPeaks.spotSizeUB,...
                    'intThresh',paramsPeaks.intThresh,...
                    'hMax',paramsPeaks.hMax,...
                    'lZero',paramsPeaks.lZero,...
                    'showGuessing',paramsAn.checkVals,...
                    'frameNumber',currFrame);
            end
            
            nFits = cellfun(@(x)size(x,1),p,'uniformoutput',0);
            
            % results selection and reorganization to legacy format
            orgFun = @(in_p,in_ers,nf,fn,pmID) cat(2,...
                fn*ones(nf,1),(1:nf)',...       % frame number, fit number
                in_p(:,1),in_ers(:,1),...       % center 1, center 1 error
                in_p(:,2),in_ers(:,2),...       % center 2, center 2 error
                in_p(:,3),in_ers(:,3),...       % s.d., s.d. error
                in_p(:,4),in_ers(:,4),...       % amplitude, amplitude error
                in_p(:,5),in_ers(:,5),...       % offset, offset error
                nan(nf,1),nan(nf,1),nan(nf,1),nan(nf,1),nan(nf,1),...
                nan(nf,1),nan(nf,1),nan(nf,1),... % 8 empty columns
                phaseMask(pmID),...             % phase mask ids
                nan(nf,1),nan(nf,1));           % empty columns
            
            % cell id assignments
            posInds = cellfun(@(x)x(:,1:2),p, 'uniformoutput', false);
            posInds = cat(1,posInds{:});
            posInds(posInds(:,1) > vSize(1) | posInds(:,2) > vSize(2) |... 
                any(posInds < 1, 2),:) = nan;
            
            allFits = cellfun(orgFun,p,ers,nFits,num2cell(1:nFrames),...
                sub2ind(vSize,posInds(:,1),posInds(:,2)), 'uniformoutput',0);
            allFits = cat(1,allFits{:});
            
            whichGood = ...
                allFits(:,7) > paramsAn.widthLB & ...
                allFits(:,7) < paramsAn.widthUB & ...
                allFits(:,21) > 0;
            
            goodFits = allFits(whichGood,:);
            
            % WRITE ANALYSIS FILE
            m.goodFits = goodFits;
            m.allFits = allFits;
            m.guesses = guesses;
            m.outPut = outPut;
            
            imshow(phaseImg,[]); hold all
            scatter(goodFits(:,5),goodFits(:,3),'.'); hold off
            saveas(gcf,[fullfile(dlocs{currMovie},dnames{currMovie}) '_goodFits.fig'])
            close
        end
    end
    
    %% Tracking
    if ~paramsAn.tracking
        try
            trfile = m.trackfile;
        catch
            trfile = zeros(0,16);
            warning(['no tracking data in analysis file number ' num2str(currMovie)])
            continue
        end
    else
        
        % if goodfits is empty, don't bother tracking.
        if isempty(goodFits)
            going = 0;
            trfile = zeros(0,16);
        else
            going = 1;
        end
        
        while going
            oldParams = [paramsTr.minMerit, paramsTr.intTime, paramsTr.gamma, ...
                paramsTr.minTrLength, paramsTr.maxStepSize, paramsTr.swh, ...
                paramsTr.delay];
            alpha=-log(oldParams(1))/oldParams(5);
            trfile=Track_3D2(goodFits,[],[],oldParams(1),...
                alpha,oldParams(3),oldParams(4),oldParams(6),...
                cParams.pixelSize,oldParams(7),oldParams(2));
            
            if isempty(trfile)
                fprintf(['No available tracks for ''',dnames{currMovie}...
                    '''. Check tracking parameters.\n']);
            else
                hastrack=unique(trfile(:,1))';
                cm=jet(numel(hastrack));
                
                imshow(phaseImg,[]); hold all
                
                for trnum = hastrack
                    plot(trfile(trfile(:,1)==trnum,5),trfile(trfile(:,1)==trnum,4),...
                        'Color',cm(hastrack==trnum,:),'linewidth',1);
                end
                title([num2str(numel(hastrack)) ' tracks'])
                hold off
                
                saveas(gcf,[fullfile(dlocs{currMovie},dnames{currMovie}),'_tracks.fig']);
            end
            
            if ~paramsAn.checkVals
                % progress to next movie if parameters are assumed correct
                going = 0;
            else
                % prompt user for parameter changes
                dlgPrompt=fieldsTr;
                dlgTitle='Tracking';
                dValues={num2str(oldParams(1)),num2str(oldParams(2)),...
                    num2str(oldParams(3)),num2str(oldParams(4)),...
                    num2str(oldParams(5)),num2str(oldParams(6)),...
                    num2str(oldParams(7))};
                
                newParams = cellfun(@str2double,inputdlg(dlgPrompt,...
                    dlgTitle,numDlgLines,dValues,opts))';
                if isempty(newParams)
                    going = 0;
                end
                
                if ~isempty(newParams)&&any(newParams~=oldParams)
                    paramsTr.minMerit = newParams(1);
                    paramsTr.intTime = newParams(2);
                    paramsTr.gamma = newParams(3);
                    paramsTr.minTrLength = newParams(4);
                    paramsTr.maxStepSize = newParams(5);
                    paramsTr.swh = newParams(6);
                    paramsTr.delay = newParams(7);
                else
                    going = 0;  % progress to next movie if parameters are left unchanged
                end
            end
            
            m.paramsTr=paramsTr;
        end
        
        % record tracking results
        m.trackfile=trfile;
    end
    
    %% misc. outputs
    % Output ViewFit Files for the Current Movie (If Selected)
    if any(ismember(paramsAn.viewFits,currMovie))||isinf(paramsAn.viewFits)
        % Output ViewFits frames if the fit file is not empty
        display(['Generating ViewFits frames for ''',...
            fullfile(dlocs{currMovie},dnames{currMovie}),'''...'])
        
        if paramsAn.phase
            bounds = [find(any(phaseMask,2),1,'first'),find(any(phaseMask,2),1,'last'),...
                find(any(phaseMask,1),1,'first'),find(any(phaseMask,1),1,'last')];
        else
            bounds = [1,size(phaseMask,1),1,size(phaseMask,2)];
        end
        
        Viewfits3([fullfile(dlocs{currMovie},dnames{currMovie}),'.bin'],...
            goodFits,guesses,trfile,cParams.frameRate,bounds);
    end
    
    if paramsAn.makeXls
        tstr=[{'Frame'} {'Molecule'} {'Amplitude'} {'+/-'} {'Offset'}...
            {'+/-'} {'X-Width'} {'+/-'} {'X Center'} {'+/-'}  {'Y Center'}...
            {'+/-'} {'Good Fit?'} {'Integral'} {'Small box width'}...
            {'Y-Width'} {'+/-'} {'Wx/Wy'} {'Z Center'} {'+/-'} {'Cell No'}...
            {'sROI'} {'tROI'}];% {'distance from center'}];
        xlswrite([dnames{currMovie}],cat(1,tstr,num2cell(goodFits)))
    end
end
end

function Viewfits3(vidloc,fitfile1,guesses,trfile,framerate,bounds)
% This code takes in raw single-molecule tif frames and fits files and put
% a square at each fit. Useful for checking to see if fitting parameters
% are right.

% INPUTS:

% OUTPUTS:

halfBoxWidth = 7;
intMax=1;

[~,nframes,~]=binGetFrames2(vidloc,1);
sz = [bounds(2)-bounds(1)+1,bounds(4)-bounds(3)+1];

wobj=VideoWriter([vidloc(1:end-4),'_Viewfits'],'Uncompressed AVI');
wobj.FrameRate=framerate;

open(wobj)
d=onCleanup(@()close(wobj));

h=waitbar(0,['writing viewfits for movie [...' vidloc(end-20:end) ']']);
c=onCleanup(@()close(h));

for ii=1:nframes
    if rem(ii,10)==0
        waitbar(ii/nframes,h)
    end
    
    img=double(binGetFrames2(vidloc,ii));
    img=img(bounds(1):bounds(2),bounds(3):bounds(4));
    
    % autoscale
    img=img-min(img(:));
    img=img/max(img(:));
    
    % Loop for all guess peaks in the current frame
    for jj = 1:size(guesses{ii},1)
        pix_x=guesses{ii}(jj,2)-bounds(3);
        pix_y=guesses{ii}(jj,1)-bounds(1);
        if pix_x<sz(2)&&pix_y<sz(1)&&pix_x>0&&pix_y>0
            img=makeBox(img,intMax/2,pix_x,pix_y,halfBoxWidth,sz);
        end
    end
    
    if any(size(img)~=sz)
        keyboard
    end
    
    % Indexing which rows in the fits file belong to the current frame
    r=find(fitfile1(:,1)==ii);
    
    % Loop for all good fits in the current frame
    for jj=1:numel(r)
        pix_x=round(fitfile1(r(jj),5))-bounds(3);
        pix_y=round(fitfile1(r(jj),3))-bounds(1);
        
        img=makeBox(img,intMax,pix_x,pix_y,halfBoxWidth,sz);
    end
    
    if any(size(img)~=sz)
        keyboard
    end
    
    if numel(trfile)>0
        tr=trfile(trfile(:,2)==ii,:);
        
        cm=jet(5);
        img=img(:,:,[1,1,1]);
        for jj=1:size(tr,1)
            pix_x=round(tr(jj,5))-bounds(3);
            pix_y=round(tr(jj,4))-bounds(1);
            img=makeBox(img,cm(mod(tr(jj,1),size(cm,1)-1)+1,:),...
                pix_x,pix_y,halfBoxWidth,sz);
        end
    end
    
    if any(size(img(:,:,1))~=sz)
        keyboard
    end
    
    %     imwrite(currentframe2, [vidloc(1:end-4),'_Viewfits.tif'], 'writemode', 'append');
    writeVideo(wobj,img)
end
close(wobj);
end

function currentframe=makeBox(currentframe,curr_fr_maxint,...
    pix_x,pix_y,half_symbol_size,sz)

for ii=1:size(currentframe,3)
    b=pix_y-half_symbol_size;
    t=pix_y+half_symbol_size;
    l=pix_x-half_symbol_size;
    r=pix_x+half_symbol_size;
    
    b(b<1)=1;
    l(l<1)=1;
    t(t>sz(1))=sz(1);
    r(r>sz(2))=sz(2);
    
    currentframe(b:t,[l,r],ii)=curr_fr_maxint(ii);
    currentframe([b,t],l:r,ii)=curr_fr_maxint(ii);
end
end

function phaseMask=selectCells(phaseMask,img)
% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

% INPUTS:

% phaseMask: Numbered phase mask of cell images returned by the 'valley.m'
% code or any other cell segmentation code

% OUTPUTS:

% phaseMask

%#ok<*AGROW>

h=figure;
c=onCleanup(@()close(h));
vSize = size(phaseMask);

subplot(121);
imshow(img,[],'Border','Tight');
subplot(122);
imshow(phaseMask~=0,[],'Border','Tight');
title(sprintf(['Click on the cells to be analyzed.\n Press ''Enter'' to'...
    ' proceed.\n Or hit enter without clicking\n on any cell to analyze ALL of them.']))

% Now move the axes slightly so that the top of the title is visible
set(h,'Units','normalized')
P=get(h,'Position');
set(h,'Position',P+[0 -0.03 0 0])
hold all

% highlight the outlines of the cells
rProp=regionprops(phaseMask,'Convexhull');
for i=1:length(rProp)
    plot(rProp(i,1).ConvexHull(:,1),rProp(i,1).ConvexHull(:,2),'c-','linewidth',2)
end

% Let the user pick cells with good shapes
keyInput=0; clicksX=[]; clicksY=[];
while keyInput~=121
    [clickX,clickY,keyInput]=ginput(1);
    plot(clickX,clickY,'m*','markersize',12);
    clicksX=cat(1,clicksX,clickX); 
    clicksY=cat(1,clicksY,clickY);
end

% If the user did click on something
if ~isempty(clicksX)
    posInds = ceil(cat(2,clicksX,clicksY));
    posInds(clicksX > vSize(1) | clicksY > vSize(2) | any(posInds < 1, 2),:) = nan;
            
    goodCells=phaseMask(sub2ind(vSize,posInds(:,1),posInds(:,2)));
    goodCells(goodCells == 0) = [];
    
    phaseMask(~ismember(phaseMask,goodCells))=0;
end
end