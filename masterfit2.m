function goodfitdata=masterfit2(mainfold,yesCheckVals,yesFitting)
% mainfold should have nd2 files or tiff stacks. if there are phasemask
% tiff files, they should correspond 1:1 with the nd2 files. the outputs
% written to disk are a 2Dtracks.fig, fits.fig, _analysis.mat, and a binary
% version of the nd2 (or tiff stack) movie.
%
% inside the _analysis.mat file, there are the good fits, the tracked
% particles, the guessed particle positions, and a human-readable output
%
% there are various usage settings below.

%% CAMERA PARAMETERS
% Pixel size in nm.
pxSize=49;

% camera capture rate in 1/s
frameRate=25;


%% PHASEMASK PARAMETERS. 
% dilate factor, low thresh, high thres, autofill, min area, max area
phaseparams=[1,0,1,0,100,10000];


%% PEAK GUESSING PARAMETERS, 
% noise size, particle size, Intensity Threshold, H-Max, lzero]
peak_guessing_params=[1,10,20,20,10];


%% TRACKING PARAMETERS
% minimum merit
trackparams(1)=0.1;

% Integration time (ms)
trackparams(2)=1/frameRate;

% gamma
trackparams(3)=3;

% maximum step size
trackparams(4)=8;

% minimum track length
trackparams(5)=3;

% speed estimation window halfsize
trackparams(6)=1;

% time delay between consecutive frames (ms)
trackparams(7)=0;


%% ANALYSIS PARAMETERS
% run tracking?
yesTracking=1;

% use phasemasks?
yesPhasemasks=0;

% if you're happy with the phasemask inside the analysis files, don't make
% another one
happyWithPhasemask=1;

% skip manual selection of actual bacterial cells. do you trust that the
% mask finding algorithm won't make mistakes?
skipSelect=0;

% plot things?
yesPlot=1;

% run 3d code?
yes3d=0;

% make pov-ray csv files?
yespovray=0;

% make xls file from goodfitdata
makeXls=0;

% remove activation?
removeactivation=0;

% background subtract? subwidth must be an integer.
yesbgroundsub=1;
subwidth=50;
offset=1000;

% Input which movie(s) will generate the ViewFits frames that are useful
% for checking guessing/fitting parameters. Use "0" if you do not want any
% movie to output the ViewFits frames. Use "inf" if you want all movies to
% generate ViewFits frames. Example, 'viewfits = [1 3]' will tell the code
% to output ViewFits frames for the 1st and the 3rd movie.
viewfits_mv=[1];

% Minimal separation of peaks (px). Of two putative peaks that are closer than
% this value, the dimmer of the two will be discarded.
min_sep=5;

% Lower and upper bounds (in pixels, except for the aspect ratio) for
% various fit parameters that will be used to determine whether a fit is
% "good".
width_lb=1;
width_ub=15;
width_error_ub=5;

% The maximum allowable separation (nm) of an input width pair and the 2
% defocusing calibration curves. Fits falling too far from the calibration
% curve will be rejected, have NaN for the z-center value and will not show
% up in the good fits .dat file. 
max_allowed_D=100;

% Correction factor for index of refraction mismatch. Usually less than 1.
indRefr_corr=0.79;

%% ------------------------------------------------------------------------
%  Load 3D Fitting Parameters and Phase Masks

if yes3d
    % Load defocusing_param, which stores results from fitting to the
    % calibration curves.
    m=matfile(fullfile(mainfold,'calibdata.mat'));
    if ~isempty(who(m))
        defocusing_param=m.defocusing_param;
        z_std_LUT=m.zuncLUT;
    else
        display('put a calibdata.mat file in the data folder, or change yes3d to 0')
        return
    end
end

% find all the relevant files
if yesFitting&&ischar(mainfold)         % fitting and folder location string given as input
    display('Select the movies.')
    [datalist,dataloc,findex]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
        'Matlab Files','multiselect','on');
    
    if findex==0
        fprintf('no data selected\n')
        return
    end
    
    if ~iscell(datalist); datalist={datalist}; end
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
    [dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);
    
    % WRITE BIN FILES FOR ALL MOVIES
    bFiles=dir([dataloc,'*.bin']);
    for ii=1:numel(dnames);
        if ~any(ismember({bFiles.name},[dnames{ii},'.bin']))
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]));
        end
    end
elseif ~yesFitting&&ischar(mainfold)    % not fitting and folder location string given as input
    display('Select the analysis files.')
    [datalist,dataloc,findex]=uigetfile([mainfold filesep '*.mat'],...
        'Matlab Files','multiselect','on');
    
    if findex==0
        fprintf('no data selected\n')
        return
    end
    
    if ~iscell(datalist); datalist={datalist}; end
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
    [dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);
    for ii=1:numel(datalist); dnames{ii}(end-8:end)=[]; end
elseif isstruct(mainfold)               % folder location string is actually a stucture
    [dlocs,dnames,~]=fileparts(mainfold.name);
    dlocs={dlocs};
    dnames={dnames};
    dataloc=dlocs{1};
end

%% WRITE BACKGROUND SUBTRACTED BINARY FILE
if yesbgroundsub&&yesFitting
    bFiles=dir([dataloc,'*_bgsub.bin']);
    for ii=1:numel(dlocs);
        if ~any(ismember({bFiles.name},[dnames{ii},'_bgsub.bin']))&&...
                ~any(ismember({bFiles.name},dnames{ii}))
            h1=waitbar(0);
            
            fid=fopen([fullfile(dlocs{ii},dnames{ii}) '_bgsub.bin'],'W');
            [~,nframes,sz]=bingetframes([fullfile(dlocs{ii},dnames{ii}),'.bin'],1,[]);
            
            waitbar(0,h1,['writing bg sub frames for ' datalist{ii}(1:end-4)]);
            for jj=1:floor(nframes/subwidth)
                waitbar(jj/floor(nframes/subwidth),h1)
                if jj~=floor(nframes/subwidth)
                    v=bingetframes([fullfile(dlocs{ii},dnames{ii}) '.bin'],...
                        1+subwidth*(jj-1):subwidth*jj,[]);
                else
                    v=bingetframes([fullfile(dlocs{ii},dnames{ii}) '.bin'],...
                        1+subwidth*(jj-1):nframes,[]);
                end
                
                fwrite(fid,uint16(bsxfun(@plus,double(v),-mean(v,3))+offset),'uint16');
            end
            fwrite(fid,sz,'double');
            fwrite(fid,nframes,'double');
            
            fclose(fid);
            
            % rename the data file names
            dnames{ii}=[dnames{ii},'_bgsub'];
            end
        else
            % rename the data file names
            dnames{ii}=[dnames{ii},'_bgsub'];
        end
    try
        close(h1)
    end
end

%% WRITE PHASEMASKS FILE FOR ALL MOVIES

% find phase contrast images
if yesPhasemasks&&ischar(mainfold)
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
end

if yesPhasemasks&&yesFitting&&yesCheckVals&&ischar(mainfold)
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
                
        mnamelist=who(m); mnamelist=mnamelist(:);
        if any(ismember(mnamelist,'phasemask'))&&happyWithPhasemask
            continue
        end
        if any(ismember(mnamelist,'phaseparams'))
            phaseparams=m.phaseparams;
        end
        
        if strcmp(pexts{ii},'.nd2')
            vidid=bfGetReader(phaselist{ii});
            img=bfGetPlane(vidid,ii);
        else
            img=imread(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
        end
        
        fprintf(['phasing file named: ' dnames{ii} '.\n'])
        
        counter=0;
        while counter<1
            phasemask=valley(img(:,:,1),phaseparams,yesCheckVals);
            
            % prompt user for parameter changes
            dlgPrompt={'dilation factor', 'lower threshold',...
                'higher threshold','autofill y/n','minimum cell area',...
                'maximum cell area'};
            dlgTitle='Phasemask';
            numDlgLines=1;
            def={num2str(phaseparams(1)),num2str(phaseparams(2)),...
                num2str(phaseparams(3)),num2str(phaseparams(4)),...
                num2str(phaseparams(5)),num2str(phaseparams(6))};
            opts.WindowStyle='normal';
            
            phaseparams_new=cellfun(@str2double,...
                inputdlg(dlgPrompt,dlgTitle,numDlgLines,def,opts))';
            
            if any(phaseparams_new~=phaseparams)
                phaseparams=phaseparams_new;
            else
                counter=counter+1;
            end
        end
        if ~skipSelect
            [~,goodcells]=select_cells(phasemask,img);
            phasemask(logical(-1*(ismember(phasemask,goodcells)-1)))=0;
        end
        m.phasemask=phasemask;
        m.phaseparams=phaseparams;
        m.peak_guessing_params=peak_guessing_params;
    end
elseif ~yesPhasemasks
    for ii=1:numel(dnames)
        [~,~,sz]=bingetframes([fullfile(dlocs{ii},dnames{ii}),'.bin'],1,[]);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phasemask=ones(sz);
        m.phasemask=phasemask;
        m.peak_guessing_params=peak_guessing_params;
    end
end

% Spatial ROI and temporal ROI. Though currently not used and function
% merely as place holders, these parameters will be employed to expand the
% functionality of the code.
sROI=1;
tROI=1;

% Loop each movie for guessing/fitting/tracking
for curr_mainfold=1:numel(dnames)
    [~,num_files,~]=bingetframes([fullfile(...
        dlocs{curr_mainfold},dnames{curr_mainfold}),'.bin'],1,[]);
    
    m=matfile([fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}),...
        '_analysis.mat'],'Writable',true);
    mnamelist=who(m); mnamelist=mnamelist(:);
    
    phasemask=m.phasemask;
    
    %% psf fitting
    goodfitdata={[]}; guesses=cell(1,num_files);
    
    skippedframes={};
    if ~any(ismember(mnamelist,'goodfitdata'))||yesFitting==1
        % Display movie folder counter
        fprintf(['Fitting this movie: ',dnames{curr_mainfold},'\n'])
        
        if yesCheckVals==1
            frameskip=0;
            
            yn=input(['input new peak guessing parameters? enter for ''no'', '...
                'the parameter vector for ''yes''.\n']);
            
            if ~isempty(yn)
                peak_guessing_params=yn;
                m.peak_guessing_params=peak_guessing_params;
            end
            
            for ii=1:num_files; % Loop through each frame of the current movie
                
                if frameskip<=ii
                    % Read the current frame
                    thisframe=bingetframes(...
                        [fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}) '.bin'],ii,[]);
                    if removeactivation
                        if sum(sum(thisframe))>tval
                            skippedframes{ii}=1;
                            continue
                        end
                    else
                        skippedframes{ii}=nan;
                    end
                    
                    % 3D Peak Fitting and z-position extraction -----------
                    % Fit putative peaks. peak guessing is in
                    % fit_asym_gauss_int2
                    [all_fitparam,all_fiterr,guesses{ii},frameskip]=...
                        fit_asym_gauss_int2(thisframe,phasemask,1,...
                        peak_guessing_params,min_sep,frameskip,ii);
                    
                    Wx=all_fitparam(:,3)*pxSize; Wy=all_fitparam(:,6)*pxSize;
                    
                    if yes3d==1
                        % Find the z-location (the output is converted to
                        % pixels here).
                        z_nm=find_astigZ_position2(Wx,Wy,defocusing_param,...
                            max_allowed_D,indRefr_corr);
                        z_px=z_nm/pxSize;
                        
                        % Find the localization uncertainty associated with
                        % each z-position using a lookup table generated
                        % from calibration. z_uncertainty_nm =
                        % find_z_uncertainty_from_LUT(z_nm', z_std_LUT);
                        % interp1 replaces the Yi-code
                        % find_z_uncertainty_from_LUT
                        z_uncertainty_nm=interp1(z_std_LUT(1,:),z_std_LUT(2,:),z_nm);
                        z_uncertainty_px=z_uncertainty_nm/pxSize;
                    else
                        z_px=zeros(size(all_fitparam,1),1);
                        z_uncertainty_px=zeros(size(all_fitparam,1),1);
                    end
                    
                    % Organize raw fit data
                    % ---------------------------------------
                    
                    % Get rid of fits (along with their z-positions) which
                    % have their x- and y-positions lying outside the image
                    valid_param_rows=find(all_fitparam(:,4)>=0&all_fitparam(:,5)>=0&...
                        all_fitparam(:,4)<=size(phasemask,2)&...
                        all_fitparam(:,5)<=size(phasemask,1));
                    
                    all_fitparam=all_fitparam(valid_param_rows,:);
                    all_fiterr=all_fiterr(valid_param_rows,:);
                    
                    % # of total fits from current frame (good or bad)
                    numfits=size(all_fitparam,1);
                    
                    if numfits~=0
                        raw_fitdata=[repmat(ii,[numfits,1]),(1:numfits)',... % col 1,2
                            all_fitparam(:,1),all_fiterr(:,1),all_fitparam(:,2),... % 3,4,5
                            all_fiterr(:,2),all_fitparam(:,3),all_fiterr(:,3),... % 6,7,8
                            all_fitparam(:,4),all_fiterr(:,4),all_fitparam(:,5),... % 9,10,11
                            all_fiterr(:,5),nan(numfits,1),all_fitparam(:,7),... % 12,13,14
                            repmat(min_sep*2+1,[numfits, 1]),all_fitparam(:,6),... % 15,16
                            all_fiterr(:,6),all_fitparam(:,3)./all_fitparam(:,6),... % 17,18
                            z_px(valid_param_rows),z_uncertainty_px(valid_param_rows),... % 19,20
                            phasemask(sub2ind(size(phasemask),... % 21
                            ceil(all_fitparam(:,5)),ceil(all_fitparam(:,4)))),... % 21 cont.
                            repmat(sROI,[numfits,1]),repmat(tROI,[numfits,1])];  % 22, 23
                        % Organize fitting data. Here, each column stores
                        % the following information:  (1) frame # (2)
                        % molecule # (3) amplitude (4) amplitude error (5)
                        % offset (6) offset error (7) x-width (8) x-width
                        % error (9) x-center (10) x-center error (11)
                        % y-center (12) y-center error (13) good fit? (NaN,
                        % TBD) (14) integral (15) small box width (16)
                        % y-width (17) y-width error (18) aspect ratio (x
                        % width/y width) (19) z-center (20) z-center error
                        % uncertainty (21) cell # (22) sROI (23) tROI
                        
                        raw_fitdata(:,13)=(all_fitparam(:,1)>0).*... % Amplitude > amp error
                            (all_fitparam(:,3)>all_fiterr(:,3)).*... % width > error
                            (all_fiterr(:,3)<width_error_ub).*...
                            (all_fiterr(:,6)<width_error_ub).*... % Width error < UB
                            (all_fitparam(:,3)>width_lb).*...
                            (all_fitparam(:,3)<width_ub).*... % LB < x-width < UB
                            (all_fitparam(:,6)>width_lb).*...
                            (all_fitparam(:,6)<width_ub).*... % LB < y-width < UB
                            raw_fitdata(:,21)>0.*... % Within cell boundaries
                            ~isnan(z_px(valid_param_rows));  % Non-NaN z-position
                        
                        goodfitdata{ii}=raw_fitdata(raw_fitdata(:,13)==1,:);
                    end
                end
            end
        elseif ~yesCheckVals&&yesFitting
            % parallel processing. cannot be used for debugging or changing
            % guessing/fitting parameters. this is identical to the for
            % loop above, except there're no comments to save space
            if ~yes3d
                defocusing_param=[]; % because parfor is weird
                z_std_LUT=[];
            end
            
            peak_guessing_params=m.peak_guessing_params;
            
            parfor ii=1:num_files
                [thisframe,~,~,tval]=bingetframes(...
                    [fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}) '.bin'],ii,[]);
                if removeactivation==1
                    if sum(sum(thisframe))>tval
                        skippedframes{ii}=1;
                        continue
                    end
                else
                    skippedframes{ii}=nan;
                end
                [all_fitparam,all_fiterr,guesses{ii}]=fit_asym_gauss_int2(...
                    thisframe,phasemask,1,peak_guessing_params,min_sep);
                Wx=all_fitparam(:,3)*pxSize; Wy=all_fitparam(:,6)*pxSize;
                if yes3d
                    % Find the z-location (the output is converted to
                    % pixels here).
                    z_nm=find_astigZ_position2(Wx,Wy,defocusing_param,...
                        max_allowed_D,indRefr_corr);
                    z_px=z_nm/pxSize;
                    z_uncertainty_nm=interp1(z_std_LUT(1,:),z_std_LUT(2,:),z_nm);
                    z_uncertainty_px=z_uncertainty_nm/pxSize;
                else
                    z_px=zeros(size(all_fitparam,1),1);
                    z_uncertainty_px=zeros(size(all_fitparam,1),1);
                end
                % the centers have to be inside the frame
                valid_param_rows=find(all_fitparam(:,4)>=0&all_fitparam(:,5)>=0&...
                    all_fitparam(:,4)<=size(phasemask,2)&...
                    all_fitparam(:,5)<=size(phasemask,1));
                all_fitparam=all_fitparam(valid_param_rows,:);
                all_fiterr=all_fiterr(valid_param_rows,:);
                numfits=size(all_fitparam,1);
                
                if numfits
                    raw_fitdata=[repmat(ii,[numfits,1]),(1:numfits)',... % col 1,2
                        all_fitparam(:,1),all_fiterr(:,1),all_fitparam(:,2),... % 3,4,5
                        all_fiterr(:,2),all_fitparam(:,3),all_fiterr(:,3),... % 6,7,8
                        all_fitparam(:,4),all_fiterr(:,4),all_fitparam(:,5),... % 9,10,11
                        all_fiterr(:,5),nan(numfits,1),all_fitparam(:,7),... % 12,13,14
                        repmat(min_sep*2+1,[numfits, 1]),all_fitparam(:,6),... % 15,16
                        all_fiterr(:,6),all_fitparam(:,3)./all_fitparam(:,6),... % 17,18
                        z_px(valid_param_rows),z_uncertainty_px(valid_param_rows),... % 19,20
                        phasemask(sub2ind(size(phasemask),... % 21
                        ceil(all_fitparam(:,5)),ceil(all_fitparam(:,4)))),... % 21 cont.
                        repmat(sROI,[numfits,1]),repmat(tROI,[numfits,1])];  % 22, 23
                    raw_fitdata(:,13)=(all_fitparam(:,1)>0).*... % Amplitude > amp error
                        (all_fitparam(:,3)>all_fiterr(:,3)).*... % width > error
                        (all_fiterr(:,3)<width_error_ub).*...
                        (all_fiterr(:,6)<width_error_ub).*... % Width error < UB
                        (all_fitparam(:,3)>width_lb).*...
                        (all_fitparam(:,3)<width_ub).*... % LB < x-width < UB
                        (all_fitparam(:,6)>width_lb).*...
                        (all_fitparam(:,6)<width_ub).*... % LB < y-width < UB
                        raw_fitdata(:,21)>0.*... % Within cell boundaries
                        ~isnan(z_px(valid_param_rows));  % Non-NaN z-position
                    
                    goodfitdata{ii}=raw_fitdata(raw_fitdata(:,13)==1,:);
                end
            end
        end
        
        goodfitdata=cat(1,goodfitdata{:});
        
        % WRITE FITS TO ANALYSIS FILE
        m.goodfitdata=goodfitdata;
        m.guesses=guesses;
        m.skippedframes=find(cat(1,skippedframes{:}));
    else
        goodfitdata=m.goodfitdata;
        guesses=m.guesses;
    end
    
    if isempty(goodfitdata)
        display(['the movie named ', dnames{curr_mainfold},...
            ' has yielded zero good fits'])
        continue
    end
    
    %% observe the fitting results
    if yespovray
        % WRITE CSV FILES The 2 files required to run POV-Ray reconstuction
        % is stored under the same folder as the fitting file from which
        % they are generated.
        GoodFitsFile_DatToCsv_PovRay3D2(...
            fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}),goodfitdata);
    end
    
    if yesPlot
        % WRITE FITS FIG Plot 3D Localization (Without Gaussian-Blur)
        Plot_3D_fits2(goodfitdata,3,49,...
            [fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}),'_fits.fig']);
    end
    %% ------------------------------------------------------------------------
    %  3D Tracking
    %  ------------------------------------------------------------------------
    
    % if tracking is disabled
    if any(ismember(mnamelist,'trackfile'))&&~yesTracking
        trfile=m.trackfile;
        
    % if tracking is enforced
    elseif yesTracking
        
        if yesPhasemasks            % if phasemasks are used
            % plot 2D tracks
            img=imread(fullfile(plocs{curr_mainfold},[pnames{curr_mainfold},...
                pexts{curr_mainfold}]),'tif');
        else                        % if phasemasks aren't used
            img=m.phasemask;
        end
        img=double(img);            % make sure the image isn't integers
        
        % magnify the movie?
        magfactor=1; %ceil(2e3/min(size(img)));
        img=kron(img,ones(magfactor));
        
        % TRACKING
        if yesPlot
            fprintf(['tracking file named: ' dnames{curr_mainfold} '.\n'])
        end
        
        counter=0; yn=[];
        while counter<2         % two sequential agreements finish this loop
            
            % if trackparams exist in the analysis file and it's the first
            % time the while loop has run and the user hasn't input
            % tracking parameters (at the end of this while loop)
            mnamelist=who(m); mnamelist=mnamelist(:);
            if any(ismember(mnamelist,'trackparams'))&&counter==0
                temptrackparams=m.trackparams;
                
                % if yescheckvals is enabled and any of the trackparams
                % saved in the analysis file differ from those at the top
                % of this function.
                if yesCheckVals&&any(temptrackparams~=trackparams)
                    useold=input('use tracking params from analysis file?');
                    if ~useold
                        trackparams=temptrackparams;
                    end
                end
                yn=trackparams;
            end
            
            % if (at the end of this while loop) the user input tracking
            % parameters
            if ~isempty(yn)&&any(yn~=trackparams)
                trackparams=yn;
            else
                counter=counter+1;
            end
            
            alpha=-log(trackparams(1))/trackparams(4);
            trfile=Track_3D2(goodfitdata,sROI,tROI,trackparams(1),...
                alpha,trackparams(3),trackparams(5),trackparams(6),...
                pxSize,trackparams(7),trackparams(2));
            
            if isempty(trfile)
                fprintf(['No available tracks for ''',dnames{curr_mainfold}...
                    '''. Check tracking parameters.\n']);
            elseif yesPlot
                trfile=trfile*magfactor;
                
                hastrack=unique(trfile(:,1))';
                cm=jet(numel(hastrack));
                
                pcolor(img); colormap('gray'); hold all
                shading flat; axis image
                
                % Loop through each track
                for trnum=hastrack % Loop through each track
                    tr_start_row=find(trfile(:,1)==trnum,1);
                    tr_end_row=find(trfile(:,1)==trnum,1,'last');
                    
                    % If the current track meets the length requirement.
                    if (tr_end_row-tr_start_row+1)>=trackparams(5)
                        plot(trfile(tr_start_row:tr_end_row,4),...
                            trfile(tr_start_row:tr_end_row,5),...
                            'Color',cm(hastrack==trnum,:),...
                            'linewidth',1);
                    end
                end
                title([num2str(numel(hastrack)) ' tracks'])
                hold off
                
                % WRITE 2D TRACKING FIG
                saveas(gcf,[fullfile(dlocs{curr_mainfold},...
                    dnames{curr_mainfold}),'_2Dtracks.fig']);
                
                if yes3d
                    % Plotting 3D tracks
                    
                    % Loop through each track
                    for trnum=max(trfile(:,1)):-1:1
                        tr_start_row=find(trfile(:,1)==trnum,1);
                        tr_end_row=find(trfile(:,1)==trnum,1,'last');
                        
                        % If the current track meets the length
                        % requirement.
                        if (tr_end_row-tr_start_row+1)>=min_tr_length
                            plot3(trfile(tr_start_row:tr_end_row,4)*pxSize,...
                                trfile(tr_start_row:tr_end_row,5)*pxSize,...
                                trfile(tr_start_row:tr_end_row,15)*pxSize,...
                                'Color','r');
                            hold all
                        end
                    end
                    grid on
                    % WRITE 3D TRACKING FIG
                    saveas(gcf,[datalist{curr_mainfold}(1:end-4),'_3Dtracks.fig']);
                    close
                end
                trfile=trfile/magfactor;
                
            end
            
            % end while loop if not checking values
            if ~yesCheckVals
                counter=2;
            else
                % prompt user for parameter changes
                dlgPrompt={'minimum merit', 'integration time',...
                    'gamma','maximum step size','minimum track length',...
                    'half size of speed boxcar','time delay'};
                dlgTitle='Tracking';
                numDlgLines=[1,length(dlgTitle)+10];
                def={num2str(trackparams(1)),num2str(trackparams(2)),...
                    num2str(trackparams(3)),num2str(trackparams(4)),...
                    num2str(trackparams(5)),num2str(trackparams(6)),...
                    num2str(trackparams(7))};
                opts.WindowStyle='normal';
                
                yn=cellfun(@str2double,inputdlg(dlgPrompt,dlgTitle,numDlgLines,def,opts))';
            end
        end
        m.trackfile=trfile;
        m.trackparams=trackparams;
    end
    
    % Output ViewFit Files for the Current Movie (If Selected)
    if any(ismember(viewfits_mv,curr_mainfold))||any(isinf(viewfits_mv))
        % Output ViewFits frames if the fit file is not empty
        if numel(goodfitdata)>0
            fprintf(['\nGenerating ViewFits frames for ''',...
                fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}),'''...\n'])
            
            % WRITE MP4 MOVIE OF FITS
            if yesbgroundsub&&~viewfits_onBGSUB
                Viewfits3([fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}(1:(end-6))),'.bin'],...
                    goodfitdata,guesses,trfile,7,frameRate);
            else
                Viewfits3([fullfile(dlocs{curr_mainfold},dnames{curr_mainfold}),'.bin'],...
                    goodfitdata,guesses,trfile,7,frameRate);
            end
        else
            fprintf(['The fit file for ''', dnames{curr_mainfold},...
                ''' does not exist or contains no data. Skip producing ',...
                'ViewFits files. \n']);
        end
    end
    
    if makeXls
        tstr=[{'Frame'} {'Molecule'} {'Amplitude'} {'+/-'} {'Offset'}...
            {'+/-'} {'X-Width'} {'+/-'} {'X Center'} {'+/-'}  {'Y Center'}...
            {'+/-'} {'Good Fit?'} {'Integral'} {'Small box width'}...
            {'Y-Width'} {'+/-'} {'Wx/Wy'} {'Z Center'} {'+/-'} {'Cell No'}...
            {'sROI'} {'tROI'}];% {'distance from center'}];
        xlswrite([dnames{curr_mainfold}],cat(1,tstr,num2cell(goodfitdata)))
    end
end
end
function Viewfits3(vidloc,fitfile1,guesses,trfile,half_symbol_size,framerate)

% This code takes in raw single-molecule tif frames and fits files and put
% a square at each fit. Useful for checking to see if fitting parameters
% are right.

% INPUTS:

% tifdirectory: Full path of the folder storing raw tif frames.

% goodfitsfile_full_path: Full path of the .dat file that stores fitting
% parameters of good fits. The 9-th and the 11-th columns should
% respectively store the x- and the y-location of each fit.

% viewfits_output_folder_path: Full path of the folder in which ViewFits
% frames are going to be stored.

% half_symbol_size: The half size (in pixels) of the symbol highlighting
% the fits. E.g., if half_symbol_size = 7, and a square is drawn, the
% square will be 15-by-15 pixels large.

% OUTPUTS:

% ViewFits frames are stored in the designated folder, with '_ViewFits'
% appended to the original frame name.

%-------------------------------------------------------------------------%
[~,nframes,sz]=bingetframes(vidloc,1,[]);
guesses=cellfun(@round,guesses,'uniformoutput',0);

wobj=VideoWriter([vidloc(1:end-4),'_Viewfits'],'MPEG-4');
wobj.FrameRate=framerate;
wobj.Quality=95;
% wobj.CompressionRatio=20; % for 'Motion JPEG 2000' codec

magfactor=ceil(150/min(sz)); % movie should have 150 pixels across

open(wobj)
h1=waitbar(0,'writing viewfits movie.');
for ii=1:nframes % Looping through each frame
    waitbar(ii/nframes,h1);
    currentframe=double(bingetframes(vidloc,ii,[]));
    
    % AUTOSCALING
    currentframe=currentframe-min(currentframe(:));
    
    % written video must be 8-bit
    currentframe=currentframe/max(currentframe(:));
    
    % Maximum intensity of the current frame.
    curr_fr_maxint=1;
    
    % Loop for all guess peaks in the current frame
    for jj=1:numel(guesses{ii})/2
        pix_y=guesses{ii}(2,jj);
        pix_x=guesses{ii}(1,jj);
        if pix_y>0&&pix_x>0&&pix_x<sz(2)&&pix_y<sz(1)
            currentframe=makebox(currentframe,curr_fr_maxint/2,...
                pix_x,pix_y,half_symbol_size,sz);
        end
    end
    
    % Indexing which rows in the fits file belong to the current frame
    r=find(fitfile1(:,1)==ii);
    
    % Loop for all good fits in the current frame
    for jj=1:numel(r)
        pix_y=round(fitfile1(r(jj),11));
        pix_x=round(fitfile1(r(jj),9));
        if pix_y>0&&pix_x>0&&pix_x<sz(2)&&pix_y<sz(1)
            currentframe=makebox(currentframe,curr_fr_maxint,...
                pix_x,pix_y,half_symbol_size,sz);
        end
    end
    
    if numel(trfile)>1
        tr=trfile(trfile(:,2)==ii,:);
        
        cm=jet(5);
        currentframe=currentframe(:,:,[1,1,1]);
        for jj=1:size(tr,1)
            pix_x=round(tr(jj,4));
            pix_y=round(tr(jj,5));
            currentframe=makebox(currentframe,cm(mod(tr(jj,1),size(cm,1)-1)+1,:),...
                pix_x,pix_y,half_symbol_size,sz);
        end
    end
    % magnification
    if magfactor~=1
        for jj=1:size(currentframe,3)
            currentframe2(:,:,jj)=kron(currentframe(:,:,jj),ones(magfactor));
        end
    else
        currentframe2=currentframe;
    end
    
    writeVideo(wobj,currentframe2)
end
if ishandle(h1)
    close(h1)
end

close(wobj);
end
function currentframe=makebox(currentframe,curr_fr_maxint,...
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
function Plot_3D_fits2(fits, spot_size, pxsize,fname)
% This is a part of the "Astig3D_master_program.m' code. It takes in the
% good fit data, and plot each localization as a sphere of fixed
% (user-specified) size. The output is for preliminary visualization only,
% and each spot is not Gaussian-blurred. The localizations belonging to
% different sROIs and tROIs are colored differently.

% The x- and y-coordinates are shifted w.r.t the smallest x- and y-values
% such that the axes labels won't go extremely large.


% INPUTS:

% fits: good fits generated by the 'Astig3D_master_program.m' code.

% spot_size: Specify how large the spheres (representing localizations) are
% in the plotted figure.

% pxsize: Pixel size in nm/px.


% OUTPUT:

% A .fig file will be generated at the same location as the input good fits
% file, with a file name the same as that of the good fits file appended
% with "_3D.fig".

% Shift the x- and y-coordinates w.r.t. the smallest values
min_x=min(fits(:,9));
min_y=min(fits(:,11));
fits(:,9)=fits(:,9)-min_x;
fits(:,11)=fits(:,11)-min_y;

% Check how many spatial ROIs (sROIs) and temporal ROIs (tROIs) are there.
% Fits of different ROIs will be colored differently.
uniq_sROI = unique(fits(:,22))';
uniq_tROI = unique(fits(:,23))';
ct_sROI = length(uniq_sROI);
ct_tROI = length(uniq_tROI);

% Pre-assign color for each ROI
couleur = hsv(ct_sROI*ct_tROI);

% Pre-assign space to store legend text
legend_text = cell(ct_sROI*ct_tROI, 1);

fits_3D=figure;

curr_ROI_ind=1; % Current ROI index for assigning colors and legend text
for curr_sROI=uniq_sROI % Loop through fits from different spatial ROIs
    for curr_tROI=uniq_tROI % Loop through fits from temporal spatial ROIs
        curr_color=couleur(curr_ROI_ind, :);
        
        % Fits from the current sROI and the current tROI
        curr_fits=fits(fits(:,22)==curr_sROI&fits(:,23)==curr_tROI,:);
        
        curr_x=curr_fits(:,9)*pxsize;
        curr_y=curr_fits(:,11)*pxsize;
        curr_z=curr_fits(:,19)*pxsize;
        
        plot3(curr_x,curr_y,curr_z,'o','MarkerSize',spot_size,...
            'MarkerFaceColor',curr_color,'MarkerEdgeColor',curr_color)
        
        legend_text{curr_ROI_ind}=...
            ['sROI ',num2str(curr_sROI),' tROI ',num2str(curr_tROI)];
        
        hold all
        curr_ROI_ind=curr_ROI_ind+1;
    end
end

xlabel('x (nm)','FontSize',16)
ylabel('y (nm)','FontSize',16)
zlabel('z (nm)','FontSize',16)

grid on
axis tight equal ij
legend(legend_text)
view([0,90])
% Save the figure
saveas(gcf, fname);
close (fits_3D);
end
function [cell_xy,good_cell]=select_cells(PhaseMask,img)

% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

% INPUTS:

% PhaseMask: Numbered phase mask of cell images returned by the 'valley.m'
% code or any other cell segmentation code

% OUTPUTS:

% cell_xy: The x- and y-coordinates clicked on by the users to select cells

% good_cell: The indices of chosen cells

cell_fig_h=figure;
subplot(121);
imshow(img,[],'Border','Tight');
subplot(122);
imshow(PhaseMask~=0,[],'Border','Tight');
title(sprintf(['Left-click on cells to be analyzed.\n Press ''Enter'' to'...
    ' proceed.\n Or click return without clicking\n on any cell to analyze ALL of them.']))

% Now move the axes slightly so that the top of the title is visible
set(cell_fig_h,'Units','normalized')
P=get(cell_fig_h,'Position');
set(cell_fig_h,'Position',P+[0 -0.03 0 0])

conv_hull=regionprops(PhaseMask,'Convexhull');
for i=1:length(conv_hull)
    hold all,
    plot(conv_hull(i,1).ConvexHull(:,1),...
        conv_hull(i,1).ConvexHull(:,2),'c-','linewidth',2)
end

key_input=0;
% Let the user pick cells with good shapes
cell_x=[]; cell_y=[];
while key_input~=121
    [curr_cell_x,curr_cell_y,key_input]=jsbginput(1);
    hold all,
    plot(curr_cell_x,curr_cell_y,'m*','markersize',12);
    cell_x=vertcat(cell_x,curr_cell_x);  %#ok<AGROW>
    cell_y=vertcat(cell_y,curr_cell_y); %#ok<AGROW>
end

if ~isempty(cell_x) % If the user did click on something
    cell_xy=[cell_x,cell_y];
    
    good_cell=PhaseMask(sub2ind(size(PhaseMask),round(cell_y),round(cell_x)));
    good_cell=good_cell(good_cell~=0);
    if ~isempty(good_cell) % If the user did click on at least 1 cell
        good_cell=unique(good_cell);
    else % If the user only clicked on background, use all cells.
        good_cell=1:max(PhaseMask(:));
    end
else % If the user hits "Enter" directly, use all cells.
    cell_xy=[];
    good_cell=1:max(PhaseMask(:));
end
close(cell_fig_h);
end