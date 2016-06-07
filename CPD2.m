function [dOut,aOut,BIC] = CPD2(mainfold,varargin)
% tracking data in the analysis files must be in pixels.

opts = optimset('Display','off');
%#ok<*PFBNS>

nY = @(m,v,x)exp(-(x-m).^2/2/v)/sqrt(2*pi*v);

if ischar(mainfold)
    %% find the relevant files
    display('Select the analysis files.')
    [datalist,dataloc,~]=uigetfile([mainfold, filesep, '*.mat'],...
        'Matlab Files','multiselect','on');
    
    if ~dataloc
        display('no data selected')
        return
    end
    
    tr=cell(1,numel(datalist)); tic
    for ii=1:numel(datalist)
        m=matfile(fullfile(dataloc,datalist{ii}));
        tr{ii}=m.trackfile;
    end; toc
else
    tr=mainfold;
    if ~iscell(tr)
        tr={tr};
    end
end

%% Default analysis parameters
anProp.nMobile = 1;     % number of diffusive populations
anProp.immBool = 0;     % presence or absence of stationary population
anProp.tFrame = .04;    % camera integration time in seconds
anProp.pixSize = .049;  % pixel size in microns
anProp.maxTau = 5;      % maximum time lag in frames
anProp.confBool = 0;    % confined or unconfined diffusion
anProp.globBool = 1;    % global or local cpd fit
anProp.overBool = 0;    % use overlapping or non-overlapping displacements?
anProp.plotBool = 0;    % plot output or not
anProp.dim = 2;         % 1d or 2d diffusion analysis
anProp.whichDim = 1;    % for 1d diffusion, which dimension to consider
anProp.rotAngle = pi/3; % for 1d diffusion, clockwise angle to rotate the trajectory
anProp.bootNum = 200;   % number of bootstraps

fNames=fieldnames(anProp);

% if any analysis parameters are included as inputs, change the analysis
% parameters mentioned
if nargin>1&&rem(nargin,2)==1
    for ii=1:2:nargin-1
        whichField = strcmp(fNames,varargin{ii});
        
        if all(~whichField)
            warning(['Check spelling. ', ...
                'Parameter change may have not occurred'])
        end
        
        eval([fNames{whichField} ' = varargin{ii+1};'])
        eval(['anProp.' fNames{whichField} ' = ' fNames{whichField},';'])
    end
elseif nargin>1
    warning('use paired inputs')
    dOut = [];
    return
end

% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% partial 1d cpd function
c1=@(x,y,p)p*erf(sqrt(x./2/y));

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% 1d unconfined msd function
m1=@(t,p)2*p(1)*t+p(2);

%% calculate and accumulate squared step sizes
for kk=1:numel(tr)
    trackNums = unique(tr{kk}(:,1))';
    
    if ~isempty(trackNums)
        for ii = trackNums
            tracks = tr{kk}(tr{kk}(:,1)==ii,[2,4,5]);
            
            % fill in the time holes with nans
            fixedTrack = nan(max(tracks(:,1)),size(tracks,2));
            fixedTrack(tracks(:,1),:) = tracks;
            
            % remove leading nans
            fixedTrack(1:find(all(isnan(fixedTrack),2)==0,1,'first')-1,:) = [];
            
            nLocs = size(fixedTrack,1);
            for jj=1:anProp.maxTau
                if anProp.overBool
                    indvec1=jj+1:nLocs;
                    indvec2=1:nLocs-jj;
                elseif ~anProp.overBool
                    indvec2=1:jj:nLocs;
                    indvec1=indvec2(1:end-1);
                    indvec2=indvec2(2:end);
                end
                
                % calculate squared step sizes
                if anProp.dim == 2
                    allSqSteps{kk,ii,jj}=nansum( (fixedTrack(indvec1,[2,3]) - ...
                        fixedTrack(indvec2,[2,3])).^2, 2);
                elseif anProp.dim == 1
                    allSqSteps{kk,ii,jj}=(fixedTrack(indvec1,[2,3]) - ...
                        fixedTrack(indvec2,[2,3])).^2;
                end
            end
        end
    end
end

if anProp.dim == 2
    sqSteps=cell(anProp.maxTau,1);
    for ii=1:anProp.maxTau
        wSteps = cat(1,allSqSteps{:,:,ii});
        sqSteps{ii}=sort(wSteps(wSteps > eps));     % nansum puts zeros where there were nans
    end
    oRanks=cellfun(@(x)linspace(0,1,numel(x))',sqSteps,'uniformoutput',0);
else
    rMat = [cos(anProp.th),-sin(anProp.th);sin(anProp.th),cos(anProp.th)];
    
    sqSteps=cell(anProp.maxTau,1);
    for ii=1:anProp.maxTau
        y=sort(cat(1,allSqSteps{:,ii}));
        rPoints=y*rMat;
        sqSteps{ii}=rPoints(:,anProp.whichDim);
    end
    
    oRanks=cellfun(@(x)linspace(0,1,size(x,1))',sqSteps,'uniformoutput',0);
end
nSteps = cellfun(@numel,sqSteps,'uniformoutput',0);
nTotal = sum(cat(1,nSteps{:}));

%% fitting function selection
funFinds = cpdFunFinder(anProp)
cpdFun = funFinds.cpdFun;
msdFun = funFinds.msdFun;
pStart = funFinds.pStart;
bounds = funFinds.bounds;
dID = funFinds.dID;
aID = funFinds.aID;

% time lag domain in seconds
tau = (1:anProp.maxTau)'*anProp.tFrame;

if anProp.globBool
    %% GLOBAL FITTING
    linCell=@(x)cat(1,x{:});
    fHandle=@(p,tau,sqSteps,ranks)linCell(...
        cellfun(@(x,y)x-y,...
        cellfun(@(x,y)cpdFun(x,y,p),...
        sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0),...
        ranks,'uniformoutput',0));
    eHandle=@(p,tau,sqSteps)cellfun(@(x,y)cpdFun(x,y,p),...
        sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0);
    
    fParams = zeros(numel(pStart{1}),anProp.bootNum); 
%     resnorm = zeros(anProp.bootNum,1);
    resids = zeros(nTotal,anProp.bootNum);
    parfor kk = 1:anProp.bootNum
        % converts track positions from pixels to microns and resamples at the same
        % time
        y=cellfun(@(x,y)sort(x(randsample(y,y,1))*anProp.pixSize.^2),sqSteps,nSteps,'uniformoutput',0);
        
        [fParams(:,kk),~,resids(:,kk)] = lsqnonlin(@(p)fHandle(p,tau,y,oRanks),...
            pStart{1},bounds{1},bounds{2},opts);
    end
    dOut = fParams(dID,:);
    aOut = fParams(aID,:);
    BIC = nTotal*log(mean(resids.^2,1)) + numel(pStart{1})*log(nTotal);
elseif ~anProp.globBool
    %% LOCAL FITTING
    cpdLB = bounds{1};
    cpdUB = bounds{2};
    msdLB = bounds{3};
    msdUB = bounds{4};
    
    fParams1a = zeros(numel(pStart{1}),anProp.maxTau,anProp.bootNum);
    aOut = zeros(numel(aID),anProp.maxTau,anProp.bootNum);
    resnorm = zeros(anProp.bootNum,anProp.maxTau);
    parfor kk = 1:anProp.bootNum
        cpdStart = pStart{1};
        msdStart = pStart{2};
        
        fParams1 = zeros(numel(cpdStart),anProp.maxTau);
        fParams2 = zeros(numel(msdStart),anProp.nMobile);
        
        y=cellfun(@(x,y)sort(x(randsample(y,y,1))*anProp.pixSize.^2),sqSteps,nSteps,'uniformoutput',0);
        
        % fit the cpd curves to get the msd values
        for mm=1:anProp.maxTau
            for jj = 1:anProp.nMobile
                cpdStart(jj) = mean(y{mm})/(10^(jj-1))*anProp.tFrame;
            end
            [fParams1(:,mm),~] = lsqcurvefit(@(p,x)cpdFun(x,p),...
                cpdStart,y{mm},oRanks{mm},cpdLB,cpdUB,opts);
        end
        
        fParams1a(:,:,kk) = fParams1;
        
        % fit msds to get diffusion coefficients
        for mm = 1 : anProp.nMobile
            fParams2(:,mm)=lsqcurvefit(@(p,x)msdFun(x,p),msdStart,tau,fParams1(mm,:)',...
                msdLB,msdUB,opts);
            %             resids{2,ii} = msdFun(tau,fParams(ii,:)) - msds(ii,:)';
        end
        dOut(:,kk) = fParams2(dID,:);
        aOut(:,:,kk) = fParams1(aID,:);
    end
    
    % calculate BIC
%     mean(-log(nY(mean(x),var(x),x)))
    BIC = nTotal*log(resnorm.^2) + numel(pStart{1})*log(nTotal);
end

%% plot fits
if anProp.plotBool
    %     figure
    %     if anProp.globBool
    %         fRes = eHandle(fParams(:,1),tau,y);
    %         for ii = 1:anProp.maxTau
    %             subplot(1,anProp.maxTau,ii)
    %             plot(y{ii},cat(2,oRanks{ii},fRes{ii}),'.')
    %             set(gca,'xscale','log')
    %         end
    %
    %     elseif ~anProp.globBool
    %         for ii = 1:anProp.maxTau
    %             subplot(1,anProp.maxTau,ii)
    %             plot(y{ii},cat(2,oRanks{ii},cpdFun(y{ii},fParams1a(1,ii,:))),'.')
    %             set(gca,'xscale','log')
    %         end
    %     end
    subplot(131)
    for ii=1:size(dOut,1)
        %     subplot(1,numel(dID),ii)
        ecdf(dOut(ii,:)); hold all
    end; hold off
    xlabel('d in microns^2/s')
    set(gca,'xscale','log')
    
    subplot(132)
    for ii = 1:size(aOut,1)
        ecdf(aOut(ii,:)); hold all
    end; hold off
    xlabel('relative population amplitudes')
    
    subplot(133)
    for ii = 1:anProp.bootNum
        plot(resids(:,ii)); hold all
    end; hold off
end
end