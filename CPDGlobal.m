function [x,BIC]=CPDGlobal(trfile,nDiffs,immPop,msdModel)
% example: x=CPDGlobal(testCPD(2),1,0,'unconfined');
% computes msds from cumulative probability step size distributions, and
% then fits the msds to a msd model to estimate the diffusion coefficients
%
% to use, run this without inputs and select the analysis files or run this
% with one input, which is the particular tracking file you want to
% analyze.
% opts=optimset('Display','off');

%% Experimental constants

% integration time in seconds
intTime=.04;

% pixel size in microns in the object plane
mpp=.049;

%% Experimental parameters
% 4 different parameters each with 2 values = 2^4 different combos

% diffusion dimensionality (1 or 2)
dim=2;

% direction of unit vector along the dimension you wish to consider (if
% dim=1), in radians
dimAngle=0;%pi/2;

% number of diffusive populations in the CPD model (1 or 2)
% nDiffs=1;

% presence or absence of immobile population (0 or 1)
% immPop=0;

% diffusion model. ('confined' or 'unconfined')
% msdModel='unconfined';

%% Algorithm parameters

% maximum time lag in frames
maxTau=10;

% use overlapping displacements to calculate the step size distributions?
yesOverlap=1;

% minimum track length (integer)
minTrLength=5;

% aux function
linCell=@(x)cat(1,x{:});

% which time-lags to plot in the cpd fit result figure
plotTau=1:maxTau;

%% Choose the functions for fitting
pStart= [.1,.5,.01,...          % msd 1 parameters D1, L, S1
        .01,.01,...             % msd 2 parameters D2, S2
        .3,.3,.01];             % cpd function parameters:
                                % amp1, amp2, immSize
cpdLB=zeros(1,numel(pStart));
cpdLB([3,5])=-inf;

cpdUB=inf(1,numel(pStart));

% restricted maximum size of diffusion coefficient in microns^2/sec
cpdUB([1,4])=10;

% restricted maximum size of confinement in microns
cpdUB([2,4])=20;

[cpdFun,msdFun,pId]=cpdFunFinder(dim,nDiffs,immPop,msdModel);
pStart=pStart(pId);
cpdLB=cpdLB(pId);
cpdUB=cpdUB(pId);
ds=ismember(pId,[1,4]);
amps=ismember(pId,[6,7]);

% gooooooooooood luck figuring this out. linCell linearizes a cell array
% into a single 1D vector. linCell and longmsd both have to be aux
% functions in the code that fHandle and eHandle are used in.
fHandle=@(p,tau,sqSteps,ranks)linCell(...
    cellfun(@(x,y)x-y,...
    cellfun(@(x,y)cpdFun(x,y,p),...
    sqSteps,num2cell(msdFun(tau,p),1)','uniformoutput',0),...
    ranks,'uniformoutput',0));
eHandle=@(p,tau,sqSteps)cellfun(@(x,y)cpdFun(x,y,p),...
    sqSteps,num2cell(msdFun(tau,p),1)','uniformoutput',0);

%% Get the locations and names of all the analysis files
if ~nargin
    [trFileName,trFileLoc,idf]=uigetfile({'*.mat'},'Select analysis files.',...
        'MultiSelect','on');
    if ~idf
        display('no files chosen.')
        return
    end
    
    % convert the name to a cell array if only one file is chosen
    if ~iscell(trFileName)
        trFileName={trFileName};
    end
    
else
    % if there's an input, don't load anything
    trFileName={'manual input'};
end

%% Calculate and accumulate squared step sizes
allSqSteps=cell(1,maxTau); counter=0;
for kk=1:numel(trFileName)
    
    % if there's no input, load the files, if there's an input, just use
    % the input as the tracking file
    if ~nargin
        m1=matfile([trFileLoc,filesep,trFileName{kk}]);
        try
            trfile=m1.trackfile;
        catch
            display([trFileName{kk} ' does not include tracking data. Skipping'])
            continue
        end
    end
    
    % loop through tracks
    for ii=unique(trfile(:,1))'
        % select current track from the array of all tracks
        trackii=trfile(trfile(:,1)==ii,:);
        
        if size(trackii,1)<minTrLength
            continue
        end
        
        counter=counter+1;
        
        % fill in the time holes with nans. these two lines are genius!
        fixedTrack=nan(max(trackii(:,2)),size(trackii,2));
        fixedTrack(trackii(:,2),:)=trackii;
        
        % remove leading nans
        fixedTrack(1:find(all(isnan(fixedTrack),2)==0,1,'first')-1,:)=[];
        
        for jj=1:maxTau
            if yesOverlap                       % overlapping frame pairs
                indvec1=jj+1:size(fixedTrack,1);
                indvec2=1:size(fixedTrack,1)-jj;
            else                                % nonoverlapping frame pairs
                indvec2=1:jj:size(fixedTrack,1);
                indvec1=indvec2(1:end-1);
                indvec2=indvec2(2:end);
            end
            
            if dim==1
                rTrack=nansum(bsxfun(@times,fixedTrack(:,4:5),...
                    [cos(dimAngle),sin(dimAngle)]),2);
                
                allSqSteps{counter,jj}=(rTrack(indvec1)-rTrack(indvec2)).^2;
            else
                % nansum because there are nans as placeholders
                allSqSteps{counter,jj}=nansum((fixedTrack(indvec1,4:5)-...
                    fixedTrack(indvec2,4:5)).^2,2);
            end
        end
    end
end

sqSteps=cell(maxTau,1);
for ii=1:maxTau
    sqSteps{ii}=sort(cat(1,allSqSteps{:,ii}));
end
ranks=cellfun(@(x)linspace(0,1,numel(x))',sqSteps,'uniformoutput',0);
tau=(1:maxTau)*intTime;

%% Fit the cpd curves to get the diffusion coefficient
x=lsqnonlin(@(p)fHandle(p,tau,sqSteps,ranks),pStart,cpdLB,cpdUB);

d=x(ds);

%% Plot the results
fRanks=eHandle(x,tau,sqSteps);
residCPD=cellfun(@(x,y)x-y,ranks,fRanks,'uniformoutput',0);

color_ind=0;
cmap=hsv(numel(plotTau));
for ii=1:maxTau
    color_ind=color_ind+1;
    
    % Plot CPDs
    subplot(50,1,1:40); 
    semilogx(sqSteps{ii},ranks{ii},'.','Color',cmap(color_ind,:),...
        'MarkerSize',8);
    set(gca,'XTickLabel',[])
    ylabel('Cumulative probability')
    hold on
    
    % Plot fitted lines
    semilogx(sqSteps{ii},fRanks{ii},'--','color','k','Linewidth',1,...
        'HandleVisibility','off')
    axis tight
    xL=get(gca,'xlim');
    set(gca,'ylim',[0,1]);
    
    % Plot residuals
    subplot(50,1,41:50)
    semilogx(sqSteps{ii},residCPD{ii},'.','Color',cmap(color_ind,:),...
        'MarkerSize',4);
    set(gca,'xlim',xL)
    set(gca,'ylim',[-.05,.05]);
    
    xlabel('Squared step size')
    hold on
end

if any(ismember(1:maxTau,plotTau))
    subplot(50,1,1:40)
    title(num2str(d))
    hold off
    
    subplot(50,1,41:50)
    hold off
end

n=sum(cellfun(@numel,sqSteps));
k=numel(x);
s2=sum(cat(1,residCPD{:}).^2)/n;

BIC=n*log(s2)+k*log(n);
end

% %% square confinement model
% function z=longmsd2d(p,x)
% % global camerasd
% d=p(1); l=p(2); tau=x;
% 
% summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);
% 
% temp=eps*ones(size(tau));
% for ii=1:2:2*400-1
%     s=summedterm(tau,d,l,ii);
%     if sum(s./temp)<1e-10
%         break
%     end
%     temp=temp+s;
% end
% z=l^2/3*(1-96/pi^4*temp)+p(3);
% end
% 
% %% square confinement model
% function z=longmsd1d(p,x)
% % global camerasd
% d=p(1); l=p(2); tau=x;
% 
% summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);
% 
% temp=eps*ones(size(tau));
% for ii=1:2:2*400-1
%     s=summedterm(tau,d,l,ii);
%     if sum(s./temp)<1e-10
%         break
%     end
%     temp=temp+s;
% end
% z=l^2/6*(1-96/pi^4*temp)+p(3);
% end