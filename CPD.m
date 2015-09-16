function [dmeas]=CPD(trfile)
% computes msds from cumulative probability step size distributions, and
% then fits the msds to a msd model to estimate the diffusion coefficients
%
% to use, run this without inputs and select the analysis files or run this
% with one input, which is the particular tracking file you want to
% analyze.
opts=optimset('Display','off');
maxTau=3;         % maximum time lag in frames
intTime=.04;        % integration time in seconds
mpp=.049;           % pixel size in microns in the object plane

% number and type of terms in the CPD fit
nDiffs=2;
immPop=0;
cpdStart=[]; %[.5,0,eps,.1,eps];        % leave blank unless you know what you're doing

% mean squared displacement model for MSD fit
msdModel='unconfined'; % 'square confinement'
msdStart=[];        % leave blank unless you know what you're doing

% use overlapping displacements?
yesOverlap=1;

% minimum track length
minTrLength=5;

% 2 or 1 dimensional diffusion. 1D isn't set up yet
dim=2;

% choose the functions for fitting
[cpdFun,cpdStart,cpdLB,cpdUB]=cpd_function(dim,nDiffs,immPop,cpdStart);
[msdFun,msdStart,msdLB,msdUB]=msd_function(dim,msdModel,msdStart);

% which taus to plot in the cpd fit result figure
plotTau=1:maxTau;   % 1:maxTau; % [];

% which msds to plot
plotMsds=1:nDiffs;  % 1:nTerms; % [];

%% get the locations and names of all the analysis files
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

%% calculate and accumulate squared step sizes
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
        fixedtrack=nan(max(trackii(:,2)),size(trackii,2));
        fixedtrack(trackii(:,2),:)=trackii;
        
        % remove leading nans
        fixedtrack(1:find(all(isnan(fixedtrack),2)==0,1,'first')-1,:)=[];
        
        for jj=1:maxTau
            if yesOverlap                       % overlapping frame pairs
                indvec1=jj+1:size(fixedtrack,1);
                indvec2=1:size(fixedtrack,1)-jj;
            else                                % nonoverlapping frame pairs
                indvec2=1:jj:size(fixedtrack,1);
                indvec1=indvec2(1:end-1);
                indvec2=indvec2(2:end);
            end
            
            % nansum because there are nans as placeholders
            allSqSteps{counter,jj}=nansum((fixedtrack(indvec1,4:5)-...
                fixedtrack(indvec2,4:5)).^2,2)*mpp^2;
        end
    end
end

sqSteps=cell(maxTau,1);
for ii=1:maxTau
    sqSteps{ii}=sort(cat(1,allSqSteps{:,ii}));
end
ranks=cellfun(@(x)linspace(0,1,numel(x))',sqSteps,'uniformoutput',0);

%% fit the cpd curves to get the msd values
color_ind=0;
for ii=1:maxTau
    
    msds(:,ii)=lsqcurvefit(cpdFun,cpdStart,sqSteps{ii},ranks{ii},...
        cpdLB,cpdUB,opts);
    %         mMdl{ii}=fitnlm(sqSteps{ii},ranks{ii},cpdFun,cpdStart);
    %     msds(:,ii)=abs(mMdl{ii}.Coefficients{:,1});
    
    cpdfxn=cpdFun(msds(:,ii),sqSteps{ii});
    residual=ranks{ii}-cpdfxn(:);
    
    % Create colormap for plotting raw CPD data
    cmap=hsv(numel(plotTau));
    
    %  plotting CPD at specified tau
    if ismember(ii,plotTau)
        color_ind=color_ind+1;
        
        % Plot CPDs
        subplot(50,2,sub2ind([2,50],ones(1,40),1:40))
        semilogx(sqSteps{ii},ranks{ii},'.','Color',cmap(color_ind,:),...
            'MarkerSize',8);
        set(gca,'XTickLabel',[])
        ylabel('Cumulative probability')
        hold on
        
        % Plot fitted lines
        semilogx(sqSteps{ii},cpdfxn,'-','color','k','Linewidth',2,...
            'HandleVisibility','off')
        
        % Plot residuals
        subplot(50,2,sub2ind([2,50],ones(1,10),41:50))
        semilogx(sqSteps{ii},residual,'.','Color',cmap(color_ind,:),...
            'MarkerSize',4);
        xlabel('Squared step size')
        hold on
    end
end

if any(ismember(1:maxTau,plotTau))
    subplot(50,2,sub2ind([2,50],ones(1,40),1:40))
    hold off
    subplot(50,2,sub2ind([2,50],ones(1,10),41:50))
    hold off
end

%% fit msds
tau=(1:maxTau)*intTime;
dmeas=nan(nDiffs,numel(msdStart));
counter=0;
for ii=nDiffs+2*immPop:size(msds,1)
    counter=counter+1;
    y=msds(ii,:);
    if ~sum(~isnan(y))
        display(['diffusive population number ', ...
            num2str(ii-nDiffs+1), ' has no msd values'])
        dmeas(counter,:)=nan(size(msdStart));
        
        continue
    end
%     dMdl{counter}=fitnlm(tau(~isnan(y)),y(~isnan(y)),msdFun,msdStart);
%     dmeas(counter,:)=dMdl{counter}.Coefficients{:,1};
%     dmeas95(counter,:)=diff(coefCI(dMdl{counter}),[],1);
    dmeas(counter,:)=lsqcurvefit(msdFun,msdStart,tau,y,...
        msdLB,msdUB,opts);
    
    if ismember(counter,plotMsds)
        subplot(122)
        scatter(tau,y,'fill')
        hold all
        plot(tau,msdFun(dmeas(counter,:),tau),'--k')
        ylabel('Mean Squared Displacement')
        xlabel('time lag')
        
        leg{2*(counter-1)+1}=num2str(dmeas(counter,1));
        leg{2*(counter-1)+2}='fit';
    end
end
legend(leg)
hold off

dmeas=abs(dmeas);

end

%% ------------------------------------------------------------------------
%  Available CPD Model Functions
%  ------------------------------------------------------------------------
function [fhandle,pstart,cpdLB,cpdUB]=cpd_function(dim,nDiffs,immPop,pstart)
switch dim
    case 2
        switch immPop
            case 1
                switch nDiffs
                    case 1
                        fhandle=@(p,sqr)1-...
                            p(1)*               exp(-sqr/(p(3)+p(2)))-...
                            (1-p(1))*           exp(-sqr/p(2));
                        if isempty(pstart)
                            pstart=[.5,.01^2,1]; % 3
                        end
                        cpdLB=[0,0,0];
                        cpdUB=[1,inf,inf];
                        
                    case 2
                        fhandle=@(p,sqr)1-...
                            p(1)*               exp(-sqr/(p(4)+p(3)))-...
                            p(2)*               exp(-sqr/(p(5)+p(3)))-...
                            (1-p(1)-p(2))*      exp(-sqr/p(3));
                        if isempty(pstart)
                            pstart=[.5,.5,eps,.1,.001]; % 4:5
                        end
                        cpdLB=[0,0,0,0,0];
                        cpdUB=[1,1,inf,inf,inf];
                end
            case 0
                switch nDiffs
                    case 1
                        fhandle=@(p,sqr)1-      exp(-sqr/p(1));
                        if isempty(pstart)
                            pstart=.1; % 1
                        end
                        cpdLB=0;
                        cpdUB=inf;
                        
                    case 2
                        fhandle=@(p,sqr)1-...
                            p(1)*               exp(-sqr/p(2))-...
                            (1-p(1))*           exp(-sqr/p(3));
                        if isempty(pstart)
                            pstart=[.5,.1,.001]; % 2:3
                        end
                        cpdLB=[0,0,0];
                        cpdUB=[1,inf,inf];
                        
                end
        end
    case 1
        error('insert 1d cpd model function')
end
end

%% ------------------------------------------------------------------------
%  Available MSD Model Functions
%  ------------------------------------------------------------------------
function [fhandle,pstart,msdLB,msdUB]=msd_function(dim,model,pstart)

switch model
    case 'unconfined'
        switch dim
            case 1
                fhandle=@(p,x) 2*p(1)*x+p(2);
                
            case 2
                fhandle=@(p,x) 4*p(1)*x+p(2);
                
        end
        if isempty(pstart)
            pstart=[.001,0];
        end
        msdLB=[0,-inf];
        msdUB=[inf,inf];
        
    case 'square confinement'
        switch dim
            case 1
                fhandle=@(p,x) longmsd(p,x)/2;
            case 2
                fhandle=@(p,x) longmsd(p,x);
        end
        if isempty(pstart)
            pstart=[1,.1,0];
        end
        msdLB=[0,0,-inf];
        msdUB=[inf,inf,inf];
        
    otherwise
        display('mismelled the confinement model type')
        fhandle=[];
        pstart=[];
end
end

%% square confinement model
function z=longmsd(p,x)
% global camerasd
l=p(1); d=p(2); tau=x;

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));
for ii=1:2:2*400-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
z=l^2/3*(1-96/pi^4*temp)+p(3);
end