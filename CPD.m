function msds=CPD(trFileName)
max_tau=5;
nMobile=1;
nImmobile=0;

% use overlapping displacements?
yesoverlap=0;

% minimum track length
minTrLength=5;

[cpdfun,pstart]=cpd_function(nMobile,nImmobile);

% options=statset('nlinfit');
% options.RobustWgtFun='bisquare';

plot_tau=[1,2,3];

if ~nargin
    [trFileName,trFileLoc,idf]=uigetfile({'*.mat'},'Select analysis files.',...
        'MultiSelect','on');
    if ~idf
        display('no files chosen.')
        return
    end
    if ~iscell(trFileName)
        trFileName={trFileName};
    end
end

temp=cell(1,max_tau); counter=0;
for kk=1:numel(trFileName)
    m1=matfile([trFileLoc,filesep,trFileName{kk}]);
    try
        trfile=m1.trfile;
    catch
        display([trFileName{kk} ' does not include tracking data. Skipping'])
        continue
    end
    
    for ii=unique(trfile(:,1))'
        counter=counter+1;
        
        trackii=trfile(trfile(:,1)==ii,:);
        for jj=1:max_tau
            if yesoverlap
                indvec1=jj:size(trackii,1);
                indvec2=1:size(trackii,1)-jj;
            else
                indvec2=1:jj:size(trackii,1);
                indvec1=indvec2(1:end-1);
                indvec2=indvec2(2:end);
            end
            
            whichbad=trfile(indvec2,2)-trfile(indvec1,2)~=jj;
            indvec1(whichbad)=[];
            indvec2(whichbad)=[];
            
            temp{counter,jj}=sum((trackii(indvec1,4:5)-trackii(indvec2,4:5)).^2,2);
        end
    end
end

sqsteps=cell(max_tau,1);
for ii=1:max_tau
    sqsteps{ii}=sort(cat(1,temp{:,ii}));
end
ranks=cellfun(@(x)linspace(0,1,numel(x))',sqsteps,'uniformoutput',0);

color_ind=0;
for ii=1:max_tau
    if numel(ranks{ii})<minTrLength
        continue
    end
    
    mdl=fitnlm(sqsteps{ii},ranks{ii},cpdfun,pstart);
    msds(ii,:)=mdl.Coefficients{:,1};
        
    cpdfxn=cpdfun(msds(ii),sqsteps{ii});
    residual=ranks{ii}-cpdfxn(:);
    
    % Create colormap for plotting raw CPD data
    cmap=hsv(numel(plot_tau));
    
    %  plotting CPD at specified tau
    if ismember(ii,plot_tau)
        color_ind=color_ind+1;
                
        % Plot CPDs
        subplot(50,1,1:40)
        semilogx(sqsteps{ii},ranks{ii},'.','Color',cmap(color_ind,:),...
            'MarkerSize',8);
        hold all
        
        % Plot fitted lines
        semilogx(sqsteps{ii},cpdfxn,'-','color','k','Linewidth',2,...
            'HandleVisibility','off')
        
        % Plot residuals
        subplot(50,1,41:50)
        semilogx(sqsteps{ii},residual,'.','Color',cmap(color_ind,:),...
            'MarkerSize',4);
        hold all
        
        subplot(50,1,1:40)
        set(gca,'XTickLabel',[])
    end
end
end
%% ------------------------------------------------------------------------
%  Available CPD Model Functions
%  ------------------------------------------------------------------------
function [fhandle,pstart]=cpd_function(nMobile,nImmobile)
switch nMobile
    case 1
        switch nImmobile
            case 0
                fhandle=@(p,sqr)1-...
                    exp(-sqr/abs(p(1)));
                pstart=1;
            case 1
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/(abs(p(2))+abs(p(3))))-...
                    (1-abs(p(1)))*exp(-sqr/abs(p(3)));
                pstart=[.5,1,.01];
        end
    case 2
        switch nImmobile
            case 0
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/abs(p(2)))-...
                    (1-abs(p(1)))*exp(-sqr/abs(p(3)));
                pstart=[.5,1,.1];
            case 1
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/(abs(p(3))+abs(p(5))))-...
                    abs(p(2))*exp(-sqr/(abs(p(4))+abs(p(5))))-...
                    (1-abs(p(1))-abs(p(2)))*exp(-sqr/(abs(p(5))));
                pstart=[.3,.3,1,.1,.01];
        end
    case 3
        switch nImmobile
            case 0
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/abs(p(3)))-...
                    abs(p(2))*exp(-sqr/abs(p(4)))-...
                    (1-abs(p(1))-abs(p(2)))*exp(-sqr/abs(p(5)));
                pstart=[.3,.3,1,.1,.01];
            case 1
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/(abs(p(4))+abs(p(7))))-...
                    abs(p(2))*exp(-sqr/(abs(p(5))+abs(p(7))))-...
                    abs(p(3))*exp(-sqr/(abs(p(6))+abs(p(7))))-...
                    (1-abs(p(1))-abs(p(2))-abs(p(3)))*exp(-sqr/(abs(p(7))));
                pstart=[.25,.25,.25,1,.1,.01,.001];
        end
end
end