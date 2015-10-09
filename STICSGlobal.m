function fCoeff=STICSGlobal(v)

nDiffs=1;
immPop=0;
diffModel='unconfined';

% 3 elements if unconfined, 5 if confined.
pStart=[.01,.01,.1,0,0];

% filters to use when calculating the correlation functions
% 1:
% 2:
% 3: pad with mean
% 4: non-overlapping time lags;
% 5: pad with zeros
sticsFilters=[3];

% maximum frame separation to consider
maxTau=20;

% magnification of the microscope in micrometers/pixel.
mpp=.049;

% time from beginning of one frame to the beginning of the next in seconds
intTime=.04;

% Choose the functions for fitting
[corrFun,msdFun,pStart,corrLB,corrUB]=sticsFunFinder(nDiffs,immPop,diffModel,pStart);

% gooooooooooood luck figuring this out. nanLinCell linearizes a cell array
% into a single 1D vector. nanLinCell and longmsd both have to be aux
% functions in the code that fHandle and eHandle are used in.
fHandle=@(p,tau,corrDom,corrVal)nanLinCell(cellfun(@(x,y)x-y,...
    cellfun(@(x,y)corrFun(x,y,p),corrDom,num2cell(msdFun(tau,p),1),...
    'uniformoutput',0),corrVal,'uniformoutput',0));
eHandle=@(p,tau,corrDom)cellfun(@(x,y)corrFun(x,y,p),...
    corrDom,num2cell(msdFun(tau,p),1),'uniformoutput',0);

%% calculate correlations
tCorr=stics3(v,[],maxTau,sticsFilters);

%% fit correlations
wSizeM=size(tCorr)*mpp;
wSizeP=size(tCorr);

% space correlations for each tau must be in a cell element by themselves
tCorr=squeeze(mat2cell(tCorr(:,:,2:end),wSizeP(1),wSizeP(2),ones(1,maxTau)))';

% space domain
[X,Y]=ndgrid(linspace(-wSizeM(1)/2,wSizeM(1)/2,wSizeP(1)),...
    linspace(-wSizeM(2)/2,wSizeM(2)/2,wSizeP(2)));
corrDom=cat(3,X,Y);
corrDom={corrDom}; 
corrDom=corrDom(1,ones(1,maxTau));

% time domain
tau=(1:maxTau)*intTime;

% fit
pStart(3)=max(tCorr{1}(:))-min(tCorr{1}(:));
display('fitting...')
fCoeff=lsqnonlin(@(p)fHandle(p,tau,corrDom,tCorr),pStart,corrLB,corrUB);

%% Plot the results
predCorr=eHandle(fCoeff,tau,corrDom);
startCorr=eHandle(pStart,tau,corrDom);
residCorr=cellfun(@(x,y)x-y,tCorr,predCorr,'uniformoutput',0);

pTau=3;
for kk=1:pTau
    subplot(3,pTau,kk); pcolor(tCorr{(kk-1)*5+1});
    axis image; shading flat; axis off
    cl=get(gca,'clim');
    
    if kk==2
        title(num2str(fCoeff(1)))
    end
    
    subplot(3,pTau,kk+pTau); pcolor(residCorr{(kk-1)*5+1});
    axis image; shading flat; axis off
    cl=cl-min(cl)+nanmin(nanmin(residCorr{(kk-1)*5+1}));
    set(gca,'clim',cl);
    
    subplot(3,pTau,kk+2*pTau); pcolor(startCorr{(kk-1)*5+1});
    axis image; shading flat; axis off
    cl=cl-min(cl)+nanmin(nanmin(startCorr{(kk-1)*5+1}));
    set(gca,'clim',cl);
end

end

% square confinement model
function z=longmsd2d(p,x)
% global camerasd
d=p(1); l=p(2); tau=x;

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
function out=nanLinCell(in)
out=cat(1,in{:});
out=out(~isnan(out));
end