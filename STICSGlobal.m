function outPut=STICSGlobal(v,nDiffs,immPop,diffModel)

% nDiffs=1;
% immPop=0;
% diffModel='confined';

% filters to use when calculating the correlation functions. the order is
% preserved during analysis.
% 1: non-overlapping image pairs time lags
% 2: replace out-mask values with in-mask mean of each frame
% 3: lowpass fourier filter
% 4: pad with mean of each frame
% 5: pad with zeros
sticsFilters=[4];

% maximum frame separation to consider
maxTau=20;
startTau=1;

% magnification of the microscope in micrometers/pixel.
mpp=.049;

% time from beginning of one frame to the beginning of the next in seconds
intTime=.04;

% time domain
tau=(startTau:maxTau)*intTime;

% rotation angle
th=0;

% Choose the functions for fitting
pStart=[.1,.5,.01,1,.01,...  % msd 1 parameters      D1, L1, S1, L2, S2
    .1,.01,.01,...          % msd 2 parameters      D2, S3, S4
    0,.01,0,0,.01,.01,.01]; % correlation function parameters:
                            % offset, amp1, xc, yc, amp2, immAmp, immSize
cpdLB([3,5,7,8,9,11,12])=-inf;
cpdLB([1,2,4,6,10,13,14,15])=0;
cpdUB=inf(1,50);

% restricted maximum size of diffusion coefficient
cpdUB([1,6])=10;

% restricted maximum size of confinement
cpdUB([2,4])=20;

[corrFun,msdFun,pId]=sticsFunFinder(nDiffs,immPop,diffModel,th);
pStart=pStart(pId);
cpdLB=cpdLB(pId);
cpdUB=cpdUB(pId);
amps=ismember(pId,[10,13,14]);

% gooooooooooood luck figuring this out. nanLinCell linearizes a cell array
% into a single 1D vector. nanLinCell and longmsd both have to be aux
% functions in the code that fHandle and eHandle are used in.
fHandle=@(p,tau,corrDom,corrVal)nanLinCell(cellfun(@(x,y)x-y,...
    cellfun(@(x,y)corrFun(x,y,p),corrDom,num2cell(abs(msdFun(tau,p)),1),...
    'uniformoutput',0),corrVal,'uniformoutput',0));
eHandle=@(p,tau,corrDom)cellfun(@(x,y)corrFun(x,y,p),...
    corrDom,num2cell(abs(msdFun(tau,p)),1),'uniformoutput',0);

%% calculate correlations
display('correlating...')
tCorr=stics3(v,[],maxTau,sticsFilters);

%% fit correlations
wSizeP=size(tCorr);

% space correlations for each tau must be in a cell element by themselves
ind1=round(wSizeP(1)/2-wSizeP(1)/4):round(wSizeP(1)/2+wSizeP(1)/4);
ind2=round(wSizeP(2)/2-wSizeP(2)/4):round(wSizeP(2)/2+wSizeP(2)/4);
tCorr=squeeze(mat2cell(tCorr(ind1,ind2,startTau+1:end),...
    numel(ind1),numel(ind2),ones(1,maxTau-startTau+1)))';

% space domain
[X,Y]=ndgrid((ind1-mean(ind1))*mpp,(ind2-mean(ind2))*mpp);
corrDom=cat(3,X,Y);
corrDom={corrDom};
corrDom=corrDom(1,ones(1,maxTau-startTau+1));

% fit
display('fitting...')

pStart(amps)=(max(tCorr{1}(:))-min(tCorr{1}(:)))/sum(amps)/8;
fCoeff=lsqnonlin(@(p)fHandle(p,tau,corrDom,tCorr),pStart,cpdLB,cpdUB);

%% Plot the results
predCorr=eHandle(fCoeff,tau,corrDom);
startCorr=eHandle(pStart,tau,corrDom);
residCorr=cellfun(@(x,y)x-y,tCorr,predCorr,'uniformoutput',0);

n=ceil(sqrt(maxTau*2));
for kk=1:maxTau-startTau+1
    subplot(n,n,kk); pcolor(tCorr{kk});
    axis image; shading flat; axis off
    cl=get(gca,'clim');
    
    subplot(n,n,kk+maxTau); pcolor(residCorr{kk});
    axis image; shading flat; axis off
    cl=cl-min(cl)+nanmin(nanmin(residCorr{kk}));
    set(gca,'clim',cl);
    
    %     subplot(n,n,kk+maxTau); pcolor(startCorr{kk});
    %     axis image; shading flat; axis off
    %     cl=cl-min(cl)+nanmin(nanmin(startCorr{kk}));
    %     set(gca,'clim',cl);
end

%% collect and label outputs
% non-existent fitting parameters exist as empty fields in outPut.
outPut.diffC1       =fCoeff(pId==1);    % diffusion coefficient 1
outPut.confX        =fCoeff(pId==2);    % confinement length 1
outPut.locNoise1    =fCoeff(pId==3);
outPut.confY        =fCoeff(pId==4);    % confinement length 2
outPut.locNoise2    =fCoeff(pId==5);
outPut.diffC2       =fCoeff(pId==6);    % diffusion coefficient 2
outPut.locNoise3    =fCoeff(pId==7);
outPut.locNoise4    =fCoeff(pId==8);
outPut.corrOffset   =fCoeff(pId==9);    % correlation offset
outPut.amp1         =fCoeff(pId==10);   % correlation amplitude 1
outPut.corrCenterX  =fCoeff(pId==11);   % correlation center 1
outPut.corrCenterY  =fCoeff(pId==12);   % correlation center 2
outPut.amp2         =fCoeff(pId==13);   % correlation amplitude 2
outPut.ampImm       =fCoeff(pId==14);   % correlation amplitude 3
outPut.immNoise     =fCoeff(pId==15);

end
function out=nanLinCell(in)
out=cat(1,in{:});
out=out(~isnan(out));
end