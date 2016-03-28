function resStruct = dataParse(data, mask)
%
% NAME:
%       dataParse
% PURPOSE:
%       Analyze single-molecule imaging data with tracking or STICS
% CATEGORY:
%       Image Processing
% CALLING SEQUENCE:
%       resStruct = dataParse(data,flag);
% INPUTS:
%       data:               x by y by t movie of fluorescence data
%       mask:               (optional) 2d binary roi selection mask
% OUTPUTS:
%       resStruct:          fitting result 'structure'
% PROCEDURE:
%       1. Peak fitting, then MSD calculation a la single molecule analysis
%       2. Correlation function calculation and width estimation
%       3. Diffusion coefficient estimation by MSD fitting
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
% NOTES:
%       This code 'dataParse.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%
%       For testing purposes, run this script:
%
%       v = dataGen('confined',1);
%       r = dataParse(v);

%% analysis parameters
tFrame = 0.05;          % camera integration time in seconds
nTau = 15;              % number of smallest time lag values to use (excluding 0)
pixelSize = .049;       % pixel size in microns
nFrames = size(data,3); % number of frames in the data
yesOverlap = 1;         % use overlapping time lags?
isConfined = 1;         % use confined MSD curve fit?
padStyle = 1;           % padding style for STICS analysis
                        % 1 : pad zeros, 2: pad mean, 3: pad mean outside mask

%% single molecule tracking analysis

% spot fitting. use parfor to suppress figure generation or manually do it.
fitP = zeros(nFrames, 6);
parfor ii = 1:nFrames
    fitP(ii,:) = gaussFit(data(:,:,ii), 1);
end

% 'tracking'. missed spots will register as nans.
tr = fitP(:,1:2) * pixelSize;

% calculate mean squared displacements, one for each time lag value
MSDs=zeros(nTau,2);
for ii=1:nTau
    if yesOverlap                       % overlapping frame pairs
        indvec1=ii+1:nFrames;
        indvec2=1:nFrames-ii;
    elseif ~yesOverlap                  % nonoverlapping frame pairs
        indvec2=1:ii:nFrames;
        indvec1=indvec2(1:end-1);
        indvec2=indvec2(2:end);
    end
    
    % mean squared displacements vs time lag
    MSDs(ii,:)=nansum((tr(indvec2,:)-tr(indvec1,:)).^2,1);
end

%% STICS analysis

% pad images
switch padStyle
    case 1
        data = padZeros(data);
    case 2
        data = padMean(data);
    case 3
        data = padMask(data,mask);
end

% time-space correlation function calculation
if yesOverlap                           % overlapping frame pairs
    famps=abs(fft(fft(fft(data,[],1),[],2),[],3)).^2;
    STCorr = fftshift(fftshift(real(ifft(ifft(ifft(famps...
        ,[],1),[],2),[],3)),1),2)/numel(famps)/mean(data(:))^2-1;
    STCorr = STCorr(:,:,2:nTau+1);
    
elseif ~yesOverlap                      % nonoverlapping frame pairs
    vFft=fft2(data);
    STCorr=zeros(size(vFft,1),size(vFft,2),nTau+1);
    
    STCorr(:,:,1)=mean(bsxfun(@times,real(fftshift(fftshift(...
        ifft2(vFft.*conj(vFft)),1),2)),1./(mean(mean(data)).^2)),3);
    
    for kk=1:nTau
        ind1 = 1:kk:vidsize(3);
        ind2 = ind1(2:end);
        ind1 = ind1(1:end-1);
        STCorr(:,:,kk+1) = mean(bsxfun(@times,real(fftshift(fftshift( ...
            ifft2(vFft(:,:,ind2).*conj(vFft(:,:,ind1))),1),2)), ...
            1./mean(mean(data(:,:,ind1)))./mean(mean(data(:,:,ind2)))),3);
    end
    STCorr=STCorr/numel(vFft(:,:,1))-1;
end

% estimate the widths of the correlation function
iMSDs=zeros(nTau,2); MSDd=iMSDs;
[x,y]=ndgrid(1:size(STCorr,1),1:size(STCorr,2));
for ii = 1:nTau
    fitP = gaussFit(STCorr(:,:,ii),0);
    iMSDs(ii,:) = fitP(3:4).^2 * pixelSize^2;
    
    % discrete variance calculation for x dimension
    pmf=sum(STCorr(:,:,ii),2);
    MSDd(ii,1) = sum(pmf.*(x(:,1)-mean(x(:))).^2) * pixelSize^2;
    
    % discrete variance calculation for x dimension
    pmf=sum(STCorr(:,:,ii),1);
    MSDd(ii,2) = sum(pmf.*(y(1,:)-mean(y(:))).^2) * pixelSize^2;
end

%% MSD Fitting

% choose fitting function
if isConfined
    f=@(p,X)sqconfMSD1D(p,X);
    pStart = [1, .1, 0];
    lb = [0, 0, -inf];
    ub = [inf, inf, inf];
else
    f=@(p,X) 2*p(2)*X+p(1);
    pStart = [0, .1];
    lb = [-inf, 0];
    ub = [inf, inf];
end

% time lag domain vectgor
tau = (1:nTau)'*tFrame;

% D from tracking
pT(:,1)=lsqcurvefit(f,pStart,tau,MSDs(:,1),lb,ub);
pT(:,2)=lsqcurvefit(f,pStart,tau,MSDs(:,2),lb,ub);

% D from STICS
pS(:,1)=lsqcurvefit(f,pStart,tau,iMSDs(:,1),lb,ub);
pS(:,2)=lsqcurvefit(f,pStart,tau,iMSDs(:,2),lb,ub);

% D from 'discrete variance'
pD(:,1)=lsqcurvefit(f,pStart,tau,MSDd(:,1),lb,ub);
pD(:,2)=lsqcurvefit(f,pStart,tau,MSDd(:,2),lb,ub);

resStruct.Dtracking = pT(2,:);
resStruct.Dstics =  pS(2,:);
resStruct.Dvar =  pD(2,:);
end

function zi=sqconfMSD1D(p,X)
l=p(1);
d=exp(p(2));
ns=p(3);
tau=X;

summedTerm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));
for ii=1:2:2*400-1
    s=summedTerm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
zi=l^2/6*(1-96/pi^4*temp)+ns;
zi(isnan(zi))=eps;
zi(isinf(zi))=eps;
end

function imStack=padMask(imStack,mask)
% replace pixels outside the mask with the average value inside the mask in
% each frame

imsize=size(imStack);
mMean=mean(reshape(imStack(mask(:,:,ones(1,imsize(3)))),[],imsize(3)));
imStack(~mask(:,:,ones(1,imsize(3))))=mMean(ones(1,sum(~mask(:))),:);
end

function imstack=padMean(imstack)
% pad the first two dimensions to double size with the mean of each image
imstack=mat2cell(imstack,size(imstack,1),size(imstack,2),ones(1,size(imstack,3)));
imm=cellfun(@(x)mean(mean(x(mask))),imstack,'uniformoutput',false);

imstack=cellfun(@(x,y)padarray(x,floor(size(x)/2),y),imstack,imm,...
    'uniformoutput',false);
imstack=cat(3,imstack{:});
end

function imstack=padZeros(imstack)
% pad the first two dimensions to double size with zeros
imstack=mat2cell(imstack,size(imstack,1),size(imstack,2),ones(1,size(imstack,3)));

imstack=cellfun(@(x,y)padarray(x,floor(size(x)/2),y),...
    imstack,repmat({0},1,1,size(imstack,3)),'uniformoutput',false);
imstack=cat(3,imstack{:});
end