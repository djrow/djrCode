function [fitPars, conf95]=gaussFit(img,findTheSpot)
%
% NAME:
%       gaussFit
% PURPOSE:
%       Fits a generalized gaussian function to 2d imaging data. This code
%       produces results in units of pixels for the center position and
%       widths.
% CATEGORY:
%       Image Processing
% CALLING SEQUENCE:
%       [fitPars, conf95] = gaussFit(img,findTheSpot);
% INPUTS:
%       img:            The two-dimensional array to be fit to a gaussian
%       findTheSpot:    (optional) 1 or 0. Default behavior is to fit an
%           ROI in the center of the image. If the spot is not near the
%           center or the image is very large, findTheSpot enables the code
%           to first roughly locate the spot and then use that location as
%           the ROI center.
% OUTPUTS:
%       mdl:            fitting result "object"
%       fitPars:        fitting coefficient vector
%       fitCI:          95% confidence interval of fitting coefficients at
%                       end of fitting
% PROCEDURE:
%       1. Peak guessing and/or data ROI selection of local area inside img
%       2. Non-linear least squares minimization for 7 (or 6) - parameter
%           Gaussian function on the ROI selected.
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
% NOTES:
%       This code 'gaussFit.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%
%       For testing purposes, run this script:
%
%       [x,y]=ndgrid(linspace(-1,1,20),linspace(-1,1,20));
%       imwrite(exp(-x.^2/2/2^2-y.^2/2/3^2)+.02*randn(size(x)),'gImg.tif','tif');
%
%       mdl=gaussFit();
%
%

%% declaring fitting predicates

% user specified width in pixels of the 'local area' to be selected from
% the data for fitting. should be odd-valued.
nPixels=11;

% freely rotating bivariate gaussian function for least squares minimization
% parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude, offset]
% xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
% yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);
% f=@(p,X) exp( -xR(X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(4)^2 + ...
%     -yR( X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(5)^2 ) *p(6) + p(7);

% fixed angle fit
% parameters: [xCenter, yCenter, xSD, ySD, amplitude, offset]
th=0;
xR=@(x,y,xc,yc)(x-xc)*cos(th)-(y-yc)*sin(th);
yR=@(x,y,xc,yc)(x-xc)*sin(th)+(y-yc)*cos(th);
f=@(p,X) exp( -xR(X(:,1), X(:,2), p(1), p(2)).^2/2/p(3)^2 + ...
    -yR( X(:,1), X(:,2), p(1), p(2)).^2/2/p(4)^2 ) *p(5) + p(6);

% bounds
lb=[-inf, -inf, 0, 0, -inf, -inf];
ub=[inf, inf, inf, inf, inf, inf];


%% data selection

if findTheSpot
    % select the local area around a bright spot in a larger image
    
    % bandpass and threshold
    LP=1;           % low pass value
    HP=10;          % high pass value
    intThresh=0.1;  % intensity threshold. set to zero and then check by inspection
    hMax=0.1;       % larger if the dynamic range of your data is larger
    lzero=4;        % this squelches a 5 pixel boundary around the filtered image
    
    bIm=bpassDJR(img, LP, HP, intThresh, lzero);
    
    % watershed algorithm
    extImg=imextendedmax(bIm,hMax);
    
    % shrink to a point. this is the estimated location of the spot
    sIm=bwmorph(extImg,'shrink',inf);
    
    % the index of the one pixel is a good guess for the particle location
    [locInds(:,1),locInds(:,2)]=find(sIm);
    
    % temporally coincdident guesses are not treated with this code.
    if size(locInds(:,1))~=1
        fitPars=nan(1,6);
        conf95=nan(1,6);
        return
    end
    
else
    % otherwise, assume the spot is in near the center of the image
    locInds=round(size(img)/2);
end

% pad the img(s) with nans (removed later).
padsize=[nPixels,nPixels];
padVal=nan;
direction='both';
img=padarray(img,padsize,padVal,direction);
locInds=locInds+nPixels;

% find the selection domain
[sDom1,sDom2]=ndgrid(locInds(1)-(nPixels-1)/2:locInds(1)+(nPixels-1)/2, ...
    locInds(2)-(nPixels-1)/2:locInds(2)+(nPixels-1)/2);
inds=sub2ind(size(img),sDom1(:),sDom2(:));

% select the data
truImg=reshape(img(inds),[nPixels,nPixels]);

%% starting parameter selection for 6-parameter Gaussian Fit

% x, y centers starting guess
pStart(1)=0;
pStart(2)=0;

% xSD, ySD in units of pixels
pStart(3)=2;
pStart(4)=2;

% amplitude, offset
mVals=[max(truImg(:)),min(truImg(:))];
pStart(5)=mVals(1)-mVals(2);
pStart(6)=mVals(2);

%% fitting the data

[x,y]=ndgrid(1:nPixels,1:nPixels);
X=cat(2,x(:),y(:)) - nPixels/2;
[fitPars, ~, residual, ~, ~, ~,jacobian] = ...
    lsqcurvefit(f,pStart,X(~isnan(truImg(:)),:),truImg(~isnan(truImg(:))),lb,ub);

% confidence intervals
conf95 = nlparci(fitPars, residual,'jacobian',jacobian);

%% plot the output
fVals=reshape(f(fitPars,X),[nPixels,nPixels]);
dVals=truImg;
sVals=reshape(f(pStart,X),[nPixels,nPixels]);

subplot(221)
title('Data')
pcolor(kron(dVals,ones(10)))
shading flat; axis image; colorbar

subplot(222)
title('starting values')
pcolor(kron(sVals,ones(10)))
shading flat; axis image; colorbar

subplot(223)
title('fit result')
pcolor(kron(fVals,ones(10)))
shading flat; axis image; colorbar

subplot(224)
title('residuals')
pcolor(kron(dVals - fVals,ones(10)))
shading flat; axis image; colorbar

% shift center back to lab frame
fitPars([1,2])=fitPars(1:2)+locInds-nPixels;
end