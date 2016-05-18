function [fitPars, conf95, g, outPut]=gaussFit(img, varargin)
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
%
%       varargin:       use paired inputs to set the property (input 1) to the
%           value (input 2) desired.
%
%       Properties:     Descriptions:
%
%       findTheSpot:    1 or 0. Default behavior is to fit an
%           ROI in the center of the image. If the spot is not near the
%           center or the image is very large, findTheSpot enables the code
%           to first roughly locate the spot and then use that location as
%           the ROI center.
%
%       plottingFlag:   1 or 0. show plotting output. default is 0.
%
%       widthGuess:     set the starting value for the width of the
%                       Gaussian in units of pixels.
%
%       nPixels         pixel width of ROI to be selected from img. default
%                       is 11. the value should be odd.
%
% OUTPUTS:
%       fitPars:        fitting coefficient vector, all units are pixels.
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
%       img = exp(-x.^2/2/2^2-y.^2/2/3^2)+.02*randn(size(x));
%       p = gaussFit(img,'widthGuess',2);
opts = optimset('Display','off');
warning('off','all')

% default parameters
params.spotSizeLB = 1;
params.spotSizeUB = 10;
params.intThresh = 1;
params.lZero = 4;
params.hMax = 1;
params.searchBool = 1;
params.showFitting = 0;
params.showGuessing = 0;
params.widthGuess = 2;
params.nPixels = 11;        % should be odd valued
params.frameNumber = 1;
params.ffSwitch = 3;        % 3 is a symmetric gaussian (5 parameters) fit
                            % 2 is a fixed angle asymmetric gaussian fit
                            % 1 is a 7 parameter asymmetric gaussian fit
fNames=fieldnames(params);

% if any sim parameters are included as inputs, change the simulation
% parameters mentioned
if nargin>1
    for ii=1:2:nargin-2
        whichField = strcmp(fNames,varargin{ii});
        
        if all(~whichField)
            warning('Check spelling. Parameter change may have not occurred.')
        end
        
        eval(['params.' fNames{whichField} ' = varargin{ii+1};'])
    end
end

switch params.ffSwitch
    case 1
        % freely rotating bivariate gaussian function for least squares minimization
        % parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude, offset]
        xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
        yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);
        f=@(p,X) exp( -xR(X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(4)^2 + ...
            -yR( X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(5)^2 ) *p(6) + p(7);
        
    case 2
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
        
    case 3
        % symmetric gaussian
        f=@(p,X) exp( -((X(:,1)-p(1)).^2 + (X(:,2)-p(2)).^2)/2/p(3).^2) * p(4) + p(5);
        lb = [-params.nPixels, -params.nPixels, 0, 0, -inf];
        ub = [params.nPixels*2, params.nPixels, inf, inf, inf];
        
end

%% rough localization of molecules
if params.searchBool
    % band pass
    bIm=bpassDJR(img, params.spotSizeLB, params.spotSizeUB, params.intThresh, params.lZero);
    
    % watershed algorithm
    extImg=imextendedmax(bIm,params.hMax);
    
    % failed watershed algorithm can result in all ones
    if all(extImg(:))
        extImg=extImg-1;
    end
    
    % shrink to a point. this is the estimated location of the spot
    sIm=bwmorph(extImg,'shrink',inf);
    
    % if shrinking the image produces rings, remove the rings
    cc = bwconncomp(sIm);
    if cc.NumObjects<sum(sIm(:))
        whichBad = cellfun(@numel,cc.PixelIdxList) > 1;
        sIm(cc.PixelIdxList{whichBad}) = 0;
    end
    
    if params.showGuessing
        subplot(221)
        imshow(img,[])
        
        subplot(222)
        imshow(bIm,[])
        
        subplot(223)
        imshow(extImg,[])
        
        subplot(224)
        imshow(sIm,[])
        
        set(gcf,'NextPlot','add');
        axes;
        h = title(['Frame number ' num2str(params.frameNumber)]);
        set(gca,'Visible','off');
        set(h,'Visible','on');
    end
    
    % the index of the one pixel is a good guess for the particle location
    [locInds(:,1),locInds(:,2)]=find(sIm);
    g=locInds;
else
    % otherwise, assume the spot is in near the center of the image
    locInds=round(size(img)/2);
end

nFits = size(locInds,1);
if nFits > 50
    warning(['too many fits in frame number ' num2str(params.frameNumber)])
    fitPars=[];
    conf95=[];
    outPut = [];
    return
end

% skip the fitting if the guessing parameters are being checked
if params.showGuessing
    nFits = 0;
end

%% fit the data
% pad the img(s) with nans (removed later).
padsize=params.nPixels([1,1]);
img=padarray(img,padsize,nan,'both');
locInds=locInds+params.nPixels;

% starting guesses
pStart(1)=0;
pStart(2)=0;
pStart(3)=params.widthGuess;

fitPars = zeros(nFits,numel(lb));
conf95 = zeros(nFits,numel(lb));

for ii = 1:nFits
    % find the selection domain
    [sDom1,sDom2]=ndgrid(locInds(ii,1)-(params.nPixels-1)/2:locInds(ii,1)+(params.nPixels-1)/2, ...
        locInds(ii,2)-(params.nPixels-1)/2:locInds(ii,2)+(params.nPixels-1)/2);
    inds=sub2ind(size(img),sDom1(:),sDom2(:));
    
    % select the data
    truImg=double(reshape(img(inds),[params.nPixels,params.nPixels]));
    
    % amplitude, offset
    mVals=[max(truImg(:)),min(truImg(:))];
    pStart(4)=mVals(1)-mVals(2);
    pStart(5)=mVals(2);
    
    % fitting the data
    [x,y]=ndgrid(1:params.nPixels,1:params.nPixels);
    X=cat(2,x(:),y(:)) - params.nPixels/2;
    [fitPars(ii,:),~,residual,~,outPut(ii),~,jacobian] = ...
        lsqcurvefit(f,pStart,X(~isnan(truImg(:)),:),truImg(~isnan(truImg(:))),lb,ub,opts);
    
    % confidence intervals
    conf95(ii,:) = diff(nlparci(fitPars(ii,:), residual,'jacobian',jacobian),1,2);
    
    % plot the output
    if params.showFitting
        fVals=reshape(f(fitPars(ii,:),X),[params.nPixels,params.nPixels]);
        dVals=truImg;
        sVals=reshape(f(pStart,X),[params.nPixels,params.nPixels]);
        
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
    end
end

if nFits>0
    % shift center back to lab frame
    fitPars(:,[1,2])=fitPars(:,1:2)+locInds-params.nPixels;
end
if ~exist('outPut','var')
    outPut = [];
end
end