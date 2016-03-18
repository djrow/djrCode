function [mdl,fitPars,fitCI,f]=gaussFit(img,findTheSpot)
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
%       [mdl, fitPars, fitCI] = gaussFit(img,findTheSpot);
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

%% parsing inputs
if ~nargin    
    % user input file location
    
    % The filterspec parameter determines the initial display of files in
    % the dialog box.
%     filterspec{:,1}={'*.mat', '*.tiff*', '*.avi'};
    % filterspec{:,2}={'Matlab data files (*.mat)', ...
    %     'tiff images or tiff image stacks (*.tif or *.tiff)', ...
    %     'avi movie files (*.avi)'};
    
    % the title parameter is a string containing the title of the dialog
    % box.
    title='Choose the data.';
    
    % the file parameter is a string containing the name to use as the
    % default selection.
    % file='*.mat';
    
    % filename is a cell array of strings if multiple filenames are
    %       selected. Otherwise, it is a string representing the selected
    %       filename.
    % pathname is is a string containing the path of the file selected in
    %       the dialog box. If the user presses Cancel, it is set to 0.
    % filterindex is the index of the filter selected in the dialog box.
    %       The indexing starts at 1. If the user presses Cancel, it is set to 0.
    % multiselect has been turned off because that capability has not yet
    %       been built into the code.
    [filename, pathname, filterindex]=uigetfile(...
        {'*.tif*','*.tif Images';'*.mat*','*.mat Matlab Data Files'}, title);
    
    img=imread(fullfile(pathname,filename),'tif');
%     for ii=1:numel(imLoc)
        % load data
%         img{ii}=imread(fullfile(pathname,imLoc{ii}));
%     end
    
    % if numel(img)==1
    manyImages=0;
    % else
    % manyImages=1;
    % end
    
    % default behavior: assume the spot is near the center of img
    findTheSpot=0;
    
elseif nargin>0
    % default behavior: assume the spot is near the center of img
    if ~exist('findTheSpot','var')
        findTheSpot=0;
    end
end

%% declaring fitting predicates

% user specified width in pixels of the 'local area' to be selected from
% the data for fitting. should be odd. just keep it odd. i'm superstitious.
nPixels=11;

% coordinate rotation functions
xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);

% if you wish to fix the angle use this portion of commented code:
% userSpecifiedValue=[...]; th=userSpecifiedValue;
% xR=@(x,y,xc,yc)(x-xc)*cos(th)-(y-yc)*sin(th);
% yR=@(x,y,xc,yc)(x-xc)*sin(th)+(y-yc)*cos(th); f=@(p,X) exp( xR(X(:,1),
% X(:,2), p(1), p(2)).^2/2/p(3)^2 + ...
%       yR( X(:,1), X(:,2), p(1), p(2)).^2/2/p(4)^2 ) *p(5) + p(6);

% freely rotating bivariate gaussian function for least squares
% minimization parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude,
% offset]
f=@(p,X) exp( -xR(X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(4)^2 + ...
    -yR( X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(5)^2 ) *p(6) + p(7);

%% data selection

% pad the img(s) with nans (to be removed later). 
% B = padarray(img,padsize,padVal,direction) pads img with values padVal, and
% in the direction specified by the string direction. By default, direction
% is 'both'.
padsize=[nPixels,nPixels];
padVal=nan;
direction='both';
oImSize=size(img);
img=padarray(img,padsize,padVal,direction);

if findTheSpot
    % select the local area around a bright spot in a larger image
    
    % bandpass and threshold
    LP=1;           % low pass value
    HP=10;          % high pass value
    intThresh=100;  % intensity threshold
    hMax=[1e2];     % larger if the dynamic range of your data is larger
    lzero=5;        % this squelches a 5 pixel boundary around the filtered image
    
    bIm=bpassDJR(img, LP, HP, intThresh, lzero);
    
    % watershed algorithm
    extImg=imextendedmax(bIm,hMax);
    
    % shrink to a point. this is the estimated location of the spot
    sIm=bwmorph(extImg,'shrink',inf);
    
    % the index of the one pixel is a good guess for the particle location
    [locInds(:,1),locInds(:,2)]=find(sIm);
    
    % temporally coincdident guesses are not treated with this code.
    if size(locs(:,1))>1
        display('[...fix this...] data quality too poor for this algorithm.')
        return
    end
    
else
    % otherwise, the 
    locs=round(size(img)/2);
end

% find selection domains
[sDom1,sDom2]=ndgrid(locs(1)-(nPixels-1)/2:locs(1)+(nPixels-1)/2, ...
    locs(2)-(nPixels-1)/2:locs(2)+(nPixels-1)/2);

% select the data
truImg=reshape(img(sDom1(:),sDom2(:)),[nPixels,nPixels]);

%% starting parameter selection

% x, y centers starting guess
pStart(1)=0;
pStart(2)=0;

% angle
pStart(3)=0;

% xSD, ySD in units of percent of width of ROI to be fit
pStart(4)=.2;
pStart(5)=.2;

% amplitude, offset
mVals=[max(truImg(:)),min(truImg(:))];
pStart(6)=mVals(1)-mVals(2);
pStart(7)=mVals(2);

%% fitting the data

[x,y]=ndgrid(linspace(-.5,.5,size(truImg,1)),linspace(-.5,.5,size(truImg,1)));
X=cat(2,x(:),y(:));
% mdl=fitnlm(X(~isnan(truImg),:),truImg(~isnan(truImg)),f,pStart);
x=lsqnonlin(@(p)nansum((img-f(p,X)).^2),pStart,lb,ub)

%% organizing outputs
% confidence interval
% fitCI=diff(coefCI(mdl),1,2);

% fitting coefficients
% fitPars=mdl.Coefficients{:,1};

% shift and scale results back to laboratory frame.
fitPars([1,2])=fitPars([1,2]).*nPixels+locs';
fitPars([4,5])=fitPars([4,5]).*nPixels;


%% plot the output




end