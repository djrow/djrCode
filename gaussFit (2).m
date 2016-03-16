function gaussFit(inputs,findTheSpot)
% 
% NAME:
%       gaussFit
% PURPOSE:
%       Fits a generalized gaussian function to 2d imaging data.
% 
% CATEGORY:
%       Image Processing
% CALLING SEQUENCE:
%       [] = gaussFit(inputs,findTheSpot);
% INPUTS:
%       inputs:         The two-dimensional array to be filtered.
%       findTheSpot: (optional)    Characteristic lengthscale of noise in pixels.
%                       Additive noise averaged over this length should
%                       vanish. May assume any positive floating value.
%                       May be set to 0 or false, in which case only the
%                       highpass "background subtraction" operation is 
%                       performed.
%
% OUTPUTS:
%               res:    filtered image.
% PROCEDURE:
%       1. Data ROI selection of local area 
%       2. Non-linear least squares minimization for 7-parameter Gaussian 
%           function on the ROI selected.
%       3. Find the confidence intervals of all the parameters and various
%           other outputs
% NOTES:
% Performs a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 
% original - lowpassed version of original => highpassed version of the
% original
% 
% Performing a lowpass and a highpass results in a bandpassed image.
% 
% Converts input to double.  Be advised that commands like 'image' display 
% double precision arrays differently from UINT8 arrays.

% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
%
%       reatly revised version DGG 5/95.
%
%       Added /field keyword JCC 12/95.
% 
%       Memory optimizations and fixed normalization, DGG 8/99.
%       Converted to Matlab by D.Blair 4/2004-ish
%
%       Fixed some bugs with conv2 to make sure the edges are
%       removed D.B. 6/05
%
%       Removed inadvertent image shift ERD 6/05
% 
%       Added threshold to output.  Now sets all pixels with
%       negative values equal to zero.  Gets rid of ringing which
%       was destroying sub-pixel accuracy, unless window size in
%       cntrd was picked perfectly.  Now centrd gets sub-pixel
%       accuracy much more robustly ERD 8/24/05
%
%
%       This code 'gaussFit.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.  

%% parsing inputs
if ~nargin
    % there are no inputs
    thereAreNoInputs=1;    
    
elseif nargin
    % there are inputs
    thereAreNoInputs=0;    
end

%% loading data

if thereAreNoInputs
    % user input file location
    imLoc=uigetfile();
    
    % load image
    img=imread(imLoc);
    
    % default behavior: estimate particle location and fit the local area
    findTheSpot=1;
    
else
    % the data is an input
    img=inputs;
    
    % default behavior: assume the gaussian shape in the 2d image data is
    % centered and has a width between 1 and 5 pixels.
    if ~exist('findTheParticle','var')
        findTheSpot=0;
    end
end

%% declaring fitting predicates

% width in pixels of the 'local area' to be selected from the data
if thereAreNoInputs
    if round(
    
    
    end
    
else
    nPixels=11;
end

[x,y]=ndgrid(linspace(-.5,.5,nPixels),linspace(-.5,.5,nPixels));
X=cat(2,x(:),y(:));

xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);

% rotating bivariate gaussian function for least squares minimization
% parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude, offset]
f=@(p,X) exp( xR(X(:,1), X(:,2), p(1), p(2), p(3) ).^2/2/p(4)^2 + ...
    yR( X(:,1), X(:,2), p(1), p(2), p(3) ).^2/2/p(5)^2 ) *p(6)+p(7);

%% data selection

% pad the data with nans to fix discretization issues
whichN=rem(size(img),2);
if whichN(1)
    img=padarray;
end
if whichN(2)
    img=padarray;
end


if findTheSpot
    % select the local area around a bright spot in a larger image
    
    % bandpass and threshold
    bIm=bpass();
    
    % watershed algorithm
    wIm=waterShed(bIm);
    
    % shrink to a point. this is the estimated location of the spot
    sIm=imshrink(wIm);
    
    p=find(sIm);
    truImgInds=ndgrid(round(p(1))-(nPixels-1)/2:round(p(1))+(nPixels-1)/2, ...
        round(p(2))-(nPixels-1)/2:round(p(2))+(nPixels-1)/2);
    
    % select the data
    truImg=img(truImgInds(:));
    
else
    % otherwise, the original data is used
    truImg=img;
end

%% starting parameter selection

% x, y centers
pStart(1)=0;
pStart(2)=0;

% angle
pStart(3)=0;

% xSD, ySD
pStart(4)=2;
pStart(5)=2;

% amplitude, offset
mVals=[max(truImg(:)),min(truImg(:))];
pStart(6)=mVals(1)-mVals(2);
pStart(7)=mVals(2);


%% fitting the data

mdl=fitnlm(X,truImg(:),f,pStart);

%% organizing outputs
% confidence interval
fitCI=diff(coefCI(mdl),1,2);

% fitting coefficients
p=mdl.Coefficients{:,1};

end
% have a nice day