function [mdl,fitPars,fitCI]=gaussFit(img,findTheSpot)
% 
% NAME:
%       gaussFit
% PURPOSE:
%       Fits a generalized gaussian function to 2d imaging data.
% CATEGORY:
%       Image Processing
% CALLING SEQUENCE:
%       [mdl, fitPars, fitCI] = gaussFit(inputs,findTheSpot);
% INPUTS:
%       img:            The two-dimensional array to be fit to a gaussian
%       findTheSpot:    (optional) 1 or 0. Default behavior is to assume the spot is near the 
%                           enter of img
% OUTPUTS:
%       mdl:            fitting result "object"
%       fitPars:        fitting coefficient vector
%       fitCI:          95% confidence interval of fitting coefficients at end of fitting
% PROCEDURE:
%       1. Peak guessing and/or data ROI selection of local area inside img
%       2. Non-linear least squares minimization for 7 (or 6) - parameter Gaussian 
%           function to the ROI selected.
%       3. Find the confidence intervals of parameters.
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
% NOTES:
%       This code 'gaussFit.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.  

%% parsing inputs
if ~nargin
    % there are no inputs
    thereAreNoInputs=1;    

    % user input file location

    % The filterspec parameter determines the initial display of files in the dialog box.
    filterspec{:,1}={'*.mat'; '*.tif*'; '*.avi'};
    filterspec{:,2}={'Matlab data files (*.mat)'; ...
                        'tiff images or tiff image stacks (*.tif or *.tiff)'; ...
                        'avi movie files (*.avi)'}

    % the title parameter is a string containing the title of the dialog box.
    title='Choose the data. You can shift-click to select multiple files.';

    % the file parameter is a string containing the name to use as the default selection.
    file='*.mat';

    % filename is a cell array of strings if multiple filenames are selected. Otherwise, it is a 
    %       string representing the selected filename. 
    % pathname is is a string containing the path of the file selected in the dialog box. If the 
    %       user presses Cancel, it is set to 0.
    % filterindex is the index of the filter selected in the dialog box. The indexing starts at 1. 
    %       If the user presses Cancel, it is set to 0.
    [filename, pathname, filterindex] = uigetfile(filterspec, title, file, 'multiselect','on')
    imLoc={filename}
    for ii=1:numel(filename)
        img{ii}=imread(imLoc{ii});      % load data
    end

    if numel(img)==1
        manyImages=0;
    else
        manyImages=1;
    end

    % default behavior: assume the spot is near the center of img
    findTheSpot=0;
    
elseif nargin>0
    % there are inputs
    thereAreNoInputs=0;

    if ~exist('findTheParticle','var')
        findTheSpot=0;
    end
end

%% declaring fitting predicates

% user specified width in pixels of the 'local area' to be selected from the data for fitting
nPixels=11;

[x,y]=ndgrid(linspace(-.5,.5,nPixels),linspace(-.5,.5,nPixels));
X=cat(2,x(:),y(:));

% rotate the coordinates
xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);

% if you wish to fix the angle use this portion of commented code:
% userSpecifiedValue=[...];
% th=userSpecifiedValue;
% xR=@(x,y,xc,yc)(x-xc)*cos(th)-(y-yc)*sin(th);
% yR=@(x,y,xc,yc)(x-xc)*sin(th)+(y-yc)*cos(th);
% f=@(p,X) exp( xR(X(:,1), X(:,2), p(1), p(2)).^2/2/p(3)^2 + ...
%       yR( X(:,1), X(:,2), p(1), p(2)).^2/2/p(4)^2 ) *p(5) + p(6);

% rotating bivariate gaussian function for least squares minimization
% parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude, offset]
f=@(p,X) exp( xR(X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(4)^2 + ...
    yR( X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(5)^2 ) *p(6) + p(7);

%% data selection

% pad the img(s) with nans (to be removed later).
% B = padarray(img,padsize,padVal,direction) pads img with values padVal, and in the direction 
% specified by the string direction. By default, direction is 'both'.
padsize=5;
padVal=nan;
direction='both';
oImSize=size(img);
img=padarray(img,padsize,padVal,direction)  

if rem(size(img,1),2)>0
    img=cat(1, img, nans(1,size(img,2))); 
end
if rem(size(img,2),2)>0
    img=cat(2, img, nans(1,size(img,1))'); 
end
nImSize=size(img);

if findTheSpot
    % select the local area around a bright spot in a larger image
    
    % bandpass and threshold
    LP=1;           % low pass value
    HP=10;          % high pass value
    T=[];
    hMax=[1e2];     % larger if the dynamic range of your data is larger
    lzero=5;        % this squelches a 5 pixel boundary around the filtered image

    bIm=bpassDJR(img, LP, HP, T, lzero);
    
    % watershed algorithm
    extImg=imextendedmax(bimg,hMax);
    
    % shrink to a point. this is the estimated location of the spot
    sIm=bwmorph(extimg,'shrink',inf);
    
    % the index of the one pixel is a good guess for the particle location
    [locs(:,1),locs(:,2)]=find(bp_img);

    if size(locs(1))>1
        display('fix this')
        return
    end

    [temp1,temp2]=ndgrid(locs(1)-(nPixels-1)/2:locs(1)+(nPixels-1)/2, ...
        locs(2)-(nPixels-1)/2:locs(2)+(nPixels-1)/2);
    
    % select the data
    truImg=reshape(img(temp1(:),temp2(:)),[nPixels,nPixels]);
    
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

mdl=fitnlm(X(~isnan(truImg),truImg(~isnan(truImg)),f,pStart));

%% organizing outputs
% confidence interval
fitCI=diff(coefCI(mdl),1,2);

% fitting coefficients
fitPars=mdl.Coefficients{:,1};

% shift results back to laboratory frame.

[...]


end