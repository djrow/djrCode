function [mdl,fitPars,fitCI]=gaussFit(img,findTheSpot)
% 
% NAME:
%       gaussFit
% PURPOSE:
%       Fits a generalized gaussian function to 2d imaging data.
% CATEGORY:
%       Image Processing
% CALLING SEQUENCE:
%       [] = gaussFit(inputs,findTheSpot);
% INPUTS:
%       inputs:         The two-dimensional array to be filtered.
%       findTheSpot: (optional)    jklhkljhlh
% OUTPUTS:
%               mdl:    fittin result "object"
%               
% PROCEDURE:
%       1. Data ROI selection of local area 
%       2. Non-linear least squares minimization for 7-parameter Gaussian 
%           function on the ROI selected.
%       3. Find the confidence intervals of all the parameters and various
%           other outputs
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
    filterspec='*.mat';

    % the title parameter is a string containing the title of the dialog box.
    title='Choose the data. You can shift-click to select multiple files.';

    % the file parameter is a string containing the name to use as the default selection.
    file='';

    % filename is a cell array of strings if multiple filenames are selected. Otherwise, 
    % it is a string representing the selected filename.

    % pathname is is a string containing the path of the file selected in the dialog box.
    % If the user presses Cancel, it is set to 0.

    % filterindex is 

    [filename, pathname, filterindex] = uigetfile(filterspec, title, file, 'multiselect','on')
    imLoc={filename}
    for ii=1:numel(filename)
        img{ii}=imread(imLoc{ii});
    end

    if numel(img)==1
        manyImages=0;
    else
        manyImages=1;
    end

    % default behavior: estimate particle location and fit the local area
    findTheSpot=1;
    
elseif nargin>0
    % there are inputs
    thereAreNoInputs=0;

    if ~exist('findTheParticle','var')
        findTheSpot=0;
    end

    
end

%% declaring fitting predicates

if thereAreNoInputs

    % user specified width in pixels of the 'local area' to be selected from the data for
    % fitting
    nPixels=11;






    
    
else
    nPixels=size(img);

    if any(rem(nPixels,2))
        whichEven=find(~rem(nPixels,2));

        % B = PADARRAY(img,padsize,padVal,direction) pads A with values padVal, and in
        % the direction specified by the string direction. By default, direction is 'both'.
        padsize=1;
        padVal=nan;
        direction='post';

        img=padarray(img,padsize,padVal,direction)    
    end
end

[x,y]=ndgrid(linspace(-.5,.5,nPixels),linspace(-.5,.5,nPixels));
X=cat(2,x(:),y(:));

% explain angle fixing


xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);

% if you wish to fix the angle use this portion of commented code:
% th=userSpecifiedValue;
% xR=@(x,y,xc,yc)(x-xc)*cos(th)-(y-yc)*sin(th);
% yR=@(x,y,xc,yc)(x-xc)*sin(th)+(y-yc)*cos(th);

% rotating bivariate gaussian function for least squares minimization
% parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude, offset]
f=@(p,X) exp( xR(X(:,1), X(:,2), p(1), p(2)).^2/2/p(4)^2 + ...
    yR( X(:,1), X(:,2), p(1), p(2)).^2/2/p(5)^2 ) *p(6)+p(7);

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
fitPars=mdl.Coefficients{:,1};

end
% have a nice day