function [v, simProps] = dataGen(boundaryCondition, blurFlag)
%
% NAME:
%       dataGen
% PURPOSE:
%       Generates space and time resolved single-molecule imaging data of a
% 			single diffuser.
% CATEGORY:
%       Data Simulation
% CALLING SEQUENCE:
%       [v, simProps] = dataGen(boundaryCondition, blurFlag);
% INPUTS:
%       boundaryCondition:  'confined' or 'unconfined'
% 		blurFlag:			1 or 0, include or exclude blur
% OUTPUTS:
%       v:            		simulated image time sequence
%       simProps:        	properties of the simulation in the Matlab structure format
% PROCEDURE:
%       1. Simulate molecular trajectory
%       2. Evaluate pixel intensities
%		3. Add noise
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
% NOTES:
%       This code 'dataGen.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%
%       For testing purposes, this line makes a movie of confined, blurry
%       diffusion.
%
%       v=dataGen('confined', 1)


%% simulation parameters
D = .1;      			% diffusion coefficient in microns^2/s
tFrame = .05;			% frame integration time in seconds
pixSize = .049;			% width of pixels in microns
psfSize = .098;			% s.d. of the psf in microns
celSize = [1,3];        % [width, length] of confinement cylinder in microns
nFrames = 1e3;			% number of frames in the simulated movie
SNR = 20;				% signal to noise ratio for added white noise (0:inf)

% algorithmically-determined image size designed to disallow edge effects
imSize = ceil(celSize/pixSize+4*ceil(psfSize/pixSize));

% check variable size
if prod([imSize(1:2),nFrames])*8/1e9>1
    warning('video is over a GB. manually pass this block if you wish to continue')
end

% minimum increment of speed in units of microns^2/subframe
dtRef = 0.0001;

% number of subframes required
nSubs = ceil(D*tFrame/dtRef);

% update the value of the diffusion coefficient since rounding may change it.
D = nSubs*dtRef/tFrame;

%% Trajectory generation
if strcmp(boundaryCondition, 'confined')
    
    % confined particle trajectory
    mLocs = zeros(3, nFrames);
    for ii = 1 : nFrames*nSubs-1
        % three 1d steps pulled from normal distribution with variance 2*dtRef
        step = sqrt(2*dtRef) * randn(3,1);
        
        candPos = mLocs(:,ii) + step;
        prevPos = mLocs(:,ii);
        
        % the ordering in the celSize vector matters because of this line:
        r = celSize(1)/2;
        
        % if the candidate position is outside of the cylinder in the x/z dimensions, reflect the step
        % against the inside of the cylinder. path length is preserved.
        if sqrt(sum(candPos([1,3]).^2)) > celSize(1)/2
            m = (candPos(3)-prevPos(3)) / (candPos(1)-prevPos(1));
            b = candPos(3)-m*candPos(1);
            
            xi=[(-m*b+sqrt(-b^2+r^2+m^2*r^2))/(1+m^2),...
                (-m*b-sqrt(-b^2+r^2+m^2*r^2))/(1+m^2)];
            yi=m*xi+b;
            
            % there are two solutions. the one closest to the candidate position is chosen. the farther
            % one is on the other side of the cell.
            whichone = (xi-candPos(1)).^2 + (yi-candPos(3)).^2 + (xi-prevPos(1)).^2 + (yi-prevPos(3)).^2;
            xip = [xi(find(whichone == min(whichone))),yi(find(whichone == min(whichone)))]; %#ok<FNDSB>
            
            normv = -xip/sqrt(sum(xip.^2));
            l = sqrt(sum((candPos([1,3])'-xip).^2));
            pf = 2*sum((prevPos([1,3])'-xip).*normv)*normv-(prevPos([1,3])'-xip);
            pf = pf/sqrt(sum(pf.^2))*l+xip;
            
            % replace x/z components of the position with the reflected x/z components
            out=candPos;
            out([1,3])=pf;
            
            candPos=out;
        end
        
        % if the candidate position is outside of the cylinder in the y dimension (cell's long axis)
        if candPos(2) < -celSize(2)/2
            candPos(2) = 2*-celSize(2)/2 - candPos(2);
        end
        if candPos(2) > celSize(2)/2
            candPos(2) = 2*celSize(2)/2 - candPos(2);
        end
        
        mLocs(:, ii+1) = candPos;
    end
elseif strcmp(boundaryCondition, 'unconfined')
    
    % unconfined particle trajectory
    mLocs=cumsum(sqrt(2*dtRef) * randn(3,nFrames*nSubs),2);
end

%% movie generation

if blurFlag
    % arrange tracks for subframe averaging
    tr_x = zeros(nSubs, nFrames);
    tr_y = zeros(nSubs, nFrames);
    for ii = 1:nFrames
        tr_x(:, ii) = mLocs(1, 1+(ii-1)*nSubs : ii*nSubs);
        tr_y(:, ii) = mLocs(2, 1+(ii-1)*nSubs : ii*nSubs);
    end
elseif ~blurFlag
    % just use the first subframe from each frame
    tr_x = mLocs(1, 1:nSubs:end);
    tr_y = mLocs(2, 1:nSubs:end);
end

% shift to positive values
tr_x=tr_x+celSize(1)/2;
tr_y=tr_y+celSize(2)/2;

% noiseless pixel intensities
v=zeros([imSize(1:2),nFrames]); cx = 0; cy = 0;
for ii=-2*psfSize:pixSize:celSize(1)+2*psfSize        % x pixel locations with padding
    cx = cx+1;
    for jj=-2*psfSize:pixSize:celSize(2)+2*psfSize    % y pixel locations with padding
        cy = cy+1;
        % symmeteric gaussian function approximation of Airy Disk
        v(cx,cy,:) = mean(exp(-((ii-tr_x).^2+(jj-tr_y).^2)/2/psfSize^2),1);
    end
    cy = 0;
end

% add white noise to achieve SNR of 6, defined as the ratio of the maximum intensity
% of a fluorescent spot to the standard deviation of the background noise.
% note: max(v(:)) is 1
v = v + 1/SNR*randn(size(v));

% simulation properties structure
simProps.DiffusionCoefficient = D;
simProps.IntegrationTime = tFrame;
simProps.SNR = SNR;
simProps.PixelWidth = pixSize;
simProps.ConfinementFlag = strcmp(boundaryCondition, 'confined');
simProps.BlurFlag = blurFlag;
simProps.NumSubframes = nSubs;
simProps.dtRef = dtRef;
end