function [tr,caseInfo]=testCPD(dType)
% use tr as the input into CPD to see if it gives you what you expect.
% 
% dType should either be a number from 1 to 8 or the caseInfo output from a
% previous run. if dType is caseInfo then all the simulation parameters are
% pulled from the dType structure.

% types of diffusion simulatable
diffusionTitles{1}='one mobile no noise.';
diffusionTitles{2}='one mobile with noise.';
diffusionTitles{3}='one mobile one immobile no noise.';
diffusionTitles{4}='one mobile one immobile with noise.';
diffusionTitles{5}='two mobile no noise.';
diffusionTitles{6}='two mobile with noise.';
diffusionTitles{7}='two mobile one immobile no noise.';
diffusionTitles{8}='two mobile one immobile with noise.';

if ~isstruct(dType)
    % if the input is just the diffusion version number
    
    % number of total steps taken
    nSteps=1e4;
    
    % diffusion coefficients in microns^2/s
    d1=1;
    d2=.01;
    
    % fraction of total steps taken by each diffusing population
    d1Frac=1/2;
    d2Frac=1/3;
    
    % camera integration time in seconds
    intTime=.04;
    
    % standard deviation of the noise in pixels
    sdNoise=1;
    
    % magnification of the microscope in microns/pixel
    mag=.049;
    
else
    % if the input is the caseInfo structure from a previous run
    nSteps=dType.nSteps;
    d1=dType.d1;
    d2=dType.d2;
    d1Frac=dType.d1Frac;
    d2Frac=dType.d2Frac;
    intTime=dType.intTime;
    sdNoise=dType.sdNoise;
    mag=dType.mag;
    
    % get the diffusion version number from
    dType=find(strcmp(dType.diffusionTitle,diffusionTitles));
end

switch dType
    case 1      % one mobile no noise
        
        % diffusive trajectory nSteps long with diffusion coefficient d1
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(nSteps,2),1)/mag;
        
    case 2      % one mobile noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(nSteps,2),1)/mag;
        
        % add noise to every track
        tr=cellfun(@(x)x+sdNoise*randn(size(x)),tr,'uniformoutput',0);
        
    case 3      % one mobile one immobile no noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(round(nSteps*d1Frac),2),1)/mag;
        
        % stationary trajectory some fraction of nSteps long
        tr{2}=zeros(round(nSteps*(1-d1Frac)),2);
        
    case 4      % one mobile one immobile noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(round(nSteps*d1Frac),2),1)/mag;
        tr{2}=zeros(round(nSteps*(1-d1Frac)),2);
        tr=cellfun(@(x)x+sdNoise*randn(size(x)),tr,'uniformoutput',0);
        
    case 5      % two mobile no noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(round(nSteps*d1Frac),2),1)/mag;
        tr{2}=cumsum(sqrt(2*d2*intTime)*...
            randn(round(nSteps*(1-d1Frac)),2),1)/mag;
        
    case 6      % two mobile noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(round(nSteps*d1Frac),2),1)/mag;
        tr{2}=cumsum(sqrt(2*d2*intTime)*...
            randn(round(nSteps*(1-d1Frac)),2),1)/mag;
        tr=cellfun(@(x)x+sdNoise*randn(size(x)),tr,'uniformoutput',0);
        
    case 7      % two mobile one immobile no noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(round(nSteps*d1Frac),2),1)/mag;
        tr{2}=cumsum(sqrt(2*d2*intTime)*...
            randn(round(nSteps*d2Frac),2),1)/mag;
        tr{3}=zeros(round(nSteps*(1-d1Frac-d2Frac)),2);
        
        
    case 8      % two mobile one immobile noise
        tr{1}=cumsum(sqrt(2*d1*intTime)*...
            randn(round(nSteps*d1Frac),2),1)/mag;
        tr{2}=cumsum(sqrt(2*d2*intTime)*...
            randn(round(nSteps*d2Frac),2),1)/mag;
        tr{3}=zeros(round(nSteps*(1-d1Frac-d2Frac)),2);
        tr=cellfun(@(x)x+sdNoise*randn(size(x)),tr,'uniformoutput',0);
end

% organize tracks for use in the CPD code
tr=cellfun(@(x,n)organize(x,n),tr,num2cell(1:numel(tr)),'uniformoutput',0);

% % plot trajectories on top of each other
% for ii=1:numel(tr)
%     plot(tr{ii}(:,4),tr{ii}(:,5))
%     hold all
% end
% axis image
% hold off
% title(diffusionTitles{dType})

% concatenate tracks into one matrix
tr=cat(1,tr{:});

% save the diffusion type info
caseInfo.diffusionTitle=diffusionTitles{dType};
caseInfo.nSteps=size(tr,1);
caseInfo.intTime=intTime;
caseInfo.sdNoise=sdNoise;
caseInfo.mag=mag;
caseInfo.d1=d1;
caseInfo.d2=d2;
caseInfo.d1Frac=d1Frac;
caseInfo.d2Frac=d2Frac;
end

function tr=organize(tr,n)
% first column is the track number, second is time, third is useless and
% 4th and 5th are the spatial coordinates of the trajectory

tr=[n*ones(size(tr,1),1),(1:size(tr,1))',(1:size(tr,1))',tr];
end