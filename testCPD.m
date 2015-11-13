function [tr,caseInfo]=testCPD(dType)
% example: x=CPDGlobal(testCPD(2),1,0,'unconfined');
% this code generates trajectories that can be used to test the CPD code
% use tr as the input into CPD to see if it gives you what you expect.
%
% The output trajectory is in units of pixels
% 
% dType should either be a number from 1 to 8 or the caseInfo output from a
% previous run. if dType is caseInfo then all the simulation parameters are
% pulled from the dType structure.

% types of diffusion simulatable
diffusionTitles{1}='one mobile, no noise.';
diffusionTitles{2}='one mobile, with noise.';
diffusionTitles{3}='one mobile, one immobile, no noise.';
diffusionTitles{4}='one mobile, one immobile, with noise.';
diffusionTitles{5}='two mobile, no noise.';
diffusionTitles{6}='two mobile, with noise.';
diffusionTitles{7}='two mobile, one immobile, no noise.';
diffusionTitles{8}='two mobile, one immobile, with noise.';
diffusionTitles{9}='one confined mobile no noise.';
diffusionTitles{10}='one confined mobile, with noise.';
diffusionTitles{11}='one confined mobile, one immobile, no noise.';
diffusionTitles{12}='one confined mobile, one immobile, with noise.';
diffusionTitles{13}='two confined mobile, no noise.';
diffusionTitles{14}='two confined mobile, with noise.';
diffusionTitles{15}='two confined mobile, one immobile, no noise.';
diffusionTitles{16}='two confined mobile, one immobile, with noise.';

if ~isstruct(dType)
    % if the input is just the diffusion version number
    
    % magnification of the microscope in microns/pixel
    mpp=.049;
    
    % number of total steps taken
    nSteps=1e3;
    
    % diffusion coefficients in microns^2/s
    d1=.1;
    d2=.1;
    
    % fraction of total steps taken by each diffusing population
    d1Frac=1/3;
    d2Frac=1/3;
    
    % confinement lengths in each dimension in microns
    cLength1=.5;
    cLength2=1;
    
    % camera integration time in seconds
    intTime=.04;
    
    % standard deviation of the noise in pixels
    sdNoise=.02;
    
else
    % if the input is the caseInfo structure from a previous run
    nSteps=dType.nSteps;
    d1=dType.d1;
    d2=dType.d2;
    d1Frac=dType.d1Frac;
    d2Frac=dType.d2Frac;
    intTime=dType.intTime;
    sdNoise=dType.sdNoise;
    mpp=dType.mpp;
    cLength1=dType.cLength1;
    cLength2=dType.cLength2;

    % get the diffusion version number from
    dType=find(strcmp(dType.diffusionTitle,diffusionTitles));
end

if dType==1||dType==2||dType==9||dType==10          % one mobile
    if dType<3                                      % unconfined
        % diffusive trajectory nSteps long with diffusion coefficient d1
        tr{1}=cumsum(sqrt(2*d1*intTime)*randn(nSteps,2),1);
    else                                            % confined
        tr{1}=confTrack([cLength1,cLength2],d1,intTime,nSteps);
    end
    
elseif dType==3||dType==4||dType==11||dType==12     % one mobile one immobile
    if dType<5                                      % unconfined
        tr{1}=cumsum(sqrt(2*d1*intTime)*randn(round(nSteps*d1Frac),2),1);
    else
        tr{1}=confTrack([cLength1,cLength2],d1,intTime,nSteps);
    end
    
    % stationary trajectory some fraction of nSteps long
    tr{2}=zeros(round(nSteps*(1-d1Frac)),2);
    
elseif dType==5||dType==6||dType==13||dType==14     % two mobile
    if dType<7                                      % unconfined
        tr{1}=cumsum(sqrt(2*d1*intTime)*randn(round(nSteps*d1Frac),2),1);
        tr{2}=cumsum(sqrt(2*d2*intTime)*randn(round(nSteps*(1-d1Frac)),2),1);
    else
        tr{1}=confTrack([cLength1,cLength2],d1,intTime,round(nSteps*d1Frac));
        tr{2}=confTrack([cLength1,cLength2],d2,intTime,round(nSteps*(1-d1Frac)));
    end
    
elseif dType==7||dType==8||dType==15||dType==16     % two mobile one immobile
    if dType<9
        tr{1}=cumsum(sqrt(2*d1*intTime)*randn(round(nSteps*d1Frac),2),1);
        tr{2}=cumsum(sqrt(2*d2*intTime)*randn(round(nSteps*d2Frac),2),1);
    else
        tr{1}=confTrack([cLength1,cLength2],d1,intTime,round(nSteps*d1Frac));
        tr{2}=confTrack([cLength1,cLength2],d2,intTime,round(nSteps*d2Frac));
    end
    
    tr{3}=zeros(round(nSteps*(1-d1Frac-d2Frac)),2);
    
end

% add noise to the even cases
if ~rem(dType,2)
    tr=cellfun(@(x)x+sdNoise*randn(size(x)),tr,'uniformoutput',0);
end

% organize tracks for use in the CPD code
tr=cellfun(@(x,n)organize(x,n),tr,num2cell(1:numel(tr)),'uniformoutput',0);

% plot trajectories on top of each other
% for ii=1:numel(tr)
%     plot(tr{ii}(:,5),tr{ii}(:,4))
%     hold all
% end
% if dType>8
%     line([0,cLength2,cLength2,0,0],[0,0,cLength1,cLength1,0],'color','k');
% end
% hold off
% axis image
% axis([-.1,cLength2+.1,-.1,cLength1+.1])
% title(diffusionTitles{dType})

% concatenate tracks into one matrix
tr=cat(1,tr{:});

% save the diffusion type info
caseInfo.diffusionTitle=diffusionTitles{dType};
caseInfo.nSteps=size(tr,1);
caseInfo.intTime=intTime;
caseInfo.sdNoise=sdNoise;
caseInfo.mag=mpp;
caseInfo.d1=d1;
caseInfo.d2=d2;
caseInfo.d1Frac=d1Frac;
caseInfo.d2Frac=d2Frac;
caseInfo.cLength1=cLength1;
caseInfo.cLength2=cLength2;
end

function tr=organize(tr,n)
% first column is the track number, second is time, third is useless and
% 4th and 5th are the spatial coordinates of the trajectory

tr=[n*ones(size(tr,1),1),(1:size(tr,1))',(1:size(tr,1))',tr];
end

function tr=confTrack(cL,d,t,n)
% generate square-confined track

tr=zeros(n,2);
tr(1,:)=cL/2;
for ii=1:n-1
    step=sqrt(2*d*t)*randn(1,2);
    candpos=tr(ii,:)+step;
    while any(candpos>=cL)||any(candpos<=0)
        w=candpos>=cL;
        candpos(w)=candpos(w)-2*(candpos(w)-cL(w));
        
        x=candpos<=0;
        candpos(x)=-candpos(x);
    end
    tr(ii+1,:)=candpos;
end
end