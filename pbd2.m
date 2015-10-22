function [cMean,tTrace]=pbd2(v,trFile)
% [cMean,tTrace]=pbd2(v,trFile) v is the 3D double version of the entire
% raw movie data trFile is the tracking file output from masterfit2 with
% dimensions: n x 16

%% parameters
% user supplied values
fitHalfSize=15;
windowHalfSize=10;
usedparameters=12345;

% constants
vSize=size(v);

%% find a pair of nanorods
% user provided guesses for the particle locations
sIm=mean(v,3);
imshow(sIm,[]);

% the positions are flipped for some reason
nrPos=round(fliplr(ginput(2)));

% fit the two nanorods
for ii=1:size(nrPos,1)
    data=sIm(nrPos(ii,1)-fitHalfSize:nrPos(ii,1)+fitHalfSize,...
        nrPos(ii,2)-fitHalfSize:nrPos(ii,2)+fitHalfSize);
    
    pFit(ii,:)=gaussfit(data,usedparameters,0)'+...
        [0,0,nrPos(ii,1),nrPos(ii,2),0];
end

% displacement vector between the two nanorods
rVec=pFit(2,3:4)-pFit(1,3:4);
% rMag=norm(rVec);
% th=acos(rVec(1)/rMag);

%% correlations of intensity through time
% preallocate arrays
trNums=unique(trFile(:,1));
tTrace=cell(max(trNums),2,2);
corr=cell(max(trNums),2);
for ii=trNums'
    % currently considered track
    wTrack=real(trFile(trFile(:,1)==ii,:));
    
    % indices of the local area around the tracked particle in the movie
    pInds(1,:)=mean(wTrack(:,4))-windowHalfSize:...
        mean(wTrack(:,4))+windowHalfSize;
    pInds(2,:)=mean(wTrack(:,5))-windowHalfSize:...
        mean(wTrack(:,5))+windowHalfSize;
    
    % if the pbd vector is negative (?) for the particle in the track
    if all(all(pInds-rVec(ones(1,2*windowHalfSize+1),:)'>=1&...
            pInds-rVec(ones(1,2*windowHalfSize+1),:)'<...
            vSize(ones(1,2*windowHalfSize+1),1:2)'&...
            pInds>=1&pInds<=vSize(ones(1,2*windowHalfSize+1),1:2)'))
        
        % data selections
        v1=v(round(pInds(1,:)-rVec(1)),round(pInds(2,:)-rVec(2)),wTrack(:,2));
        v2=v(round(pInds(1,:)),round(pInds(2,:)),wTrack(:,2));
        
        % save the data
        tTrace(ii,1:2,1)={v1,v2};
        
        % remove the offset of the data and normalize its standard deviation
        v1=bsxfun(@plus,v1,-mean(v1,3)); v1=bsxfun(@times,v1,1./std(v1,0,3));
        v2=bsxfun(@plus,v2,-mean(v2,3)); v2=bsxfun(@times,v2,1./std(v2,0,3));
        
        % compute the convolution of one pixel's time trace with its
        % counterpart in the shifted data selection
        for jj=1:2*windowHalfSize+1
            for kk=1:2*windowHalfSize+1
                corr{ii,1}(jj,kk,:)=conv(squeeze(v1(jj,kk,:)),...
                    squeeze(v2(jj,kk,:)),'same');
            end
        end
    end
    
    % if the pbd vector is positive (?) for the particle in the track
    if all(all(pInds+rVec(ones(1,2*windowHalfSize+1),:)'>=1&...
            pInds+rVec(ones(1,2*windowHalfSize+1),:)'<...
            vSize(ones(1,2*windowHalfSize+1),1:2)'&...
            pInds>=1&pInds<=vSize(ones(1,2*windowHalfSize+1),1:2)'))
        
        % data selections
        v1=v(round(pInds(1,:)+rVec(1)),round(pInds(2,:)+rVec(2)),wTrack(:,2));
        v2=v(round(pInds(1,:)),round(pInds(2,:)),wTrack(:,2));
        
        % save the data
        tTrace(ii,1:2,2)={v1,v2};
        
        % remove the offset of the data and normalize its standard deviation
        v1=bsxfun(@plus,v1,-mean(v1,3)); v1=bsxfun(@times,v1,1./std(v1,0,3));
        v2=bsxfun(@plus,v2,-mean(v2,3)); v2=bsxfun(@times,v2,1./std(v2,0,3));
        
        % compute the convolution of one pixel's time trace with its
        % counterpart in the shifted data selection
        for jj=1:2*windowHalfSize+1
            for kk=1:2*windowHalfSize+1
                corr{ii,2}(jj,kk,:)=conv(squeeze(v1(jj,kk,:)),...
                    squeeze(v2(jj,kk,:)),'same');
            end
        end
    end
end
cMean=cellfun(@(x)nanmean(x(:)),corr);

return
%% put the cursor here and press ctrl-enter after making sure ii exists
ii=ii+1;
subplot(221);
pcolor(mean(tTrace{ii,1,1},3));
title(['1: ' num2str(cSum(ii,1))]); shading flat; axis image;
subplot(222);
pcolor(mean(tTrace{ii,2,1},3));
title('1: original'); shading flat; axis image;
subplot(223);
pcolor(mean(tTrace{ii,1,2},3));
title(['2: ' num2str(cSum(ii,2))]); shading flat; axis image;
subplot(224);
pcolor(mean(tTrace{ii,2,2},3));
title('1: original'); shading flat; axis image
end