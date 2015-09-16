function tr=testCPD(dVer)
% use tr as the input into CPD to see if it gives you what you expect

nFrames=1e3;
d1=1;          % microns^2/s
d2=.01;        % microns^2/s
intTime=.04;    % seconds
sdNoise=.01;      % pixels
mag=.049;       % microns/pixel

switch dVer
    case 1      % one component no noise
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames,2),1)/mag;               % 1 mobile
        tr=organize(tr);
        
    case 2      % one component noise
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames,2),1)/mag;               % 1 mobile
        tr=tr+sdNoise*randn(nFrames,2);                                     % add noise
        tr=organize(tr);
        
    case 3      % one mobile one immobile no noise
        p1Frac=round(nFrames*.1)/nFrames;
        
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames*p1Frac,2),1)/mag;        % 1 mobile
        tr(nFrames*p1Frac+1:nFrames,:)=...
            tr(nFrames*p1Frac*ones(1,round(nFrames*(1-p1Frac))),:);         % 1 immobile
        tr=organize(tr);
        
    case 4      % one mobile one immobile noise
        p1Frac=round(nFrames*.01)/nFrames;
        
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames*p1Frac,2),1)/mag;        % 1 mobile
        tr(nFrames*p1Frac+1:nFrames,:)=...
            tr(nFrames*p1Frac*ones(1,round(nFrames*(1-p1Frac))),:);         % 1 immobile
        tr=tr+sdNoise*randn(nFrames,2);                                     % add noise
        tr=organize(tr);
        
    case 5      % two mobile no noise
        p1Frac=round(nFrames*.999)/nFrames;
        
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames*p1Frac,2),1)/mag;        % 1 mobile
        tr=[tr;cumsum([tr(end,:);sqrt(2*d2*intTime)*...
            randn(round(nFrames*(1-p1Frac)-1),2)],1)/mag];                  % 2 mobile
        tr=organize(tr);
        
    case 6      % two mobile noise
        p1Frac=round(nFrames*.9)/nFrames;
        
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames*p1Frac,2),1)/mag;        % 1 mobile
        tr=[tr;cumsum([tr(end,:);sqrt(2*d2*intTime)*...
            randn(round(nFrames*(1-p1Frac)-1),2)],1)/mag];                  % 2 mobile
        tr=tr+sdNoise*randn(size(tr,1),2);                                  % add noise
        tr=organize(tr);
        
    case 7      % two mobile one immobile no noise
        p1Frac=round(nFrames*.4)/nFrames;
        p2Frac=round(nFrames*.4)/nFrames;
        
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames*p1Frac,2),1)/mag;        % 1 mobile
        tr=[tr;cumsum([tr(end,:);sqrt(2*d2*intTime)*...
            randn(round(nFrames*p2Frac-1),2)],1)/mag];                      % 2 mobile
        tr(nFrames*(p1Frac+p2Frac)+1:nFrames,:)=...
            tr(nFrames*(p1Frac+p2Frac)*...
            ones(1,round(nFrames*(1-p1Frac-p2Frac))),:);                    % 1 immobile
        tr=organize(tr);
        
        
    case 8      % two mobile one immobile noise
        p1Frac=round(nFrames*.4)/nFrames;
        p2Frac=round(nFrames*.4)/nFrames;
        
        tr=cumsum(sqrt(2*d1*intTime)*randn(nFrames*p1Frac,2),1)/mag;        % 1 mobile
        tr=[tr;cumsum([tr(end,:);sqrt(2*d2*intTime)*...
            randn(round(nFrames*p2Frac-1),2)],1)/mag];                      % 2 mobile
        tr(nFrames*(p1Frac+p2Frac)+1:nFrames,:)=...
            tr(nFrames*(p1Frac+p2Frac)*...
            ones(1,round(nFrames*(1-p1Frac-p2Frac))),:);                    % 1 immobile
        tr=tr+sdNoise*randn(size(tr,1),2);                                  % add noise
        tr=organize(tr);
        
end

plot(tr(:,4),tr(:,5)); axis image
end

function tr=organize(tr)
tr=[ones(size(tr,1),1),(1:size(tr,1))',(1:size(tr,1))',tr];
end