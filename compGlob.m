% compare all versions of cpd models local vs global

%% data generation parameters

% signal to noise ratio defined as the ratio of the RMS step size to the RMS noise value
SNR = [0,10];

% diffusion coefficients
D = [0,.001,.01,.1,.3,1];

% population ratio multiplicative increment
dMult = 5;

% population ratios
p{1}=[.5;.25];
p{2}=[.3,.3;.2,.2];
p{3}=[.25,.25,.25;.1,.1,.1];

%% analysis parameters

% diffusion dimensionality
d = [1,2];

% number of diffusive populations
nDiffs = [1,2];

% presence or absence of immobile population
immPop = [0,1];

% number of time lags considered
nTau = [3,10];

% one mobile
for hh = 2:numel(D)
    for ii = 1:numel(SNR)
        % generate data
        tr = trackGen('D',D(hh),'weights',1);
        
        % add noise
        nTr = tr + SNR(ii)*randn(size(tr));
        
        for kk = 1:nDiffs
            for mm = 1:numel(immPop)
                % global fit
                resGlo = CPD2(nTr,'globBool',1,'nMobile',kk,'immBool',immPop(mm));
                
                % local fit
                resLoc = CPD2(nTr,'globBool',0,'nMobile',kk,'immBool',immPop(mm));
            end
        end
    end
end

% two mobile
for hh = 1:numel(D)
    for ii = 1:numel(SNR)
        for jj = 1:numel(p{1})
            % generate data
            tr = trackGen('D',[D(hh),D(hh)*dMult],'weights',[p{1}(:,jj),1-p{1}(:,jj)]);
            
            % add noise
            nTr = tr + SNR(ii)*randn(size(tr));
            
            for kk = 1:nDiffs
                for mm = 1:numel(immPop)
                    % global fit
                    resGlo = CPD2(nTr,'globBool',1,'nMobile',jj,'immBool',immPop(kk));
                    
                    % local fit
                    resLoc = CPD2(nTr,'globBool',0,'nMobile',jj,'immBool',immPop(kk));
                end
            end
        end
    end
end

% one mobile, one immobile
for hh = 1:numel(D)
    for ii = 1:numel(SNR)
        for jj = 1:numel(p{2})
            % generate data
            tr = trackGen('D',[0,D(hh)],'weights',[p{2}(:,jj),1-p{2}(:,jj)]);
            
            % add noise
            nTr = tr + SNR(ii)*randn(size(tr));
            
            for kk = 1:nDiffs
                for mm = 1:numel(immPop)
                    % global fit
                    resGlo = CPD2(nTr,'globBool',1,'nMobile',jj,'immBool',immPop(kk));
                    
                    % local fit
                    resLoc = CPD2(nTr,'globBool',0,'nMobile',jj,'immBool',immPop(kk));
                end
            end
        end
    end
end

% two mobile, one immobile
for hh = 1:numel(D)
    for ii = 1:numel(SNR)
        for jj = 1:numel(p{3})
            % generate data
            tr = trackGen('D',[0,D(hh),D(hh)*dMult],'weights',[p{3}(:,jj),1-p{3}(:,jj)]);
            
            % add noise
            nTr = tr + SNR(ii)*randn(size(tr));
            
            for kk = 1:nDiffs
                for mm = 1:numel(immPop)
                    % global fit
                    resGlo = CPD2(nTr,'globBool',1,'nMobile',jj,'immBool',immPop(kk));
                    
                    % local fit
                    resLoc = CPD2(nTr,'globBool',0,'nMobile',jj,'immBool',immPop(kk));
                end
            end
        end
    end
end