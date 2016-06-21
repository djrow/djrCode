function outStruct = cpdFunFinder2(anProp)

dim = anProp.dim;
nDiffs = anProp.nMobile;
immBool = anProp.immBool;
confBool = anProp.confBool;
globBool = anProp.globBool;

% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% partial 1d cpd function
c1=@(x,y,p)p*erf(sqrt(x./2/y));

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% 1d unconfined msd function
m1=@(t,p)2*p(1)*t+p(2);

pStart = [.1, 0, .05,eps, .01,eps, .0025,eps, ...
    .25, .25, .25];

LB=zeros(1,numel(pStart));
UB=inf(1,numel(pStart));
UB(9:11) = 1;

if dim == 2 && immBool == 0 && confBool == 0 && globBool == 1
    
    switch nDiffs
        case 1
            msdFun=@(tau,p) ...
                m2(tau,p([1,2]));
            cpdFun=@(x,y,p)1-...
                c2(x,y(1),1);
            pID = 1:2;
            
        case 2
            msdFun=@(tau,p)cat(2,...
                m2(tau,p([1,2])),...
                m2(tau,p([3,2])));
            cpdFun=@(x,y,p)1-...
                c2(x,y(1),p(4))-...
                c2(x,y(2),1-p(4));
            pID = [1:3,9];
            
        case 3
            msdFun=@(tau,p)cat(2,...
                m2(tau,p([1,2])),...
                m2(tau,p([3,2])),...
                m2(tau,p([4,2])));
            cpdFun=@(x,y,p)1-...
                c2(x,y(1),p(5))-...
                c2(x,y(2),p(6))-...
                c2(x,y(2),1-p(5)-p(6));
            pID = [1:3,5,9:10];
            
        case 4
            msdFun=@(tau,p)cat(2,...
                m2(tau,p([1,2])),...
                m2(tau,p([3,2])),...
                m2(tau,p([4,2])),...
                m2(tau,p([5,2])));
            cpdFun=@(x,y,p)1-...
                c2(x,y(1),p(6))-...
                c2(x,y(2),p(7))-...
                c2(x,y(3),p(8))-...
                c2(x,y(4),1-p(6)-p(7)-p(8));
            pID = [1:3,5,7,9:11];
            
    end
    
    pStart={pStart(pID)};
    bounds=[{LB(pID)},{UB(pID)}];
    dID = find(ismember(pID,1:2:7));
    aID = find(ismember(pID,9:11));
    
    outStruct.cpdFun = cpdFun;
    outStruct.msdFun = msdFun;
    outStruct.pStart = pStart;
    outStruct.bounds = bounds;
    outStruct.dID = dID;
    outStruct.aID = aID;
else
    warning('unsupported parameters')
    outStruct = [];
end


end