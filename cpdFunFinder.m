function outStruct = cpdFunFinder(anProp)

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

pStart = [.1, .5, .0016, ...       % msd 1 parameters D1, L, S1
    .01, .0016, ...                % msd 2 parameters D2, S2
    .001, .0016, ...               % msd 3 parameters D3, S3
    .45, .45, .3, 0];            % cpd function parameters amp1, amp2, amp3, immSize

switch globBool
    case 1
        %% global fitting
        LB=zeros(1,numel(pStart));
        UB=inf(1,numel(pStart));
        UB([8:10]) = 1;
        switch dim
            case 2
                switch immBool
                    case 0 % no immobile term
                        switch nDiffs
                            case 1 % 1 diffuser
                                switch confBool
                                    case 0 % 1 unconfined diffuser in 2D
                                        msdFun=@(tau,p)...
                                            m2(tau,p(1:2));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),1);
                                        pID=[1,3];
                                        
                                    case 1 % 1 confined diffuser in 2D
                                        msdFun=@(tau,p)...
                                            confMSD2(tau,p(1:3));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),1);
                                        pID=1:3;
                                end
                            case 2 % 2 diffusers
                                switch confBool
                                    case 0 % 2 unconfined diffusers in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            m2(tau,p(1:2)),...
                                            m2(tau,p(3:4)));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),p(5))-...
                                            c2(x,y(2),1-p(5));
                                        pID=[1,3,4,5,8];
                                        
                                    case 1 % 2 confined diffusers in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            confMSD2(tau,p(1:3)),...
                                            confMSD2(tau,p([4,2,5])));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),p(6))-...
                                            c2(x,y(2),1-p(6));
                                        pID=[1:5,8];
                                end
                            case 3
                                switch confBool
                                    case 0 % 3 unconfined diffusers in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            m2(tau,p(1:2)),...
                                            m2(tau,p(3:4)),...
                                            m2(tau,p(5:6)));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),p(7))-...
                                            c2(x,y(2),p(8))-...
                                            c2(x,y(3),1-p(7)-p(8));
                                        pID=[1,3:7,8:9];
                                    case 1 % 3 confined diffusers in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            confMSD2(tau,p(1:3)),...
                                            confMSD2(tau,p([4,2,5])),...
                                            confMSD2(tau,p([6,2,7])));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),p(8))-...
                                            c2(x,y(2),p(9))-...
                                            c2(x,y(3),1-p(8)-p(9));
                                        pID=[1:7,8:9];
                                end
                                
                            case 4
                                switch confBool
                                    case 0
                                        msdFun=@(tau,p)cat(2,...
                                            m2(tau,p(1:2)),...
                                            m2(tau,p(3:4)),...
                                            m2(tau,p(5:6)),...
                                            m2(tau,p(7:8)));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1),p(9))-...
                                            c2(x,y(2),p(10))-...
                                            c2(x,y(3),p(11))-...
                                            c2(x,y(4),1-p(9)-p(10)-p(11));
%                                         pID=[1,3:7,,8:9];
                                    case 1
                                        
                                end
                        end
                    case 1 % 1 immobile term
                        switch nDiffs
                            case 1 % 1 diffuser
                                switch confBool
                                    case 0 % 1 unconfined diffuser with 1 imm in 2D
                                        msdFun=@(tau,p)...
                                            m2(tau,p(1:2));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y+p(4),p(3))-...
                                            c2(x,p(4),1-p(3));
                                        pID=[1,3,8,11];
                                    case 1 % 1 confined diffuser with 1 imm in 2D
                                        msdFun=@(tau,p)...
                                            confMSD2(tau,p(1:3));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y+p(5),p(4))-...
                                            c2(x,p(5),1-p(4));
                                        pID=[1:3,8,11];
                                end
                            case 2 % 2 diffusers
                                switch confBool
                                    case 0 % 2 unconfined diffusers with 1 immobile population in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            m2(tau,p(1:2)),...
                                            m2(tau,p(3:4)));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1)+p(7),p(5))-...
                                            c2(x,y(2)+p(7),p(6))-...
                                            c2(x,p(7),1-p(5)-p(6));
                                        pID=[1,3:5,8,9,11];
                                    case 1 % 2 confined diffusers with 1 immobile population in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            confMSD2(tau,p(1:3)),...
                                            confMSD2(tau,p([4,2,5])));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1)+p(8),p(6))-...
                                            c2(x,y(2)+p(8),p(7))-...
                                            c2(x,p(8),1-p(6)-p(7));
                                        pID=[1:5,8,9,11];
                                end
                            case 3
                                switch confBool
                                    case 0 % 3 unconfined diffusers with 1 immobile population in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            m2(tau,p(1:2)),...
                                            m2(tau,p(3:4)),...
                                            m2(tau,p(5:6)));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1)+p(10),p(7))-...
                                            c2(x,y(2)+p(10),p(8))-...
                                            c2(x,y(3)+p(10),p(9))-...
                                            c2(x,p(10),1-p(7)-p(8)-p(9));
                                        pID=[1,3:7,8:11];
                                    case 1 % 3 confined diffusers with 1 immobile population in 2D
                                        msdFun=@(tau,p)cat(2,...
                                            confMSD2(tau,p(1:3)),...
                                            confMSD2(tau,p([4,2,5])),...
                                            confMSD2(tau,p([6,2,7])));
                                        cpdFun=@(x,y,p)1-...
                                            c2(x,y(1)+p(11),p(8))-...
                                            c2(x,y(2)+p(11),p(9))-...
                                            c2(x,y(3)+p(11),p(10))-...
                                            c2(x,p(11),1-p(8)-p(9)-p(10));
                                        pID=[1:5,8:11];
                                end
                        end
                end
            case 1
                switch immBool
                    case 0 % no immobile term
                        switch nDiffs
                            case 1 % 1 diffuser
                                switch confBool
                                    case 0 % 1 unconfined diffuser in 1D
                                        msdFun=@(tau,p)...
                                            m1(tau,p(1:2));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1),1);
                                        pID=[1,3];
                                    case 1 % 1 confined diffuser in 1D
                                        msdFun=@(tau,p)...
                                            confMSD1(tau,p(1:3));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1),1);
                                        pID=1:3;
                                end
                            case 2 % 2 diffusers
                                switch confBool
                                    case 0 % 2 unconfined diffusers in 1D
                                        msdFun=@(tau,p)cat(2,...
                                            m1(tau,p(1:2)),...
                                            m1(tau,p(3:4)));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1),p(5))-...
                                            c1(x,y(2),1-p(5));
                                        pID=[1,3:5,8];
                                    case 1 % 2 confined diffusers in 1D
                                        msdFun=@(tau,p)cat(2,...
                                            confMSD1(tau,p(1:3)),...
                                            confMSD1(tau,p([4,2,5])));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1),p(6))-...
                                            c1(x,y(2),1-p(6));
                                        pID=[1:5,8];
                                end
                            case 3
                                switch confBool
                                    case 0 % 3 unconfined diffusers in 1D
                                        msdFun=@(tau,p)cat(2,...
                                            m1(tau,p(1:2)),...
                                            m1(tau,p(3:4)),...
                                            m1(tau,p(5:6)));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1),p(7))-...
                                            c1(x,y(2),p(8))-...
                                            c1(x,y(3),1-p(7)-p(8));
                                        pID=[1,3:9];
                                        
                                    case 1 % 3 confined diffusers in 1D
                                        msdFun=@(tau,p)cat(2,...
                                            confMSD1(tau,p(1:3)),...
                                            confMSD1(tau,p([4,2,5]),...
                                            confMSD1(tau,p([6,2,7]))));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1),p(8))-...
                                            c1(x,y(2),p(9))-...
                                            c1(x,y(3),1-p(8)-p(9));
                                        pID=1:9;
                                end
                                
                        end
                    case 1 % 1 immobile term
                        switch nDiffs
                            case 1 % 1 diffuser
                                switch confBool
                                    case 0 % 1 unconfined diffuser with 1 immobile population in 1D
                                        msdFun=@(tau,p)...
                                            m1(tau,p(1:2));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y+p(4),p(3))-...
                                            c1(x,p(4),1-p(3));
                                        pID=[1,3,8,11];
                                    case 1 % 1 confined diffuser with 1 immobile population in 1D
                                        msdFun=@(tau,p)...
                                            confMSD1(tau,p(1:3));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y+p(5),p(4))-...
                                            c1(x,p(5),1-p(4));
                                        pID=[1:3,8,11];
                                end
                            case 2 % 2 diffusers
                                switch confBool
                                    case 0 % 2 unconfined diffusers with 1 immobile population in 1D
                                        msdFun=@(tau,p)cat(1,...
                                            m1(tau,p(1:2)),...
                                            m1(tau,p(3:4)));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1)+p(7),p(5))-...
                                            c1(x,y(2)+p(7),p(6))-...
                                            c1(x,p(7),1-p(5)-p(6));
                                        pID=[1,3:5,8,9,11];
                                    case 1 % 2 confined diffusers with 1 immobile population in 1D
                                        msdFun=@(tau,p)cat(1,...
                                            confMSD1(tau,p(1:3)),...
                                            confMSD1(tau,p([4,2,5])));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1)+p(8),p(6))-...
                                            c1(x,y(2)+p(8),p(7))-...
                                            c1(x,p(8),1-p(6)-p(7));
                                        pID=[1:5,8,9,11];
                                end
                            case 3
                                switch confBool
                                    case 0 % 3 unconfined diffusers with 1 immobile population in 1D
                                        msdFun=@(tau,p)cat(1,...
                                            m1(tau,p(1:2)),...
                                            m1(tau,p(3:4)),...
                                            m1(tau,p(5:6)));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1)+p(10),p(7))-...
                                            c1(x,y(2)+p(10),p(8))-...
                                            c1(x,y(3)+p(10),p(9))-...
                                            c1(x,p(10),1-p(7)-p(8)-p(9));
                                        pID=[1,3:11];
                                    case 1 % 3 confined diffusers with 1 immobile population in 1D
                                        msdFun=@(tau,p)cat(1,...
                                            confMSD1(tau,p(1:3)),...
                                            confMSD1(tau,p([4,2,5])),...
                                            confMSD1(tau,p([6,2,7])));
                                        cpdFun=@(x,y,p)...
                                            c1(x,y(1)+p(11),p(8))-...
                                            c1(x,y(2)+p(11),p(9))-...
                                            c1(x,y(3)+p(11),p(10))-...
                                            c1(x,p(11),1-p(8)-p(9)-p(10));
                                        pID=1:11;
                                end
                                
                        end
                end
        end
        
        pStart={pStart(pID)};
        bounds=[{LB(pID)},{UB(pID)}];
        dID = find(ismember(pID,[1,4,6]));
        aID = find(ismember(pID,8:10));
    case 0
        %% local fitting
        dID = 1;
        switch confBool
            case 0 % conf
                msdStart = pStart([1,3]);
                msdLB = [0, -inf];
                msdUB = [inf, inf];
                switch dim
                    case 1 % dim
                        msdFun = @(tau,p)m1(tau,p(1:2));
                        switch immBool
                            case 0 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1),1);
                                        cpdStart = pStart(1);
                                        cpdLB = 0;
                                        cpdUB = inf;
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1),p(3))-...
                                            c1(x,p(2),1-p(3));
                                        cpdStart = [pStart([1,4]), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1),p(4))-...
                                            c1(x,p(2),p(5))-...
                                            c1(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4,6]), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                end
                            case 1 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1)+p(2),p(3))-...
                                            c1(x,p(2),1-p(3));
                                        cpdStart = [pStart(1), pStart(8), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1)+p(3),p(4))-...
                                            c1(x,p(2)+p(3),p(5))-...
                                            c1(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4]), pStart(8), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1)+p(4),p(5))-...
                                            c1(x,p(2)+p(4),p(6))-...
                                            c1(x,p(3)+p(4),p(7))-...
                                            c1(x,p(4),1-p(5)-p(6)-p(7));
                                        cpdStart = [pStart([1,4,6]), pStart(8), .25, .25, .25];
                                        cpdLB = [0, 0, 0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, inf, 1, 1, 1];
                                end
                        end
                    case 2 % dim
                        msdFun = @(tau,p)m2(tau,p(1:2));
                        switch immBool
                            case 0 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1),1);
                                        cpdStart = pStart(1);
                                        cpdLB = 0;
                                        cpdUB = inf;
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1),p(3))-...
                                            c2(x,p(2),1-p(3));
                                        cpdStart = [pStart([1,4]), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1),p(4))-...
                                            c2(x,p(2),p(5))-...
                                            c2(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4,6]), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                end
                            case 1 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1)+p(2),p(3))-...
                                            c2(x,p(2),1-p(3));
                                        cpdStart = [pStart(1), pStart(8), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1)+p(3),p(4))-...
                                            c2(x,p(2)+p(3),p(5))-...
                                            c2(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4]), pStart(8), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1)+p(4),p(5))-...
                                            c2(x,p(2)+p(4),p(6))-...
                                            c2(x,p(3)+p(4),p(7))-...
                                            c2(x,p(4),1-p(5)-p(6)-p(7));
                                        cpdStart = [pStart([1,4,6]), pStart(8), .25, .25, .25];
                                        cpdLB = [0, 0, 0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, inf, 1, 1, 1];
                                end
                        end
                end
            case 1 % conf
                msdStart = pStart(1:3);
                msdLB = [0, 0, -inf];
                msdUB = [inf, inf, inf];
                switch dim
                    case 1 % dim
                        msdFun = @(tau,p)confMSD1(tau,p(1:3));
                        switch immBool
                            case 0 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1),1);
                                        cpdStart = pStart(1);
                                        cpdLB = 0;
                                        cpdUB = inf;
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1),p(3))-...
                                            c1(x,p(2),1-p(3));
                                        cpdStart = [pStart([1,4]), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1),p(4))-...
                                            c1(x,p(2),p(5))-...
                                            c1(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4,6]), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                end
                            case 1 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1)+p(2),p(3))-...
                                            c1(x,p(2),1-p(3));
                                        cpdStart = [pStart(1), pStart(8), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1)+p(3),p(4))-...
                                            c1(x,p(2)+p(3),p(5))-...
                                            c1(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4]), pStart(8), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c1(x,p(1)+p(4),p(5))-...
                                            c1(x,p(2)+p(4),p(6))-...
                                            c1(x,p(3)+p(4),p(7))-...
                                            c1(x,p(4),1-p(5)-p(6)-p(7));
                                        cpdStart = [pStart([1,4,6]), pStart(8), .25, .25, .25];
                                        cpdLB = [0, 0, 0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, inf, 1, 1, 1];
                                end
                        end
                    case 2 % dim
                        msdFun = @(tau,p)confMSD2(tau,p(1:3));
                        switch immBool
                            case 0 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1),1);
                                        cpdStart = pStart(1);
                                        cpdLB = 0;
                                        cpdUB = inf;
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1),p(3))-...
                                            c2(x,p(2),1-p(3));
                                        cpdStart = [pStart([1,4]), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1),p(4))-...
                                            c2(x,p(2),p(5))-...
                                            c2(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4,6]), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                end
                            case 1 % immPop
                                switch nDiffs
                                    case 1
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1)+p(2),p(3))-...
                                            c2(x,p(2),1-p(3));
                                        cpdStart = [pStart(1), pStart(8), .5];
                                        cpdLB = [0, 0, 0];
                                        cpdUB = [inf, inf, 1];
                                    case 2
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1)+p(3),p(4))-...
                                            c2(x,p(2)+p(3),p(5))-...
                                            c2(x,p(3),1-p(4)-p(5));
                                        cpdStart = [pStart([1,4]), pStart(8), .3, .3];
                                        cpdLB = [0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, 1, 1];
                                    case 3
                                        cpdFun=@(x,p)1-...
                                            c2(x,p(1)+p(4),p(5))-...
                                            c2(x,p(2)+p(4),p(6))-...
                                            c2(x,p(3)+p(4),p(7))-...
                                            c2(x,p(4),1-p(5)-p(6)-p(7));
                                        cpdStart = [pStart([1,4,6]), pStart(8), .25, .25, .25];
                                        cpdLB = [0, 0, 0, 0, 0, 0, 0];
                                        cpdUB = [inf, inf, inf, inf, 1, 1, 1];
                                end
                        end
                end
        end
        pStart=[{cpdStart},{msdStart}];
        bounds=[{cpdLB},{cpdUB},{msdLB},{msdUB}];
        aID = nDiffs+immBool+1:numel(cpdStart);
        aID(aID == 0) = [];
end

outStruct.cpdFun = cpdFun;
outStruct.msdFun = msdFun;
outStruct.pStart = pStart;
outStruct.bounds = bounds;
outStruct.dID = dID;
outStruct.aID = aID;
end

%% 2d square confinement model
function z=confMSD2(tau,p)
% global camerasd
d=p(1); l=p(2);

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));
for ii=1:2:2*400-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
z=l^2/3*(1-96/pi^4*temp)+p(3);
end

%% 1d square confinement model
function z=confMSD1(tau,p)
% global camerasd
d=p(1); l=p(2);

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));
for ii=1:2:2*400-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
z=l^2/6*(1-96/pi^4*temp)+p(3);
end