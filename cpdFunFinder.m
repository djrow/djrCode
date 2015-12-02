function [cpdFun,msdFun,pId]=cpdFunFinder(dim,nDiffs,immPop,dType)
% this code creates the correct anonymous functions to be used by cpdGlobal
% given the dimensionality (1 or 2), the number of diffusive populations (1
% or 2), the presence or absence of an immobile population (0 or 1), and
% the diffusion type ('unconfined' or 'confined').

% partial cpd function
c1=@(x,y,p)p*exp(-x./y);

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% 1d unconfined msd function
m1=@(t,p)2*p(1)*t+p(2);

switch dim
    case 2
        %% 2D diffusion
        switch immPop
            case 0 % no immobile term
                switch nDiffs
                    case 1 % 1 diffuser
                        switch dType
                            case 'unconfined' % 1 unconfined diffuser in 2D
                                msdFun=@(tau,p)...
                                    m2(tau,p(1:2));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y,1);
                                pId=[1,3];
                            case 'confined' % 1 confined diffuser in 2D
                                msdFun=@(tau,p)...
                                    confMSD2(p(1:3),tau);
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y,1);
                                pId=1:3;
                        end
                    case 2 % 2 diffusers
                        switch dType
                            case 'unconfined' % 2 unconfined diffusers in 2D
                                msdFun=@(tau,p)cat(1,...
                                    m2(tau,p(1:2)),...
                                    m2(tau,p(3:4)));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1),p(5))-...
                                    c1(x,y(2),1-p(5));
                                pId=[1,3:6];
                            case 'confined' % 2 confined diffusers in 2D
                                msdFun=@(tau,p)cat(1,...
                                    confMSD2(tau,p(1:3)),...
                                    confMSD2(tau,p([4,2,5])));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1),p(6))-...
                                    c1(x,y(2),1-p(6));
                                pId=1:6;
                        end
                end
            case 1 % 1 immobile term
                switch nDiffs
                    case 1 % 1 diffuser
                        switch dType
                            case 'unconfined' % 1 unconfined diffuser with 1 immobile population in 2D
                                msdFun=@(tau,p)...
                                    m2(tau,p(1:2));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y+p(4),p(3))-...
                                    c1(x,p(4),1-p(3));
                                pId=[1,3,6,8];
                            case 'confined' % 1 confined diffuser with 1 immobile population in 2D
                                msdFun=@(tau,p)cat(1,...
                                    confMSD2(tau,p(1:3)),...
                                    confMSD2(tau,p([4,2,5])));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y+p(7),p(6))-...
                                    c1(x,p(7),1-p(6));
                                pId=[1:6,8];
                        end
                    case 2 % 2 diffusers
                        switch dType
                            case 'unconfined' % 2 unconfined diffusers with 1 immobile population in 2D
                                msdFun=@(tau,p)cat(1,...
                                    m2(tau,p(1:2)),...
                                    m2(tau,p(3:4)));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1)+p(7),p(5))-...
                                    c1(x,y(2)+p(7),p(6))-...
                                    c1(x,p(7),1-p(5)-p(6));
                                pId=[1,3:8];
                            case 'confined' % 2 confined diffusers with 1 immobile population in 2D
                                msdFun=@(tau,p)cat(1,...
                                    confMSD2(tau,p(1:3)),...
                                    confMSD2(tau,p([4,2,5])));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1)+p(8),p(6))-...
                                    c1(x,y(2)+p(8),p(7))-...
                                    c1(x,p(8),1-p(6)-p(7));
                                pId=1:8;
                        end
                end
        end
    case 1
        %% 1D diffusion
        switch immPop
            case 0 % no immobile term
                switch nDiffs
                    case 1 % 1 diffuser
                        switch dType
                            case 'unconfined' % 1 unconfined diffuser in 1D
                                msdFun=@(tau,p)...
                                    m1(tau,p(1:2));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y,1);
                                pId=[1,3];
                            case 'confined' % 1 confined diffuser in 1D
                                msdFun=@(tau,p)...
                                    confMSD1(p(1:3),tau);
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y,1);
                                pId=1:3;
                        end
                    case 2 % 2 diffusers
                        switch dType
                            case 'unconfined' % 2 unconfined diffusers in 1D
                                msdFun=@(tau,p)cat(1,...
                                    m2(tau,p(1:2)),...
                                    m2(tau,p(3:4)));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1),p(5))-...
                                    c1(x,y(2),1-p(5));
                                pId=[1,3:6];
                            case 'confined' % 2 confined diffusers in 1D
                                msdFun=@(tau,p)cat(1,...
                                    confMSD2(tau,p(1:3)),...
                                    confMSD2(tau,p([4,2,5])));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1),p(6))-...
                                    c1(x,y(2),1-p(6));
                                pId=1:6;
                        end
                end
            case 1 % 1 immobile term
                switch nDiffs
                    case 1 % 1 diffuser
                        switch dType
                            case 'unconfined' % 1 unconfined diffuser with 1 immobile population in 1D
                                msdFun=@(tau,p)...
                                    m1(tau,p(1:2));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y+p(4),p(3))-...
                                    c1(x,p(4),1-p(3));
                                pId=[1,3,6,8];
                            case 'confined' % 1 confined diffuser with 1 immobile population in 1D
                                msdFun=@(tau,p)cat(1,...
                                    confMSD1(tau,p(1:3)),...
                                    confMSD1(tau,p([4,2,5])));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y+p(7),p(6))-...
                                    c1(x,p(7),1-p(6));
                                pId=[1:6,8];
                        end
                    case 2 % 2 diffusers
                        switch dType
                            case 'unconfined' % 2 unconfined diffusers with 1 immobile population in 1D
                                msdFun=@(tau,p)cat(1,...
                                    m1(tau,p(1:2)),...
                                    m1(tau,p(3:4)));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1)+p(7),p(5))-...
                                    c1(x,y(2)+p(7),p(6))-...
                                    c1(x,p(7),1-p(5)-p(6));
                                pId=[1,3,4:8];
                            case 'confined' % 2 confined diffusers with 1 immobile population in 1D
                                msdFun=@(tau,p)cat(1,...
                                    confMSD1(tau,p(1:3)),...
                                    confMSD1(tau,p([4,5,2])));
                                cpdFun=@(x,y,p)1-...
                                    c1(x,y(1)+p(8),p(6))-...
                                    c1(x,y(2)+p(8),p(7))-...
                                    c1(x,p(8),1-p(6)-p(7));
                                pId=1:8;
                        end
                end
        end
end
end

%% square confinement model
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

%% square confinement model
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