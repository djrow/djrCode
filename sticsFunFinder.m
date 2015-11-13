function [sticsFun,msdFun,pId]=sticsFunFinder(nDiffs,immPop,dType,th)
g1=@(x,y,p)p(1)*exp(-((x(:,:,1)-p(2)).^2+...
    (x(:,:,2)-p(3)).^2)./(2*y))./(2*pi*y);

switch immPop
    case 0 % no immobile term
        switch nDiffs
            case 1 % 1 diffuser
                switch dType
                    case 'unconfined' % 1 unconfined diffuser
                        msdFun=@(tau,p)...
                            2*p(1)*tau+p(2);
                        sticsFun=@(x,y,p) p(3)+...
                            g1(x,y(1),p(4:6));
                        
                        pId=[1,3,9:12];
                    case 'confined' % 1 confined diffuser
                        msdFun=@(tau,p)cat(1,...
                            longmsd1d(p(1:3),tau),...
                            longmsd1d(p([1,4,5]),tau));
                        sticsFun=@(x,y,p) p(6)+...
                            g2(x,y(1:2),p(7:9),th);
                        
                        pId=[1:5,9:12];
                end
            case 2 % 2 diffusers
                switch dType
                    case 'unconfined' % 2 unconfined diffusers
                        msdFun=@(tau,p) cat(1,...
                            2*p(1)*tau+p(2),...
                            2*p(3)*tau+p(4));
                        sticsFun=@(x,y,p) p(5)+...
                            g1(x,y(1),p(6:8))+...
                            g1(x,y(2),p([9,7,8]));
                        
                        pId=[1,3,6,8,9:13];
                    case 'confined' % 2 confined diffusers
                        msdFun=@(tau,p) cat(1,...
                            longmsd1d(p(1:3),tau),...
                            longmsd1d(p([1,4,5]),tau),...
                            longmsd1d(p([6,2,7]),tau),...
                            longmsd1d(p([6,4,8]),tau));
                        sticsFun=@(x,y,p) p(9)+...
                            g2(x,y(1:2),p(10:12),th)+...
                            g2(x,y(3:4),p([13,11,12]),th);
                        
                        pId=[1:8,9:12,13];
                end
        end
    case 1 % 1 immobile term
        switch nDiffs
            case 1 % 1 diffuser
                switch dType
                    case 'unconfined' % 1 unconfined diffuser with 1 immobile population
                        msdFun=@(tau,p)...
                            2*p(1)*tau+p(2);
                        sticsFun=@(x,y,p) p(3)+...
                            g1(x,y,p(4:6))+...
                            g1(x,p(8),p([7,5,6]));
                        
                        pId=[1,3,9:12,13,15];
                    case 'confined' % 1 confined diffuser with 1 immobile population
                        msdFun=@(tau,p)cat(1,...
                            longmsd1d(p(1:3),tau),...
                            longmsd1d(p([1,4,5]),tau));
                        sticsFun=@(x,y,p) p(6)+...
                            g2(x,y,p(7:9),th)+...
                            g1(x,p(11),p([10,8,9]));
                        
                        pId=[1:5,9:12,14,15];
                end
            case 2 % 2 diffusers
                switch dType
                    case 'unconfined' % 2 unconfined diffusers with 1 immobile population
                        msdFun=@(tau,p) cat(1,...
                            2*p(1)*tau+p(2),...
                            2*p(3)*tau+p(4));
                        sticsFun=@(x,y,p) p(5)+...
                            g1(x,y,p(6:8))+...
                            g1(x,y,p([9,7,8]))+...
                            g1(x,p(11),p([10,7,8]));
                        
                        pId=[1,3,6,8,9:15];
                    case 'confined' % 2 confined diffusers with 1 immobile population
                        msdFun=@(tau,p) cat(1,...
                            longmsd1d(p(1:3),tau),...
                            longmsd1d(p([1,4,5]),tau),...
                            longmsd1d(p([6,2,7]),tau),...
                            longmsd1d(p([6,4,8]),tau));
                        sticsFun=@(x,y,p) p(9)+...
                            g2(x,y(1:2),p(10:12),th)+...
                            g2(x,y(3:4),p([13,11,12]),th)+...
                            g1(x,p(15),p([14,11,12]));
                        
                        pId=1:15;
                end
        end
end
end

%% rotated asymmetric gaussian model
function z=g2(x,y,p,th)
xp(:,:,1)=(x(:,:,1)-mean(mean(x(:,:,1))))*cos(th)+(x(:,:,2)-mean(mean(x(:,:,2))))*sin(th);
xp(:,:,2)=-(x(:,:,1)-mean(mean(x(:,:,1))))*sin(th)+(x(:,:,2)-mean(mean(x(:,:,2))))*cos(th);

z=p(1)*exp(-(xp(:,:,1)-p(2)).^2/2/y(1)-(xp(:,:,2)-p(3)).^2./2/y(2))./(2*pi*sqrt(y(1)*y(2)));
end
%% square confinement model
function z=longmsd1d(p,x)
% global camerasd
d=abs(p(1)); l=p(2); tau=x;

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