function [p,p95,nrmse]=gaussfit(data,usedparameters,plotres,ang)
% hw=10;
% data=data(floor(size(data,1)/2-hw)+1:floor(size(data,1)/2+hw)+1,...
%     floor(size(data,2)/2-hw)+1:floor(size(data,2)/2+hw)+1);
% perc=.25;
% drange=nanmax(data(:))-nanmin(data(:));
% data(data<drange*perc+min(data(:)))=nan;

% tmax=max(data(:))-min(data(:));
% find(data==data<tmax/2+tmax/1e2&data>tmax/2-tmax/1e2)

% pstart=[(nanmax(data(:))-nanmin(data(:)))/2,nanmin(data(:)),size(data)/2,20,20,0,nanmax(data(:))/2,10];
pstart=[nanmax(data(:))-1,1,0,0,10,20,0,nanmax(data(:))/2,10];

[x,y]=ndgrid(1:size(data,1),1:size(data,2));
x=x-1-floor(size(data,1)/2);
y=y-1-floor(size(data,2)/2);

xdata=cat(2,x(:),y(:));

options=statset('nlinfit');

if usedparameters==1234567;
    % all parameters
    fun=@(p,x)p(1)*exp(-((x(:,1)-p(3))*cos(p(7))+(x(:,2)-p(4))*sin(p(7))).^2/p(5)^2).*...
        exp(-(-(x(:,1)-p(3))*sin(p(7))+(x(:,2)-p(4))*cos(p(7))).^2/p(6)^2)+p(2);
    pstart=pstart([1,2,3,4,5,6,7]);
    
elseif usedparameters==123456;
    % no angle
    if ~exist('ang','var')
        ang=0;
    end
    fun=@(p,x)p(1)*exp(-((x(:,1)-p(3))*cos(ang)+(x(:,2)-p(4))*sin(ang)).^2/p(5)^2).*...
        exp(-(-(x(:,1)-p(3))*sin(ang)+(x(:,2)-p(4))*cos(ang)).^2/p(6)^2)+p(2);
    pstart=pstart([1,2,3,4,5,6]);
    
elseif usedparameters==12345;
    % symmetric
    fun=@(p,x)p(1)*exp(-((x(:,1)-p(3))).^2/p(5)^2).*...
        exp(-((x(:,2)-p(4))).^2/p(5)^2)+p(2);
    pstart=pstart([1,2,3,4,5]);
    
elseif usedparameters==125;
    % centered, symmetric
    fun=@(p,x)p(1)*exp(-x(:,1).^2/2/p(3)^2).*...
        exp(-x(:,2).^2/2/p(3)^2)+p(2);
    pstart=pstart([1,2,5]);
    
elseif usedparameters==134567
    % offset of 0
    fun=@(p,x)p(1)*exp(-((x(:,1)-p(2))*cos(p(6))+(x(:,2)-p(3))*sin(p(6))).^2/p(4)^2).*...
        exp(-(-(x(:,1)-p(2))*sin(p(6))+(x(:,2)-p(3))*cos(p(6))).^2/p(5)^2);
    pstart=pstart([1,3,4,5,6,7]);
    
elseif usedparameters==15
    % offset of 0, centered, symmetric
    fun=@(p,x)p(1)*exp(-x(:,1).^2/p(2)^2).*...
        exp(-x(:,2).^2/p(2)^2);
    pstart=pstart([1,5]);
    
elseif usedparameters==156
    % offset of 0, centered
    if ~exist('ang','var')
        ang=0;
    end
    
    fun=@(p,x)p(1)*exp(-(x(:,1)*cos(ang)+x(:,2)*sin(ang)).^2/p(2)^2).*...
        exp(-(-x(:,1)*sin(ang)+x(:,2)*cos(ang)).^2/p(3)^2);
    pstart=pstart([1,5,6]);
    
elseif usedparameters==1256
    % offset of 0, centered
    if ~exist('ang','var')
        ang=0;
    end
    
    fun=@(p,x)p(1)*exp(-(x(:,1)*cos(ang)+x(:,2)*sin(ang)).^2/p(3)^2).*...
        exp(-(-x(:,1)*sin(ang)+x(:,2)*cos(ang)).^2/p(4)^2)+p(2);
    pstart=pstart([1,2,5,6]);
    
elseif usedparameters==125689
    % offset of 0, centered
    if ~exist('ang','var')
        ang=0;
    end
    
    fun=@(p,x)abs(p(1))*exp(-(x(:,1)*cos(ang)+x(:,2)*sin(ang)).^2/p(3)^2).*...
        exp(-(-x(:,1)*sin(ang)+x(:,2)*cos(ang)).^2/p(4)^2)+p(2)+...
        abs(p(5))*exp(-(x(:,1).^2+x(:,2).^2)/p(6)^2);
    pstart=pstart([1,2,5,6,8,9]);
    
end
nonnanvals=~isnan(data); nonnanvals=nonnanvals(:);
% [p,~,res,~,~,~,jacob]=lsqcurvefit(fun,pstart,reshape(xdata(nonnanvals(:,[1,1])),[sum(nonnanvals),2]),...
%     data(nonnanvals),lb,ub,options);
mdl=fitnlm(reshape(xdata(nonnanvals(:,[1,1])),sum(nonnanvals),[]),...
    data(nonnanvals),fun,pstart,'Options',options);

p=mdl.Coefficients{:,1};
p95=diff(coefCI(mdl),1,2);
nrmse=mdl.RMSE/numel(data)/p(1);

if plotres==1
    funval=reshape(fun(p,xdata),size(x));
    startfunval=reshape(fun(pstart,xdata),size(x));
    
    subplot(131); pcolor(data);
    title('data')
    axis image; shading flat
    cl=get(gca,'clim');
    
    subplot(132); pcolor(startfunval);
    title('startfit')
    axis image; shading flat
    
    subplot(133); pcolor(data-funval);
    title('residuals')
    axis image; shading flat
    cl=cl-min(cl)+min(data(:)-funval(:));
    set(gca,'clim',cl);
    
    keyboard
end
end