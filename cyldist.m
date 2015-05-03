function [m,p,raddist,leng]=cyldist(varargin)
% The best line, will pass through m and have direction cosines in p, so
% it can be expressed parametrically as x=m(1)+p(1)*t, y=m(2)+p(2)*t,
% and z=m(3)+p(3)*t, where t is the distance along the line from the mean
% point at m. also, leng is the length of the cell determined by the
% difference between the largest positive and negative projections along
% the center line.
%
% size(varargin{1})=[npoints,3];
% or
% size(varargin{1:3})=[npoints,1];

umperpixel=.049;

if nargin==1
    x=varargin{1}(:,1);
    y=varargin{1}(:,2);
    try
        z=varargin{1}(:,3);
    catch
        z=zeros(size(x));
    end
elseif nargin==3
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
end

[nx,mx]=size(x);
[ny,my]=size(y);
[nz,mz]=size(z);

if (mx~=1)||(my~=1)||(mz~=1)||(ny~=nx)||(nz~=nx)
    error('The arguments must be the same length.')
end

ponts=[x,y,z]*umperpixel;                   % scale position vectors

m=mean(ponts,1);                            % center of mass
m(3)=(max(ponts(:,3))+min(ponts(:,3)))/2;   % fix z center

% center all the position vectors
w=[ponts(:,1)-m(1),ponts(:,2)-m(2),ponts(:,3)-m(3)];

a=(1/nx)*(w'*w);                         % sum of projections for x, y, and z
[u,~,~]=svd(a);                         % eigenvector for this matrix
p=u(:,1)';                              % get eigenvector for largest eigenvalue,
% which is the direction of largest
% displacement

% fix the center line to the x-y plane
p(3)=0;
p=p/norm(p);

% project (all) positions onto the center line
proj=sum(bsxfun(@plus,ponts,-m).* ...
    bsxfun(@times,ones(size(ponts,1),3),p),2);

% pythagoras tells us the length of the third side of a right triangle
raddist=sqrt(sum((bsxfun(@plus,ponts,-m)).^2,2)- ...
    abs(proj).^2);

% the center line
t=linspace(min(proj),max(proj));
leng=t(end)-t(1);
x1=m(1)+p(1)*t;
y1=m(2)+p(2)*t;
z1=m(3)+p(3)*t;



% make a figure
% subplot(1,2,1)
% scatter3(ponts(:,1),ponts(:,2),ponts(:,3),'.')
% hold on
% line(x1,y1,z1,'color','r');
% hold off
% axis image

% subplot(1,2,2)
% [binnedvals,bincenters]=hist(raddist,nbins);
% normfactor=2*pi*(max(proj)-min(proj))*(bincenters.^2-[0,bincenters(1:end-1)].^2);
% 
% plot(bincenters,binnedvals./normfactor);
% axis tight
% xlabel('radial distance (um)')
% ylabel('particle density (1/um^3)')
end