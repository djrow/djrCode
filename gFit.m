function [all_fitparam,all_fiterr,guesses,frameskip]=gFit(img,phmask,fProps)
% sym,bpparams,min_sep,frameskip,framenum)

if isempty(phmask)
    phmask=ones(size(img));
end

sym=fProps.symmetry;

img=double(img);
bpImg=guessPeaks(img,fProps);
bpImg=bpImg.*logical(phmask);

initwidth=2;        % Initial width (both x- and y-directions)

% 2*fWidth+1 is the width of the fitting window
fWidth=5;

if sym
    gFun=@(p,x)p(2)+p(1).*exp(-(x(:,1)-p(4)).^2/(2*p(3).^2)-...
        (x(:,2)-p(5)).^2./(2*p(3).^2));
elseif sym==2
    gFun=@(p,x)p(2)+p(1).*exp(-(x(:,1)-p(4)).^2/(2*p(3).^2)-...
        (x(:,2)-p(5)).^2/(2*p(6).^2));
end

[r,c,~]=find(bpImg);
numguess=size(r,1);

% Pre-alocate space for outputs
all_fitparam=zeros(numguess,7); % Addtional column stores intensity integral
all_fiterr=zeros(numguess,6);
xcentroid=zeros(1,numguess);
ycentroid=zeros(1,numguess);

%% Prepare for fitting (draw fitting box, determine centroid locations,
%  and linearize data)
for ii=1:numguess % Loop through each putative peak

    currx=c(ii);
    curry=r(ii);
    
    % Here the "bottom" is actually the top of the frame due to the
    % direction of the y-coordinate.
    smallbottom=curry-fWidth;
    smalltop=curry+fWidth;
    smallleft=currx-fWidth;
    smallright=currx+fWidth;
    
    % Make sure all boundaries of the small box lie inside the frame.
    if smallbottom<1; smallbottom=1; end
    if smalltop>size(img,1); smalltop=size(img,1); end
    if smallleft<1; smallleft=1; end
    if smallright>size(img,2); smallright=size(img,2); end
    
    tempImg=img(smallbottom:smalltop,smallleft:smallright);
    
    % mats for centroid calculation
    [x,y]=meshgrid(smallleft:smallright,smallbottom:smalltop);
    X=cat(2,x(:),y(:));
    
    % Calculate the centroid to use as inital center guess
    smallboxcounts=sum(tempImg(:));
    xcentroid(ii)=sum(sum(double(tempImg).*x))/smallboxcounts;
    ycentroid(ii)=sum(sum(double(tempImg).*y))/smallboxcounts;
    
    % Initialize fitting parameters:
    if sym==1
        initparam=double([max(tempImg(:))-min(tempImg(:)),...
            median(tempImg(:)),initwidth,xcentroid(ii),ycentroid(ii)]);
    elseif sym==2
        initparam=double([max(tempImg(:))-min(tempImg(:)),...
            median(tempImg(:)),initwidth,xcentroid(ii),ycentroid(ii),initwidth]);
        % [Amplitude, Offset, X-Width, Xcenter, Ycenter, Y-Width]
    end
    if initparam(2)>=initparam(1)
        initparam(2)=min(tempImg(:));
    end
    
    % Linearize for nlinfit
    [linear_r,linear_c]=ind2sub(size(tempImg),1:numel(tempImg));
    linearsmallboxdata=linear_c'+smallleft-1;
    linearsmallboxdata(:,2)=linear_r+smallbottom-1;
    
    %% PSF fitting with nlinfit
    [fitparam,residuals,jacobian]=nlinfit(linearsmallboxdata(:,1:2),...
        tempImg(:),gFun,initparam,options);
    
    % Estimating errors of fitting parameters
    fiterr=diff(nlparci(fitparam,residuals,'jacobian',jacobian),1,2);
    
    % Gather fit parameters for the input frame. Note that here to keep the
    % output for both symmetric and asymmetric gaussian fitting having the same
    % number of columns (for ease of post-processing), the 6th column (Y-Width)
    % of 'all_fitparam' for the symmetric gaussia fitting will have the same
    % value as the 3rd column (X-width). The same applies to 'all_fiterr'.
    if sym==2
        % fitting parameters and intensity integral. If you want to get the
        % intensity under the fitted Gaussian curve, use the "trapz" command.
        all_fitparam(ii,:)=[fitparam,sum(tempImg(:))];
        all_fiterr(ii,:)=fiterr;
    elseif sym==1
        all_fitparam(ii,:)=[fitparam,fitparam(:,3),sum(tempImg(:))];
        all_fiterr(ii,:)=[fiterr;fiterr(3,:)];
    end
end
guesses=cat(1,xcentroid,ycentroid);
end

function  bp_img=guessPeaks(img,fProps)
%-------------------------------------------------------------------------%
% INPUTS:
% img: One fluorescent image frame

% LP, HP, T, H_max: Band-pass filter parameters

% min_sep: If two putative peaks are separated by a distance smaller than
% this value (in pixels), only the brigher one will be used for the
% fitting. This value also represent the location intensities that will be
% used for the PSF fitting.

% show_bp_img: Set it to "1" to show the band-passed result with putative
% peaks highlighted (good for debugging), otherwise set it to "0". Default
% value = 0;

% OUTPUTS:
% bp_img: Binary image containing the guessed peaks
%-------------------------------------------------------------------------%


bpparams=fProps.o;
min_sep=fProps.l;
showRes=fProps.showRes;

LP=bpparams(1);
HP=bpparams(2);
T=bpparams(3);
hMax=bpparams(4);
lzero=bpparams(5);
isize=size(img);

bimg=bpass(img,LP,HP,T,lzero);

extimg=imextendedmax(bimg,hMax);

if all(extimg)
    bp_img=false(isize);
else
    bp_img=bwmorph(extimg,'shrink',inf);
end

if showRes==0
    subplot(2,2,1); imshow(img,[])
    title(['frame number ' num2str(framenum)])
    
    subplot(2,2,2); imshow(bimg,[])
    title('band passed img')
    
    subplot(2,2,3); imshow(extimg,[])
    title('extended maxima')
    
    subplot(2,2,4); imshow(bp_img,[])
    title('shrunk img')
    
    temp=input('next frame (enter), or frame number to skip to: ');
    if temp>framenum
        frameskip=temp;
    end
end

[guess_peak_r,guess_peak_c]=find(bp_img);

% Locations (rows and columns) of remaining pixels in the filtered image
% along with the intensities at corresponding locations in the original 
% frame "thisframe". 
[x,y]=meshgrid(-2:2,-2:2);
whichpix={[x(:),y(:)]};
xs=num2cell(guess_peak_r); 
ys=num2cell(guess_peak_c);
whichpix=cellfun(@(w,x,y) [w(:,1)+x,w(:,2)+y],...
    whichpix(ones(numel(guess_peak_r),1)),xs,ys,'uniformoutput',0);
for ii=1:numel(whichpix)
    whichpix{ii}(whichpix{ii}(:,1)>isize(1),1)=isize(1);
    whichpix{ii}(whichpix{ii}(:,2)>isize(2),2)=isize(2);
    whichpix{ii}(whichpix{ii}(:,1)<1,1)=1;
    whichpix{ii}(whichpix{ii}(:,2)<1,2)=1;
end

cisize={isize}; cimg={double(img)};
meanb=cellfun(@(x,y,z) mean(x(sub2ind(y,z(:,1),z(:,2)))),...
    cimg(ones(numel(guess_peak_c),1)),cisize(ones(numel(guess_peak_c),1)),...
    whichpix);

rci=cat(2,guess_peak_r,guess_peak_c,meanb);

% If peaks from guess_I are too close (separation distances smaller than
% half box size), only keep the one with the highest intensity.
if size(rci,1)>1
    pairs=nchoosek(1:size(rci,1),2);
    dists=sqrt((rci(pairs(:,1),1)-rci(pairs(:,2),1)).^2+...
        (rci(pairs(:,1),2)-rci(pairs(:,2),2)).^2);
    
    tooclose=pairs(dists<min_sep&dists>1e-8,:);
    utooclose=unique(tooclose);
    
    while numel(tooclose)>0
        goodinds=zeros(1,numel(utooclose));
        for ii=1:numel(utooclose)
            temp=logical(sum(tooclose==utooclose(ii),2));
            p=unique(tooclose(temp(:,[1,1])));
            [~,ind]=max(rci(p,3));
            goodinds(ii)=p(ind);
        end
        goodinds=[unique(goodinds),find(ismember(1:size(rci,1),utooclose)==0)];
        rci=rci(goodinds,:);
        
        if size(rci,1)>1
            pairs=nchoosek(1:size(rci,1),2);
            dists=sqrt((rci(pairs(:,1),1)-rci(pairs(:,2),1)).^2+...
                (rci(pairs(:,1),2)-rci(pairs(:,2),2)).^2);
            
            tooclose=pairs(dists<min_sep&dists>1e-8,:);
            utooclose=unique(tooclose);
        else
            tooclose=[];
        end
    end
end

bp_img=false(isize);
bp_img(sub2ind(isize,rci(:,1),rci(:,2)))=true;
end