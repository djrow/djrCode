function [all_fitparam,all_fiterr,guesses,frameskip]=...
    fit_asym_gauss_int2(f_img,phmask,sym,bpparams,min_sep,frameskip,framenum)
warning off all
%-------------------------------------------------------------------------%
% Fit intensities to an asymmetric gaussian function

% INPUTS:
% f_img: single frame of original fluorescent image.

% bp_img: band-passed image of 'f_img' returned by 'guess_peaks.m'.

% boxsize: Half box size for fitting pixel intensities to PSF. Default
% value is 7, meaning a 15-by-15 (px) box is used for each putative peak
% pixel.

% sym: Analytical form of PSF. Set to "1" for symmetric gaussian, and set
% to "2" for asymmetrical gaussian (default).

% OUTPUTS all_fitparam: Fitting parameters for all molecules from the
% current frame, as defined in the fitting function. The 7 fitting
% parameters are: Amplitude, Offset, X Width, X Center, Y Center, Y
% Width(will be equal to X Width for symmetric Gaussian fitting) and
% intensity sum. The number of rows in 'fitparam' is equal to the number of
% putative peaks (good and bad).

% all_fiterr: Statitstical errors (95% confidence) of fitting parameters

%-------------------------------------------------------------------------%
if nargin==5
    frameskip=inf;
    framenum=0;
end

f_img=double(f_img);
[bp_img,frameskip]=guess_peaks2(f_img,bpparams,min_sep,frameskip,framenum);
bp_img=bp_img.*logical(phmask);

initwidth=2;        % Initial width (both x- and y-directions)
options=statset('MaxFunEvals',5000,'MaxIter',5000);

gauss2dfunction_nlinfit_asym=@(p,x) p(2)+...
    p(1).*exp(-(x(:,1)-p(4)).^2/(2*p(3).^2)-(x(:,2)-p(5)).^2/(2*p(6).^2));
gauss2dfunction_nlinfit_sym=@(p,x) p(2)+...
    p(1).*exp(-(x(:,1)-p(4)).^2/(2*p(3).^2)-(x(:,2)-p(5)).^2./(2*p(3).^2));
testfun=@(p,x,y)p(2)+p(1).*...
    exp(-(x-p(4)).^2/(2*p(3).^2)-(y-p(5)).^2/(2*p(6).^2));

[r,c,~]=find(bp_img);
numguess=size(r,1);

% Pre-alocate space for outputs
all_fitparam=zeros(numguess,7); % Addtional column stores intensity integral
all_fiterr=zeros(numguess,6);
xcentroid=zeros(1,numguess);
ycentroid=zeros(1,numguess);

%% Prepare for fitting (draw fitting box, determine centroid locations,
%  and linearize data)
for ii=1:numguess % Loop through each putative peak
    if frameskip>framenum
        skipshow=1;
    else
        skipshow=0;
    end
    currx=c(ii);
    curry=r(ii);
    
    % Here the "bottom" is actually the top of the frame due to the
    % direction of the y-coordinate.
    smallbottom=curry-min_sep;
    smalltop=curry+min_sep;
    smallleft=currx-min_sep;
    smallright=currx+min_sep;
    
    % Make sure all boundaries of the small box lie inside the frame.
    if smallbottom<1; smallbottom=1; end
    if smalltop>size(f_img,1); smalltop=size(f_img,1); end
    if smallleft<1; smallleft=1; end
    if smallright>size(f_img,2); smallright=size(f_img,2); end
    
    smallboxdata=f_img(smallbottom:smalltop,smallleft:smallright);
    
    % mats for centroid calculation
    [x,y]=meshgrid(smallleft:smallright,smallbottom:smalltop);
    
    % Calculate the centroid to use as inital center guess
    smallboxcounts=sum(smallboxdata(:));
    xcentroid(ii)=sum(sum(double(smallboxdata).*x))/smallboxcounts;
    ycentroid(ii)=sum(sum(double(smallboxdata).*y))/smallboxcounts;
    
    % Initialize fitting parameters:
    if sym==1
        initparam=double([max(smallboxdata(:))-min(smallboxdata(:)),...
            median(smallboxdata(:)),initwidth,xcentroid(ii),ycentroid(ii)]);
    elseif sym==2
        initparam=double([max(smallboxdata(:))-min(smallboxdata(:)),...
            median(smallboxdata(:)),initwidth,xcentroid(ii),ycentroid(ii),initwidth]);
        % [Amplitude, Offset, X-Width, Xcenter, Ycenter, Y-Width]
    end
    if initparam(2)>=initparam(1)
        initparam(2)=min(smallboxdata(:));
    end
    
    % Linearize for nlinfit
    [linear_r,linear_c]=ind2sub(size(smallboxdata),1:numel(smallboxdata));
    % Note that I intentionally switched the column and row indices here such
    % that the location subscripts increase column-wise first and then
    % row-wise. This goes along with the subsequent line where 'smallboxdata'
    % is transposed first and then reshaped(linearized).
    linearsmallboxdata=linear_c'+smallleft-1;
    linearsmallboxdata(:,2)=linear_r+smallbottom-1;
    
    %% PSF fitting with nlinfit
    if sym==2 % Asymmetric Gaussian fitting
        [fitparam,residuals,jacobian]=nlinfit(linearsmallboxdata(:,1:2),...
            smallboxdata(:),gauss2dfunction_nlinfit_asym,initparam,options);
    elseif sym==1 % Symmetric Gaussian fitting
        [fitparam,residuals,jacobian]=nlinfit([linearsmallboxdata(:,1),...
            linearsmallboxdata(:,2)],linearsmallboxdata(:,3),...
            gauss2dfunction_nlinfit_sym,initparam,options);
    end
    
    %       % check fits?
    if skipshow==0
        display(['The second image should be just the background noise in '...
            'the first one.'])
        subplot(1,1,1)
        subplot(2,1,1); a=pcolor(smallboxdata);
        set(a,'edgecolor','none')
        title('data')
        cm=get(gca,'clim'); axis image
        subplot(2,1,2); a=pcolor(smallboxdata-testfun(fitparam,x,y));
        set(a,'edgecolor','none')
        cm=cm-min(cm)+min(min(smallboxdata-testfun(fitparam,x,y)));
        set(gca,'clim',cm); axis image
        display('These two columns should be similar: ')
        display(permute(cat(1,fitparam,initparam),[3,1,2]));
        title('residuals')
        temp=input('next frame (enter), or frame number to skip to: ');
        if temp>framenum
            frameskip=temp;
        end
    end
    
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
        all_fitparam(ii,:)=[fitparam,sum(smallboxdata(:))];
        all_fiterr(ii,:)=fiterr;
    elseif sym==1
        all_fitparam(ii,:)=[fitparam,fitparam(:,3),sum(smallboxdata(:))];
        all_fiterr(ii,:)=[fiterr,fiterr(:,3)];
    end
end
guesses=cat(1,xcentroid,ycentroid);
end

function  [bp_img,frameskip]=guess_peaks2(img,bpparams,min_sep,frameskip,framenum)
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
if nargin==3
    skipshow=1;
    frameskip=0;
    framenum=0;
elseif frameskip>framenum
    skipshow=1;
else
    skipshow=0;
end

LP=bpparams(1);
HP=bpparams(2);
T=bpparams(3);
H_max=bpparams(4);
lzero=bpparams(5);

isize=size(img);

bimg=bpass(img,LP,HP,T,lzero);
extimg=imextendedmax(bimg,H_max);
if all(extimg)
    bp_img=false(isize);
else
    bp_img=bwmorph(extimg,'shrink',inf);
end

if skipshow==0% ||numel(find(bp_img))>20
    display(['edit the bpparams vector in the masterfit file to change '...
        'these images.'])
    display('enlist the help of yi or david the first time.')
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
xs=num2cell(guess_peak_r); ys=num2cell(guess_peak_c);
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

% % Get rid of the artificial pixel (if picked up by band-passing) at the
% % center of the frame.
% if size(rci,1)==1
%     rci=rci(rci(:,1)~=round(isize(1)/2)+1&...
%         rci(:,2)~=round(isize(2)/2),:);
% end

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