function timecorr=stics3(imstack,mask,maxtau,whichfun)
% calculate the full space-time correlation of an image stack.
% size(imstack)=[x, y, frame number]
%
% mask is a binary image of pixels to use in each frame of the movie
%
% dim is 1 or 2, if the movie is of one dimensional or two dimensional
% diffusion

%%

imstack=double(imstack);
vidsize=[size(imstack),1];

if ~exist('whichfun','var')
    whichfun=[];
end

if ~exist('mask','var')
    mask=ones(vidsize(1:2));
end
mask=logical(mask);

vidmean=mean(imstack(:));

fun{1}=@(x,y)pixremfun(x,y);    % remove pixels outside mask
fun{2}=@(x,y)fremfun(x,y);      % remove stationary objects via fft

% apply filters in arbitrary order
for ii=whichfun(whichfun~=4)
    imstack=fun{ii}(imstack,mask);
end

imstack=imstack-mean(imstack(:));

% do all the time correlations
if any(whichfun==4)     % nonoverlapping
    vfft=fft2(imstack,2*vidsize(1)+1,2*vidsize(2)+1);
    
    timecorr=zeros(size(vfft,1),size(vfft,2),maxtau);
    for kk=1:maxtau
        ind1=1:kk:vidsize(3);
        ind2=ind1(2:end);
        ind1=ind1(1:end-1);
        timecorr(:,:,kk)=mean(real(fftshift(fftshift(...
            ifft2(vfft(:,:,ind2).*conj(vfft(:,:,ind1))),1),2)),3);
    end
    timecorr=timecorr/numel(vfft(:,:,1));
    
else                    % overlapping
    famps=abs(fft(fft(fft(imstack,2*vidsize(1)+1,1),2*vidsize(2)+1,2),[],3)).^2;
    timecorr=fftshift(fftshift(real(ifft(ifft(ifft(famps...
        ,[],1),[],2),[],3)),1),2)/numel(famps);
    timecorr(:,:,1)=[];
end

end

function imstack=pixremfun(imstack,mask)
% display('removing noise pixels')

nframes=size(imstack,3);
mmean=imstack(mask(:,:,ones(1,nframes)));

imstack(~mask(:,:,ones(1,nframes)))=mean(mmean(:));
end

function imstack=fremfun(imstack,mask) %#ok<INUSD>
% display('removing immobile artifacts')
% plot(squeeze(mean(mean(imstack,1),2))); hold all

% removing the time-invariant term sets the minimum frequency in the movie
% to fmin, in frames per second.
% fmin=1*2/movielengthinseconds;    % 1/15 1/s if the movie is 30 seconds long

% which means that things that stay stationary for longer than 1/fmin will
% be removed.
% tstatmax=1/fmin;                  % 15 seconds

immean=mean(imstack(:));
% imsize=size(imstack);

imstack=fft(double(imstack-immean),2*size(imstack,3)+1,3);

% plot(linspace(0,2/.01,imsize(3)),squeeze(mean(mean(imstack(:,:,1:imsize(3))))))

imstack(:,:,1)=0;

imstack=real(ifft(imstack,[],3))+immean;
imstack=imstack(:,:,1:(size(imstack,3)-1)/2);
% plot(squeeze(mean(mean(imstack,1),2))-1)
end