function timecorr=stics3(imstack,mask,maxtau,whichfun)
% calculate the full space-time correlation of an image stack.
% size(imstack)=[x, y, frame number]
%
% mask is a binary image of pixels to use in each frame of the movie
%
% dim is 1 or 2, if the movie is of one dimensional or two dimensional
% diffusion

%%

vidsize=[size(imstack),1];

if ~exist('whichfun','var')
    whichfun=[];
end

if ~exist('maxtau','var')
    maxtau=20;
end

if isempty(mask)
    mask=true(vidsize(1:2));
end

fun{1}=@(x,y)pixremfun(x,y);    % remove pixels outside mask
fun{2}=@(x,y)fremfun(x,y);      % remove stationary objects via fft
fun{3}=@(x,y)padmean(x,y);      % pad every image with its mean

% apply filters in arbitrary order
for ii=whichfun(whichfun~=4)
    imstack=fun{ii}(imstack,mask);
end

% do all the time correlations
if any(whichfun==4)     % nonoverlapping
    vfft=fft2(imstack);
    
    timecorr=zeros(size(vfft,1),size(vfft,2),maxtau+1);
    
    timecorr(:,:,1)=mean(bsxfun(@times,real(fftshift(fftshift(...
        ifft2(vfft.*conj(vfft)),1),2)),...
        1./mean(mean(imstack))./mean(mean(imstack))),3);
    for kk=1:maxtau
        ind1=1:kk:vidsize(3);
        ind2=ind1(2:end);
        ind1=ind1(1:end-1);
        timecorr(:,:,kk+1)=mean(bsxfun(@times,real(fftshift(fftshift(...
            ifft2(vfft(:,:,ind2).*conj(vfft(:,:,ind1))),1),2)),...
            1./mean(mean(imstack(:,:,ind1)))./mean(mean(imstack(:,:,ind2)))),3);
    end
    timecorr=timecorr/numel(vfft(:,:,1))-1;
    
else                    % overlapping
        
    famps=abs(fft(fft(fft(imstack,[],1),[],2),[],3)).^2;
    timecorr=fftshift(fftshift(real(ifft(ifft(ifft(famps...
        ,[],1),[],2),[],3)),1),2)/numel(famps)/mean(imstack(:))^2-1;
    timecorr=timecorr(:,:,1:maxtau+1);
end

end

function imstack=pixremfun(imstack,mask)
imsize=size(imstack);

% replace pixels outside the mask with the average value inside the mask in
% each frame
mmean=mean(reshape(imstack(mask(:,:,ones(1,imsize(3)))),[],imsize(3)));
imstack(~mask(:,:,ones(1,imsize(3))))=mmean(ones(1,sum(~mask(:))),:);
end

function imstack=fremfun(imstack,mask) %#ok<INUSD>
% removing the time-invariant term sets the minimum frequency in the movie
% to fmin, in frames per second.
% fmin=1*2/movielengthinseconds;    % 1/15 1/s if the movie is 30 seconds long

imsize=size(imstack);

imstack=fft(imstack-mean(imstack(:)),2*imsize(3)+1,3);

% plot(linspace(0,2/.01,imsize(3)),squeeze(mean(mean(abs(imstack(:,:,1:imsize(3)))))))

imstack(:,:,1)=0;

imstack=real(ifft(imstack,[],3));
imstack=imstack(1:imsize(1),1:imsize(2),1:imsize(3));
% plot(squeeze(mean(mean(imstack,1),2))-1)
end

function imstack=padmean(imstack,mask)
imstack=mat2cell(imstack,size(imstack,1),size(imstack,2),ones(1,size(imstack,3)));
imm=cellfun(@(x)mean(mean(x(mask))),imstack,'uniformoutput',false);

imstack=cellfun(@(x,y)padarray(x,floor(size(x)/2),y),imstack,imm,...
    'uniformoutput',false);
imstack=cat(3,imstack{:});
end