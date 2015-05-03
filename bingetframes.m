function [frames,nframes,wsize,threshval]=bingetframes(vidloc,framenumbers,rect)
% pull out frames from a binary movie file created with loadnd2.m
% framenumbers can be any vector of integers, and rect is [xmin,ymin,xmax,ymax].
% either or both can be omitted to include all frames or all of every frame,...
% respectively. omit by use of '[]'

fid=fopen(vidloc,'r');                 % open the file for reading
if fid==-1;
    display('no such file exists')
    frames=[];
    nframes=[];
    wsize=[];
    return
end

fseek(fid,-32,'eof');                   % move to the location of the info
threshval=fread(fid,1,'double');        % read the graythresh value
wsize=fread(fid,2,'double');            % read frame size
nframes=fread(fid,1,'double');          % read number of frames

if exist('rect','var')
    if isempty(rect)                        % if no ROI size specified
        rect=[1,1,wsize(1),wsize(2)];       % read all pixels.
    end
else
    rect=[1,1,wsize(1),wsize(2)];
end

% ram space required kept constant
bucketsize=ceil(1./prod([rect(3)-rect(1),rect(4)-rect(2)])/3e-7);

% if numel(framenumbers)==1               % if only beginning frame specified
%     framenumbers=framenumbers:nframes;
% elseif isempty(framenumbers)            % if no framelist specified
%     framenumbers=1:nframes;             % read all frames.
% end

if ~exist('framenumbers','var')||isempty(framenumbers)
    framenumbers=1:nframes;
end
framenumbers(framenumbers>nframes)=[];  % remove the specified framenumbers
                                        % that exceed the max
                                        
% make cells that contain the ROI size and frame size
[x1,x2]=meshgrid(rect(2):rect(4),rect(1):rect(3));
cellinds=cell(1,numel(framenumbers));
fshape=cell(1,numel(framenumbers));
cellinds(:)={(sub2ind(wsize,x2(:),x1(:)))};
fshape(:)={([rect(3)-rect(1)+1,rect(4)-rect(2)+1])};

frames=cell(1,numel(framenumbers));frames(:)={0};   % initialize variables
framestemp=cell(1,bucketsize); framestemp(:)={0};
counter=0; counter2=0;
for ii=framenumbers
    counter=counter+1;
    fseek(fid,prod(wsize)*(ii-1)*2,'bof');          % move to image location
    framestemp{counter}=fread(fid,wsize','uint16'); % read image
    
    if counter==bucketsize                          % dump the bucket
        frames(counter2*bucketsize+1:(counter2+1)*bucketsize)=cellfun(...
            @(x,y) x(y),framestemp,cellinds(1:bucketsize),'uniformoutput',false);
        counter=0;counter2=counter2+1;              % increase bucket number
    end
end; fclose(fid);                                   % close file

% dump the last semi-full bucket
frames(counter2*bucketsize+1:counter2*bucketsize+counter)=...
    cellfun(@(x,y) x(y),framestemp(1:counter),...
    cellinds(1:counter),'uniformoutput',false);

% reshape the images
frames=cellfun(@(x,y) reshape(x,y),frames,fshape,'uniformoutput',false);

% change to unsigned integers
frames=uint16(cat(3,frames{:}));

wsize=wsize';
end