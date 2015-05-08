function writebin(filename,imstack)
% inputs are strings of the full location (with filename and extension) of
% an nd2 file. in that location a file with the video data
% will be saved as a bin file with the same name. memalloc is the maximum
% number of megabytes of system memory in megabytes that you want to use.
% default is 100.
%
% use vid=bingetframes(saveloc,frames,roi); to get a particular selection
% of frames out of the saved binary file.
%
% if filename is an array of images, then it just writes that to a bin file

%% write a preloaded data matrix to the hard drive
if nargin>1

    fid=fopen(filename,'W');
    
    msize=size(imstack);
    for ii=1:msize(3)
        fwrite(fid,imstack(:,:,ii),'uint16');
    end
    brit=squeeze(sum(sum(imstack)));
    brit2=brit-min(brit);
    brit2=brit2/max(brit2);
    brit2=graythresh(brit2)*(max(brit)-min(brit))+min(brit);
    
    fwrite(fid,brit2,'double');
    fwrite(fid,msize(1:2),'double');
%     fwrite(fid,nprevious+msize(3),'double');
    fwrite(fid,msize(3),'double');
        
    fclose(fid);
    return
end

%% write a nd2 or tiff file to the hard drive
saveloc=[filename(1:end-4) '.bin'];
fid=fopen(saveloc,'W');             % make and open the file to be saved to
if strcmp(filename(end-2:end),'nd2')
%     vidid=bfGetReader(filename);
    vidid=bfGetReader(filename);        % this is an accessory program you'll need to
                                        % have in your matlab path
    display(['writing the file named ' saveloc])
    for ii=1:1e6                        % ain't no one got more than a million frames
        try                             % get the image frame
            frame=bfGetPlane(vidid,ii);
            brit(ii)=squeeze(sum(sum(frame)));
        catch
            display([filename ' done'])
            break                       % stop if there are none left
        end
        fwrite(fid,frame,'uint16');     % write the image to the growing binary file
    end
    framesize=size(frame);
    framenumber=ii-1;                   % index of last frame
else % is a tiff or tif
    info=imfinfo(filename);
    display(['writing the file named ' saveloc])
    framenumber=numel(info);
    for ii=1:framenumber
        frame=imread(filename,ii,'Info',info);
        brit(ii)=squeeze(sum(sum(frame)));
        
        fwrite(fid,frame,'uint16');     % write the image to the growing binary file
    end
    framesize=size(frame);
end

brit2=brit-min(brit);
brit2=brit2/max(brit2);
brit2=graythresh(brit2)*(max(brit)-min(brit))+min(brit);

fwrite(fid,brit2,'double');
fwrite(fid,framesize,'double');     % xwidth is at -24 bytes from end of file, ywidth is at -16
fwrite(fid,framenumber,'double');   % number of frames is at -8 from end of file

fclose(fid);

% this is how to get the frame size and number of frames

% fid=fopen(saveloc);
% fseek(fid,-24,'eof');
% framesize=fread(fid,2,'double');
% nframes=fread(fid,1,'double');

end