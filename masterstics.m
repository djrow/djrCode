function masterstics(checkvals,mainfold)

if ~exist('mainfold','var')
    mainfold=pwd;
end

% run stics?
yesstics=0;

% filters: 1: 2: 3: 4: 5: non-overlapping
sticsfilters=[];

% select cells from phase mask?
phasemasks=1;

% skip the selection of cells from the phase mask
skipselect=0;

% parameters for phase mask finding. dilate factor, low thresh, high thres,
% autofill, min area, max area
phaseparams=[2,0,.9,1,500,1000];

% Parameters for peak guessing, in the format of [noise size, particle
% size, Intensity Threshold, H-Max]. usually [1,10,2e3,1e4]
peak_guessing_params=[1,10,5000,10000];

% Minimal separation of peaks (px). Putative peaks that are closer than
% this value will be discarded unless it is the brightest one compared to
% its neighbors. Also, this is the half box size of the fitting region,
% i.e., if 'min_sep = 10', pixel intensities within a 21-by-21 square
% region will be used in PSF fitting.
min_sep=10;

% Lower and upper bounds (in pixels, except for the aspect ratio) for
% various fit parameters that will be used to determine whether a fit is
% "good". See comments of the 'determine_goodfits.m' code for more details.
width_lb=1;
width_ub=15;
aspect_ratio_ub=5;
width_error_ub=5;

% Pixel size in nm.
pxsize=49;

% camera capture rate in 1/s
framerate=100;

% maximum frame separation to consider for msd
maxtau=10;

mpp=.049;

if yesstics==1
    display('Select the movies.')
    [datalist,dataloc,findex]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
        'Matlab Files','multiselect','on');
    
    if findex==0
        fprintf('no data selected\n')
        return
    end
    
    if ~iscell(datalist); datalist={datalist}; end
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
    [dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);
    
    % WRITE BIN FILES FOR ALL MOVIES
    for ii=1:numel(dnames);
        if strcmp(dexts{ii},'.nd2')||strcmp(dexts{ii},'.tif')||strcmp(dexts{ii},'.tiff')
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]));
        end
    end
else
    fprintf('select the analysis files.\n')
    [datalist,dataloc,findex]=uigetfile([mainfold filesep '*.mat'],...
        'Matlab Files','multiselect','on');
    
    if findex==0
        fprintf('no data selected\n')
        return
    end
    
    if ~iscell(datalist); datalist={datalist}; end
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
    [dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);
end

% WRITE PHASEMASKS FILE FOR ALL MOVIES
if phasemasks==1&&yesstics==1
    display('Select the phase images.')
    [phaselist,phaselistloc,findex]=uigetfile([mainfold filesep...
        '*.nd2;*.tif*;*.bin'],'Matlab Files','multiselect','on');
    if findex==0
        fprintf('no phase images selected. try again.\n')
        return
    end
    
    if ~iscell(phaselist); phaselist={phaselist}; end
    for ii=1:numel(phaselist); phaselist{ii}=[phaselistloc phaselist{ii}]; end
    [plocs,pnames,pexts]=cellfun(@fileparts,phaselist,'uniformoutput',false);
    
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        if strcmp(pexts{ii},'.nd2')||strcmp(pexts{ii},'.tif')||strcmp(pexts{ii},'.tiff')
            writebin(phaselist{ii});
            img=bingetframes(fullfile(plocs{ii},pnames{ii},'.bin'));
        else
            img=bingetframes(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
        end
        
        mnamelist=who(m); mnamelist=mnamelist(:);
        if any(cellfun(@strcmp,mnamelist,repmat({'phaseparams'},...
                [numel(mnamelist),1])))
            phaseparams=m.phaseparams;
            
            yn=input(['input new phase image parameters? enter for ''no'', '...
                'the parameters for ''yes''.\n']);
        end
        
        if exist('yn','var')&&~isempty(yn)
            phaseparams=yn;
        end
        
        fprintf(['phasing file named: ' dnames{ii} '.\n'])
        phasemask=valley(img(:,:,1),phaseparams,checkvals);
        
        if skipselect==0
            [~,goodcells]=select_cells(phasemask);
            phasemask(logical(-1*(ismember(phasemask,goodcells)-1)))=0;
        end
        m.phasemask=phasemask;
        m.phaseparams=phaseparams;
    end
elseif yesstics==1
    for ii=1:numel(dnames)
        [~,~,sz]=bingetframes([fullfile(dlocs{ii},dnames{ii}),'.bin'],1,[]);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phasemask=ones(sz);
        m.phasemask=phasemask;
    end
end

for jj=1:numel(dnames)  % Loop each movie for sticsing
    if yesstics==0
        break
    end
    
    % Display movie folder counter
    fprintf(['This movie: ',dnames{jj},'\n'])
    
    m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat'],'Writable',true);
    phmask=m.phasemask;
    ncells=unique(phmask); ncells(ncells==0)=[];
    
    for ii=ncells(:)'
        if yesstics==0
            break
        end
        
        sphmask=phmask==ii;
        
        roiinds=zeros(sum(double(sphmask(:))),2);
        [roiinds(:,1),roiinds(:,2)]=find(sphmask);
        roi=[min(roiinds(:,1)),min(roiinds(:,2)),...
            max(roiinds(:,1)),max(roiinds(:,2))];
        
        % use phmask to assign the correct portion of phmask to another
        % variable, pm2
        sphmask=sphmask(roi(1):roi(3),roi(2):roi(4));
        
        [~,n]=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],1,roi);
        
        % max working data set size: 1 gigabyte/(8 bytes/double)
        maxdoubles=1e9/8;
        
        % movie is quadroupled in size due to fft padding
        potsize=n*(roi(3)-roi(1))*(roi(4)-roi(2))*4;
        
        % number of movie subdivisions
        nbins=ceil(potsize/maxdoubles);
        
        % number of frames in each subdivision
        binsize=floor(n/nbins);
        
        if nbins>1
            
            timecorr=cell(1,nbins);
            for kk=1:nbins
                fnums=1+binsize*(kk-1):binsize*kk;
                v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],fnums,roi);
                timecorr{kk}=stics3(v,sphmask,10,[1,2,4]);
            end
            
            timecorr=mean(cat(4,timecorr{:}),4);
        else
            
            v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],[],roi);
            timecorr=stics3(v,sphmask,10,[1,2,4]);
        end
        
        % use cyldist to find the cell's orientation and length
        [~,p,~,l]=cyldist(roiinds);
        thet(ii)=atan(p(1)/p(2));
        leng(ii)=l-phaseparams(1)*2*mpp;
        
        % fit the series to gaussians
        for kk=1:maxtau
            pgauss(ii,kk,:)=gaussfit(timecorr(:,:,kk),156,0,-thet(ii));
        end
    end
    
    % WRITE RESULTS TO ANALYSIS FILE
    m.msds=pgauss;
    m.thet=thet;
    m.leng=leng;
    m.sticsfilters=sticsfilters;
end

% %% fit the average msd curve
for jj=1:numel(dnames)
    m=matfile(fullfile(dlocs{jj},dnames{jj}),'Writable',true);
    
    pgauss=m.msds;
    msds{jj}=pgauss(:,:,[2,3]).^2*2;
end

totsmcgoats=cat(1,msds{:});

% the ratio of these two msd curves is constant and proportional to 
end

function [cell_xy,good_cell]=select_cells(PhaseMask)

% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

% INPUTS:

% PhaseMask: Numbered phase mask of cell images returned by the 'valley.m'
% code or any other cell segmentation code

% OUTPUTS:

% cell_xy: The x- and y-coordinates clicked on by the users to select cells

% good_cell: The indices of chosen cells

cell_fig_h=figure;
imshow(PhaseMask~=0,[],'Border','Tight');
title(sprintf(['Left-click on cells to be analyzed.\n Press ''Enter'' to proceed.\n'...
    ' Or click return without clicking\n on any cell to analyze ALL of them.']))

% Now move the axes slightly so that the top of the title is visible
set(cell_fig_h,'Units','normalized')
P=get(cell_fig_h,'Position');
set(cell_fig_h,'Position',P+[0 -0.03 0 0])

conv_hull=regionprops(PhaseMask,'Convexhull');
for i=1:length(conv_hull)
    hold all,
    plot(conv_hull(i,1).ConvexHull(:,1),...
        conv_hull(i,1).ConvexHull(:,2),'c-','linewidth',2)
end

key_input=0;
% Let the user pick cells with good shapes
cell_x=[]; cell_y=[];
while key_input~=121
    [curr_cell_x,curr_cell_y,key_input]=jsbginput(1);
    hold all,
    plot(curr_cell_x,curr_cell_y,'m*','markersize',12);
    cell_x=vertcat(cell_x,curr_cell_x);  %#ok<AGROW>
    cell_y=vertcat(cell_y,curr_cell_y); %#ok<AGROW>
end

if ~isempty(cell_x) % If the user did click on something
    cell_xy=[cell_x,cell_y];
    
    good_cell=PhaseMask(sub2ind(size(PhaseMask),round(cell_y),round(cell_x)));
    good_cell=good_cell(good_cell~=0);
    if ~isempty(good_cell) % If the user did click on at least 1 cell
        good_cell=unique(good_cell);
    else % If the user only clicked on background, use all cells.
        good_cell=1:max(PhaseMask(:));
    end
else % If the user hits "Enter" directly, use all cells.
    cell_xy=[];
    good_cell=1:max(PhaseMask(:));
end
close(cell_fig_h);
end