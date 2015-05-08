function output=masterstics(checkvals,mainfold)

output.msdparams=[];
if ~exist('mainfold','var')
    mainfold=pwd;
end

% select cells from phase mask?
phasemasks=1;

% skip the selection of cells from the phase mask
skipselect=0;

% run stics?
yesstics=0;

% run gaussian fitting?
% yesgauss=1;

% filters: 1: 2: 3: 4: 5: non-overlapping
sticsfilters=[1,2,4];

% which fitting parameters for the correlation gaussian fit?
% 156 is independently variable amplitude and two widths.
% with offset 0 and specified rotation angle.
gfittype=156;

% parameters for phase mask finding. dilate factor, low thresh, high thres,
% autofill, min area, max area
phaseparams=[2,0,.9,1,500,1000];

% maximum frame separation to consider for msd
maxtau=10;

mpp=.049;
% inttime=.04;

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
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}(1:end-13)]; end
    [dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);
end

% WRITE PHASEMASKS FILE FOR ALL MOVIES
if phasemasks==1&&yesstics==1&&checkvals==1
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
        
        fprintf(['phasing file named: ' dnames{ii} '.\n'])
        counter=0;
        while counter<2
            
            mnamelist=who(m); mnamelist=mnamelist(:);
            if any(cellfun(@strcmp,mnamelist,repmat({'phaseparams'},...
                    [numel(mnamelist),1])))&&counter==0
                phaseparams=m.phaseparams;
            end
            
            if exist('yn','var')&&~isempty(yn)
                phaseparams=yn;
                counter=0;
            else
                counter=counter+1;
            end
            
            phasemask=valley(img(:,:,1),phaseparams,checkvals);
            
            ppar.dilate_factor=phaseparams(1);
            ppar.low_threshold=phaseparams(2);
            ppar.high_threshold=phaseparams(3);
            ppar.autofill_bool=phaseparams(4);
            ppar.min_area=phaseparams(5);
            ppar.max_area=phaseparams(6);
            
            display(ppar)
            yn=input(['input new phase image parameters?\n enter for ''no'', '...
                'the parameters for ''yes''. \n two sequential ''nos'' means '...
                'it''s good to go.\n']);
        end
        if skipselect==0
            [~,goodcells]=select_cells(phasemask,img);
            phasemask(logical(-1*(ismember(phasemask,goodcells)-1)))=0;
        end
        m.phasemask=phasemask;
        m.phaseparams=phaseparams;
    end
elseif yesstics==1&&checkvals==1
    for ii=1:numel(dnames)
        [~,~,sz]=bingetframes([fullfile(dlocs{ii},dnames{ii}),'.bin'],1,[]);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phasemask=ones(sz);
        m.phasemask=phasemask;
    end
end

for jj=1:numel(dnames)  % Loop each movie for sticsing
    if checkvals==1
        return          % no need to run any more code
    end
    if yesstics==0
        break           % continue to the fitting shenanigans
    end
    
    % Display movie folder counter
    fprintf(['Working on this movie: ',dnames{jj},'\n'])
    
    m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat'],'Writable',true);
    phmask=m.phasemask;
    ncells=unique(phmask); ncells(ncells==0)=[];
    
    timecorrsave=cell(1,numel(ncells));
    for ii=ncells(:)'
        
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
        
        fprintf(['cross-correlating cell ', num2str(ii), ' of ' num2str(max(ncells)), '.\n'])
        if nbins>1
            
            timecorr=cell(1,nbins);
            for kk=1:nbins
                fnums=1+binsize*(kk-1):binsize*kk;
                v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],fnums,roi);
                timecorr{kk}=stics3(v,sphmask,10,sticsfilters);
            end
            
            timecorr=mean(cat(4,timecorr{:}),4);
        else
            
            v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],[],roi);
            timecorr=stics3(v,sphmask,10,[1,2,4]);
        end
        
        % use cyldist to find the cell's orientation and length
        [~,p,~,l]=cyldist(roiinds);
        thet(ii)=-atan(p(1)/p(2));
        leng(ii)=l-phaseparams(1)*2*mpp;
        
        timecorrsave{ii}=timecorr;
    end
    
    % WRITE RESULTS TO ANALYSIS FILE
    m.tcorr=timecorrsave;
    m.thet=thet;
    m.leng=leng;
    m.sticsfilters=sticsfilters;
end

for ii=1:numel(dnames)                % loop through movies
    m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
    
    tcorr=m.tcorr;
    thet=m.thet;
    
    pgauss=zeros(numel(tcorr),size(tcorr{1},3),numel(num2str(gfittype)));
    for jj=1:numel(tcorr)             % loop through cells in a movie
        workingcorr=tcorr{jj};
        if isempty(workingcorr)
            continue
        end
        for kk=1:size(workingcorr,3)  % loop through time lags
            [pgauss(jj,kk,:),~,nrmse(jj,kk)]=gaussfit(workingcorr(:,:,kk),gfittype,0,thet(jj));
        end
    end
    
    m.pgauss=pgauss;
    m.nrmse=nrmse;
end

%% fit the average msd curve
fprintf('fitting gaussians.\n')
msds=cell(1,numel(dnames));
for jj=1:numel(dnames)
    m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat']);
    pgauss=m.pgauss;
    nrmse=m.nrmse;
    
    missingcells=all(all(pgauss==0,2),3);
    pgauss(find(missingcells),:,:)=[];
    
    msds{jj}=pgauss(:,:,[2,3]).^2*2*mpp^2;
end
meanmsds=squeeze(mean(cat(1,msds{:})));

inttime=input('enter the integration time in seconds: ');
tau=(1:maxtau)*inttime;
cs=[.5,1];
for ii=1:2
    pstart=[cs(ii),1,.0384];
    mdl=fitnlm(tau,meanmsds(1:maxtau,ii),@confmodel,pstart);
    msdparams(ii,:)=abs(mdl.Coefficients{:,1});
    
    subplot(1,2,ii)
    plot(tau,meanmsds(1:maxtau,ii));
    hold all
    plot(tau,confmodel(msdparams(ii,:),tau))
    hold off
end

output.msdparams=msdparams;
end

function z=confmodel(p,x)
% global camerasd
d=abs(p(2)); l=p(1);
tau=x;

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));% counter=0;
for ii=1:2:2*200-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
% display(counter)
z=l^2/6*(1-96/pi^4*temp)+p(3); %.0384;
end

function [cell_xy,good_cell]=select_cells(PhaseMask,im)

% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

% INPUTS:

% PhaseMask: Numbered phase mask of cell images returned by the 'valley.m'
% code or any other cell segmentation code

% OUTPUTS:

% cell_xy: The x- and y-coordinates clicked on by the users to select cells

% good_cell: The indices of chosen cells

cell_fig_h=figure;
subplot(121)
imshow(im,[])
subplot(122)
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