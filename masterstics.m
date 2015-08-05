function output=masterstics(checkvals,mainfold)

output.msdp=[];

if ~exist('mainfold','var')
    mainfold=pwd;
end

% select cells from phase mask?
phasemasks=0;

% skip the selection of cells from the phase mask
skipselect=0;

% run stics?
yesstics=0;

% run gaussian fitting?
yesgauss=0;

% show gaussian fits?
showgauss=1;

% filters: 1: 2: 3: 4: non-overlapping time lags; 5:
sticsfilters=[1,3];

% which fitting parameters for the correlation gaussian fit?
% 156 is independently variable amplitude and two widths.
% with offset 0 and specified rotation angle.
gfittype=123456;
widths=find(num2str(gfittype)=='5'|num2str(gfittype)=='6');

% parameters for phase mask finding. dilate factor, low thresh, high thres,
% autofill, min area, max area
phaseparams=[2,0,.9,1,500,1000];

% maximum frame separation to consider for msd
maxtau=20;

mpp=.049;
% inttime=.04;

%% locate and prepare data
[datalist,dataloc,findex]=uigetfile([mainfold filesep ...
    '*.nd2;*.tif*;*.bin;*.xls*'],...
    'Matlab Files','multiselect','on');

if findex==0
    fprintf('no data selected\n')
    return
end

if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);

infomat=zeros(numel(dlocs),5);
if all(cellfun(@(x)~isempty(regexp(x,regexptranslate('wildcard','*xls*'), 'once')),dexts))
    [infomat,textdat]=xlsread(fullfile(dlocs{1},[dnames{1},dexts{ii}]));
    textdat=textdat(2:end,1:3);
    
    dlocs=textdat(:,2);
    dnames=textdat(:,1);
end
if yesstics==1
    % WRITE BIN FILES FOR ALL MOVIES
    for ii=1:numel(dnames);
        if strcmp(dexts,'.nd2')||strcmp(dexts,'.tif')||strcmp(dexts,'.tiff')
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts]));
        end
    end
end

%% WRITE PHASEMASKS FILE FOR ALL MOVIES
if phasemasks&&yesstics&&checkvals
    if all(cellfun(@(x) regexp(x,regexptranslate('wildcard','*xls*')),dexts))
        plocs=dlocs;
        pnames=dnames;
        pext='.tif';
        
    else
        display('Select the phase images.')
        [phaselist,phaselistloc,findex]=uigetfile([mainfold filesep...
            '*.nd2;*.tif*;*.bin'],'Matlab Files','multiselect','on');
        if findex==0
            fprintf('no phase images selected. try again.\n')
            return
        end
        
        if ~iscell(phaselist); phaselist={phaselist}; end
        for ii=1:numel(phaselist); phaselist{ii}=[phaselistloc phaselist{ii}]; end
        [plocs,pnames,pext]=cellfun(@fileparts,phaselist,'uniformoutput',false);
    end
    
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        if strcmp(pext,'.tif')||strcmp(pext,'.tiff')
            img=imread(fullfile(plocs{ii},[pnames{ii},pext]),1,...
                'Info',imfinfo(fullfile(plocs{ii},[pnames{ii},pext])));
            
        elseif strcmp(pext,'.nd2')
            vidid=bfGetReader(fullfile(plocs{ii},[pnames{ii},pext]));
            img=bfGetPlane(vidid,ii);
        end
        
        mnamelist=who(m); mnamelist=mnamelist(:);
        if any(cellfun(@strcmp,mnamelist,repmat({'phaseparams'},...
                [numel(mnamelist),1])))
            phaseparams=m.phaseparams;
        end
        
        fprintf(['phasing file named: ' dnames{ii} '.\n'])
        counter=0;
        while counter<2
            
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
        if ~skipselect
            [~,goodcells]=select_cells(phasemask,img);
            phasemask(logical(-1*(ismember(phasemask,goodcells)-1)))=0;
        end
        m.phasemask=phasemask;
        m.phaseparams=phaseparams;
    end
elseif yesstics&&checkvals
    for ii=1:numel(dnames)
        [~,~,sz]=bingetframes([fullfile(dlocs{ii},dnames{ii}),'.bin'],1,[]);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phasemask=ones(sz);
        m.phasemask=phasemask;
    end
end

if checkvals==1
    return          % no need to run any more code
end

if yesstics
    %% perform STICS
    parfor jj=1:numel(dnames)  % Loop each movie for sticsing
        
        % Display movie folder counter
        fprintf(['Working on this movie: ',dnames{jj},'\n'])
        
        m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat'],'Writable',true);
        phmask=m.phasemask;
        ncells=unique(phmask); ncells(ncells==0)=[];
        
        timecorrsave=cell(1,sum(ncells>0));
        thet=zeros(1,sum(ncells>0));
        leng=zeros(1,sum(ncells>0));
        nframes=zeros(1,sum(ncells>0));
        sphmask=cell(1,sum(ncells>0));
        counter=0;
        for ii=ncells(:)'
            counter=counter+1;
            
            sphmask{counter}=phmask==ii;
            
            roiinds=zeros(sum(double(sphmask{counter}(:))),2);
            [roiinds(:,1),roiinds(:,2)]=find(sphmask{counter});
            roi=[min(roiinds(:,1)),min(roiinds(:,2)),...
                max(roiinds(:,1)),max(roiinds(:,2))];
            
            % use phmask to assign the correct portion of phmask to another
            % variable, pm2
            sphmask{counter}=sphmask{counter}(roi(1):roi(3),roi(2):roi(4));
            
            [~,n]=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],1,roi);
            
            % max working data set size: 1 gigabyte/(8 bytes/double)
            maxdoubles=1e8/8;
            
            % movie is quadroupled in size due to fft padding
            potsize=n*(roi(3)-roi(1))*(roi(4)-roi(2))*4;
            
            % number of movie subdivisions
            nbins=ceil(potsize/maxdoubles);
            
            % number of frames in each subdivision
            binsize=floor(n/nbins);
            
            fprintf(['cross-correlating cell ', num2str(counter), ' of ' num2str(max(ncells)), '.\n'])
            if nbins>1
                
                timecorr=cell(1,nbins);
                for kk=1:nbins
                    fnums=1+binsize*(kk-1):binsize*kk;
                    v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],fnums,roi);
                    timecorr{kk}=stics3(double(v),sphmask{counter},maxtau,sticsfilters);
                end
                
                timecorr=mean(cat(4,timecorr{:}),4);
            else
                
                v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],[],roi);
                timecorr=stics3(double(v),imdilate(sphmask{counter},ones(2)),maxtau,sticsfilters);
            end
            
            % use cyldist to find the cell's orientation and length
            [~,p,~,l]=cyldist(roiinds);
            thet(counter)=-atan(p(1)/p(2));
            leng(counter)=l-phaseparams(1)*2*mpp;
            
            timecorrsave{counter}=timecorr;
            nframes(counter)=n;
        end
        
        % WRITE RESULTS TO ANALYSIS FILE
        m.tcorr=timecorrsave;
        m.thet=thet;
        m.leng=leng;
        m.sticsfilters=sticsfilters;
        m.sphmask=sphmask;
        m.nframes=nframes;
    end
end
counter=0;
%% fit the correlations
if yesgauss
    for ii=1:numel(dnames)                % loop through movies
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        fprintf(['Fitting the correlations from this movie: ',dnames{ii},'\n'])
        
        tcorr=m.tcorr;
        thet=m.thet;
        sphmask=m.sphmask;
        sphmask=sphmask(cellfun(@(x)~isempty(x),sphmask));
        
        pgauss=zeros(numel(tcorr),maxtau,numel(num2str(gfittype)));
        p95=pgauss;
        nrmse=zeros(numel(tcorr),maxtau);
        for jj=find(cellfun(@(x)~isempty(x),tcorr)) % loop through cells in a movie
            counter=counter+1;
            
            % pad the phase mask to make it the size of the correlation func
            %         phmask=m.phasemask;
            bphmask=padarray(sphmask{jj},floor((size(tcorr{jj}(:,:,1))-size(sphmask{jj}))/2));
            if size(bphmask,1)<size(tcorr{jj},1)
                bphmask=padarray(bphmask,[1,0],'pre');
            end
            if size(bphmask,2)<size(tcorr{jj},2)
                bphmask=padarray(bphmask,[0,1],'pre');
            end
            
            res=zeros([size(tcorr{jj}(:,:,1)),maxtau]);
            for kk=1:size(tcorr{jj},3)-1  % loop through time lags
                workingcorr=tcorr{jj}(:,:,kk+1);
                workingcorr(~bphmask)=nan;
                
                [pgauss(jj,kk,:),p95(jj,kk,:),nrmse(jj,kk),res(:,:,kk)]=...
                    gaussfit(workingcorr,gfittype,thet(jj));
                
            end
            ampandwidths=find(num2str(gfittype)=='1'|num2str(gfittype)=='5'|...
                num2str(gfittype)=='6');
            pgauss(jj,:,ampandwidths)=abs(pgauss(jj,:,ampandwidths));
            
            if showgauss
                im1=tcorr{jj}(:,:,2); im1(~bphmask)=nan;
                im2=tcorr{jj}(:,:,maxtau); im2(~bphmask)=nan;
                r1=res(:,:,1); r1(~bphmask)=nan;
                r2=res(:,:,maxtau); r2(~bphmask)=nan;
                
                subplot(321); pcolor(im1);
                title('data')
                axis image; shading flat
                cl=get(gca,'clim');
                
                subplot(322); pcolor(r1);
                title('residuals')
                axis image; shading flat
                cl=cl-min(cl)+nanmin(r1(:));
                set(gca,'clim',cl);
                
                subplot(323); pcolor(im2);
                title('data')
                axis image; shading flat
                cl=cl-min(cl)+nanmin(im2(:));
                set(gca,'clim',cl);
                
                subplot(324); pcolor(r2);
                title('residuals')
                axis image; shading flat
                cl=cl-min(cl)+nanmin(r2(:));
                set(gca,'clim',cl);
                
                subplot(3,2,5:6)
                scatter(1:size(pgauss,2),squeeze(pgauss(jj,:,6).^2*2*.049^2),'fill')
                axis tight
            end
        end
        
        m.pgauss=pgauss;
        m.nrmse=nrmse;
        m.p95=p95;
    end
end

%% select data

% load data
msds=cell(1,numel(dlocs));
cellLengths=cell(1,numel(dlocs));
infomatCell=cell(1,numel(dlocs));
nrmse=cell(1,numel(dlocs));
p95=cell(1,numel(dlocs));
nframes=cell(1,numel(dlocs));

fprintf('loading msd data\n')
parfor jj=1:numel(dlocs)
    m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat']);
    msds{jj}=m.pgauss.^2*2*mpp^2;
    cellLengths{jj}=m.leng;
    nrmse{jj}=m.nrmse;
    p95{jj}=m.p95;
    nframes{jj}=m.nframes;
    infomatCell{jj}=infomat(jj*ones(1,size(nframes{jj},2)),:);
    
%     missingcells=all(all(pgauss==0,2),3);
%     ltemp(missingcells)=[];
%     pgauss(find(missingcells),:,:)=[];
%     ntemp(missingcells)=[];
    
%     msds{jj}=pgauss(:,:,widths).^2*2*mpp^2;
%     cellLengths{jj}=ltemp;
%     nframes{jj}=ntemp;
end

cellLengths=cat(2,cellLengths{:})';
nrmse=cat(1,nrmse{:});
p95=cat(1,p95{:});
nframes=cat(2,nframes{:})';

infomat=cat(1,infomatCell{:});
infomat=cat(2,infomat,cellLengths);

msds=cat(1,msds{:});

%% fit msds
maxtau=15;
msdp=zeros(size(infomat,1),3,1);
dest=zeros(size(infomat,1),1);
msdp95=zeros(size(infomat,1),3,1);
fprintf('fitting msds\n')
parfor jj=1:floor(size(infomat,1))
    inttime=infomat(jj,1)*.001;
    tau=(1:maxtau)*inttime;
    pstart=[cellLengths(jj)-.5,5,.0384];
    y=msds(jj,1:maxtau,2);
    
    lb=[0,0,0]; ub=[100,200,100];
    
    %     mdl=fitnlm(tau,y,@confmodel,pstart);
    %     msdp(jj,:,1)=abs(mdl.Coefficients{:,1});
    
    [x,~,resid,~,~,~,jaco]=lsqcurvefit(@confmodel,pstart,tau,y,lb,ub);
    
    msdp95(jj,:)=diff(nlparci(x,resid,'jacobian',jaco),1,2);
    
    msdp(jj,:)=x;
    
    dest(jj,1)=(y(2)-y(1))/4/inttime;
    %     msdp95(jj,:,1)=diff(coefCI(mdl),[],2);
    
%     scatter(tau,y,'fill'); hold all
%     plot(tau,confmodel(msdp(jj,:),tau),'linewidth',2); hold off
%     plot(tau,confmodel(pstart,tau),'--'); hold off
end

% pstart=[2,5,.0384];
% mdl=fitnlm(tau(~isnan(meanmsd)),meanmsd(~isnan(meanmsd)),@confmodel,pstart);
% meanmsdp=abs(mdl.Coefficients{:,1});

output.msdp=msdp;
output.dEstimate=dest;
output.msdp95=msdp95;
output.nframes=nframes;
output.cellLengths=cellLengths;
output.inttime=infomat(:,1)*.001;
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