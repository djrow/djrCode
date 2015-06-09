function [cali_data,zuncLUT]=Z_calibrate_astig2(calibfolder,checkparams)
% calibfolder is the location of the folder with the nd2 movies of
% stationary particles at known z positions. if this function is run
% multiple times, the particles will not be fit, but the defocusing
% parameters and z-LUT will be remade using the fits saved in calibdata.mat

%% ------------------------------------------------------------------------
%  User-defined parameters
%  ------------------------------------------------------------------------
forcefit=0;

origdir=pwd;
cd(calibfolder)
z_profile = 22.10:0.05:23.10; % start z: z-increment: end-z

% Parameters for peak guessing, in the format of [LP, HP, Intensity
% Threshold, H-Max].
peak_guessing_params=[1,10,5e2,1e3];

% Minimal separation of peaks (px). Putative peaks that are closer than
% this value will be discarded unless it is the brightest one compared to
% its neighbors. Also, this is the half box size of the fitting region,
% i.e., if 'min_sep = 10', pixel intensities within a 21-by-21 square
% region will be used in PSF fitting.
min_sep=15;

% Width threshold in determined whether or not a data point is to be used
% in generating the calibration curve. The input should be a 4-member
% vector specifying [lower width threshold, upper width threshold, aspect
% ratio upper threshold, width error upper threshold].
width_T=[1,15,4,1]; % Unit: pixels

% A 4-element vector specifying the initial guesses for the parameters of
% the defocusing equation: x-offset, y-offset (both in unit of nm), the
% 3rd-order correction and the 4th-order correction (both unit-less).
defocusing_param_guess=[-200,100,500,.1,.05];

d=1000; % Focal depth of microscope (nm)

% (Used in estimation of z-position uncertainties, unit: nm) The maximum
% allowable separation of an input width pair and the 2 calibration curves.
% Fits falling too far from the calibration curve will be rejected and have
% NaN for the z-center value. It is highly recommended to use the same
% number here and in cell imaging experiments where the defocusing curve is
% searched for z-positions of real signals.
max_allowed_D=100;

pxsize=49; % Pixel size (nm/px)

% (Used in estimation of z-position uncertainties, unit: nm) The desired
% z-sampling-inverval of the "uncertainty VS z-position" profile (lookup
% table) returned by the 'zUncertaintyLUT.m' code. The uncertainties
% associated with z-positions not experimentally sampled are interpolated
% linearly using exising uncertainty values at experimentally sampled
% (during calibration) z-positions. Set to NaN to skip interpolation.
z_uncertainty_curve_interval=10;

%% ------------------------------------------------------------------------
%  Fitting PSF
%  ------------------------------------------------------------------------

nd2names=dir('*.nd2');
for ii=1:numel(nd2names)
    ie=dir([nd2names(ii).name(1:end-4) '.bin']);
    if numel(ie)~=1
        writebin(nd2names(ii).name)
    end
end

binnames=dir('*.bin');
goodfits=[]; allfits=[];
m=matfile('calibdata.mat','writable',true);
mnamelist=who(m); mnamelist=mnamelist(:);
if any(cellfun(@strcmp,mnamelist,repmat({'savedgoods'},...
        [numel(mnamelist),1])))==0||forcefit==1
    h1=waitbar(0,'fitting stuff');
    for ii=1:numel(binnames) % Loop through each movie
        waitbar(ii/numel(binnames),h1);
        curr_z=z_profile(ii);
        [vid,num_fr]=bingetframes(binnames(ii).name,[],[]);
        if checkparams==1
            frameskip=[];
            for jj=1:num_fr % Loop through each frame
                [all_fitparam,all_fiterr,~,frameskip]=fit_asym_gauss_int2(vid(:,:,jj),...
                    ones(size(vid(:,:,jj))),2,peak_guessing_params,min_sep,frameskip,jj);
                
                % # of total fits from current frame (good or bad)
                numfits=size(all_fitparam,1);
                
                allfitdata=[repmat(jj,[numfits,1]),repmat(curr_z,[numfits,1]),...
                    all_fitparam(:,1),all_fiterr(:,1),all_fitparam(:,2),all_fiterr(:,2),...
                    all_fitparam(:,3),all_fiterr(:,3),all_fitparam(:,4),all_fiterr(:,4),...
                    all_fitparam(:,5),all_fiterr(:,5),nan(numfits,1),...
                    all_fitparam(:,7),repmat(min_sep*2+1,[numfits,1]),...
                    all_fitparam(:, 6),all_fiterr(:,6),all_fitparam(:,3)./all_fitparam(:,6)];
                
                % determine good fits based on the following criterion: (1)
                % Amplitude > Amplitude error, (2) X- and Y-widths both larger than their respective
                % statistical errors, (3) Statistical errors of X- and Y-widths are both
                % smaller than a predefined value (4) X- and Y-widths are both within a
                % preset range (not too narrow or too wide) and (5) Aspect ratio (Wx / Wy
                % and Wy / Wx) is within a preset range.
                allfitdata(:,13)=(allfitdata(:,3)>allfitdata(:,4)).*...
                    (allfitdata(:,7)>allfitdata(:,8)).*(allfitdata(:,16)>allfitdata(:,17))...
                    .*(allfitdata(:,8)<width_T(4)).*(allfitdata(:,17)<width_T(4))...
                    .*(allfitdata(:,7)>width_T(1)).*(allfitdata(:,7)<width_T(2))...
                    .*(allfitdata(:,16)>width_T(1)).*(allfitdata(:,16)<width_T(2))...
                    .*(allfitdata(:,18)<width_T(3)).*((1/allfitdata(:,18))<width_T(3))';
                
                % All good fits from the CURRENT FRAME.
                goodfitdata=allfitdata(allfitdata(:,13)==1,:);
                
                savedgoods{ii}{jj}=goodfitdata;
                savedall{ii}{jj}=allfitdata;
            end
        else
            parfor jj=1:num_fr % Loop through each frame
                [all_fitparam,all_fiterr]=fit_asym_gauss_int2(vid(:,:,jj),...
                    ones(size(vid(:,:,jj))),2,peak_guessing_params,min_sep);
                numfits=size(all_fitparam,1);
                allfitdata=[repmat(jj,[numfits,1]),repmat(curr_z,[numfits,1]),...
                    all_fitparam(:,1),all_fiterr(:,1),all_fitparam(:,2),all_fiterr(:,2),...
                    all_fitparam(:,3),all_fiterr(:,3),all_fitparam(:,4),all_fiterr(:,4),...
                    all_fitparam(:,5),all_fiterr(:,5),nan(numfits,1),...
                    all_fitparam(:,7),repmat(min_sep*2+1,[numfits,1]),...
                    all_fitparam(:, 6),all_fiterr(:,6),all_fitparam(:,3)./all_fitparam(:,6)];
                allfitdata(:,13)=(allfitdata(:,3)>allfitdata(:,4)).*...
                    (allfitdata(:,7)>allfitdata(:,8)).*(allfitdata(:,16)>allfitdata(:,17))...
                    .*(allfitdata(:,8)<width_T(4)).*(allfitdata(:,17)<width_T(4))...
                    .*(allfitdata(:,7)>width_T(1)).*(allfitdata(:,7)<width_T(2))...
                    .*(allfitdata(:,16)>width_T(1)).*(allfitdata(:,16)<width_T(2))...
                    .*(allfitdata(:,18)<width_T(3)).*((1/allfitdata(:,18))<width_T(3))';
                goodfitdata=allfitdata(allfitdata(:,13)==1,:);
                temp{jj}=goodfitdata;
                temp2{jj}=allfitdata;
            end
            savedall{ii}=temp;
            savedgoods{ii}=temp2;
        end
    end
    savedgoods=cat(2,savedgoods{:});
    savedgoods=cat(1,savedgoods{:});
    savedall=cat(2,savedall{:});
    savedall=cat(1,savedall{:});    
    m.savedgoods=savedgoods;
    m.savedall=savedall;
    
    if ishandle(h1)
        close(h1)
    end
else
    savedgoods=m.savedgoods;
    savedall=m.savedall;
end

%% ------------------------------------------------------------------------
%  Compile X- and Y-widths and aspect ratio versus z-readout
%  ------------------------------------------------------------------------

% Determine the average X- and Y-widths and its associated standard error
% for each unique Z readout.
mean_x_width=zeros(1,length(z_profile));
mean_y_width=zeros(1,length(z_profile));
stderr_x_width=zeros(1,length(z_profile));
stderr_y_width=zeros(1,length(z_profile));

for ii=1:numel(z_profile)
    curr_z=z_profile(ii);
    curr_z_data=savedgoods(abs(savedgoods(:,2)-curr_z)<0.00001,:);
    curr_z_xwidth=curr_z_data(:,7);
    curr_z_ywidth=curr_z_data(:,16);
    num_pts=size(curr_z_data,1);
    
    %     mean_x_width(ii)=mean(curr_z_xwidth)*pxsize;
    %     mean_y_width(ii)=mean(curr_z_ywidth)*pxsize;
    
    curr_z_xwidtherr=curr_z_data(:,8);
    curr_z_ywidtherr=curr_z_data(:,17);
    mean_x_width(ii)=sum(curr_z_xwidth./curr_z_xwidtherr.^2)/...
        sum(1./curr_z_xwidtherr.^2)*pxsize;
    mean_y_width(ii)=sum(curr_z_ywidth./curr_z_ywidtherr.^2)/...
        sum(1./curr_z_ywidtherr.^2)*pxsize;
    
    stderr_x_width(ii)=pxsize*std(curr_z_xwidth)/num_pts;
    stderr_y_width(ii)=pxsize*std(curr_z_ywidth)/num_pts;
end

% Compile calibration data and save it as a matrix file
% The calibration data is saved into a 5-row matrix following this format:
% [z-readouts; mean_x_width; mean_y_width; stderr_x_width; stderr_y_width];
% All expressed in nm!
cali_data=[z_profile*1000;mean_x_width;mean_y_width;stderr_x_width;...
    stderr_y_width];

% Export the defocusing curve fit parameters (stored as a 2-by-3 matrix in
% a structure called 'defocusing_param.fitparam') as a .mat file into the
% same folder where the calibration data is saved. The 3 columns of the
% 'defocusing_param.fitparam' store the [focal-plane-offset (nm), 3rd-order
% correction, 4th-order correction] respectively. Also exported are the
% focal depth used in the fitting (as another field called
% 'defocusing_param.focal_depth'), the width at the central focal plane (as
% 'defocusing_param.center_w'), the original piezo-stage readout (nm) for
% the central plane (as 'defocusing_param.fit_defocusing_eq') and the
% 3-member vector storing the z-sampling range (min, max and interval, as
% 'defocusing_param.z_range'). figure saved called calibration_fit.fig
defocusing_param=fit_defocusing_eq2(cali_data,defocusing_param_guess,d,'all');

%  Estimate localization uncertainties in z-positions. figure saved called
%  z-LUT something .fig
zuncLUT=zUncertaintyLUT2(savedgoods,defocusing_param, ...
    max_allowed_D,pxsize,z_uncertainty_curve_interval);

m.cali_data=cali_data;
m.defocusing_param=defocusing_param;
m.zuncLUT=zuncLUT;
cd(origdir)
end