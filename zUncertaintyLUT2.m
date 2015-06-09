function zuncLUT=zUncertaintyLUT2(savedgoods,defocusing_param,...
    max_allowed_D,pxsize,z_uncertainty_curve_interval)

% This is a part of the 'Z_calibrate_astig.m' code which fits calibration
% data and generates calibration curves for 3D single-molecule imaging
% (astigmatism-based). Because it is difficult to propagate the
% localization uncertainties of z-positions in real imaging experiments, we
% use the localization uncertainties of calibration data at different
% (piezo-stage controlled) z-positions to approximate uncertainties
% associated with z-positions in real imaging experiments.

% Specifically in this code, we use the widths of x- and y-fits (Wx, Wy)
% from calibration samples to determine their z-positions (despite already
% knowing them from piezo-stage readouts) using the defocusing curve that
% is generated from the same set of calibration data. We then calculate the
% standard deviation of the distribution of z-positions for fits from
% the same known z-position. We repeats this procedure for fits from
% different (piezo-stage controlled) z-positions. Eventually we'll have a
% z-uncertainty VS z-position (sampled and interpolated) profile.

% The LUT produced by this code will be used by the
% 'find_z_uncertainty_from_LUT.m' code in the 'Astig3D_master_program.m'
% master code. Because we do not do extrapolation in the LUT, any
% z-localization falling out of the calibration range will be assigned a
% z-position-uncertainty that is equal to the first or the last value in
% the LUT. This, however, should not be a problem if implemented in the
% 'Astig3D_master_program.m' code because bad fits (e.g., those w/ x- and
% y-widths too far from the 3D calibration curve) will be removed by the
% 'determine_goodfits.m' code.

% INPUTS:

% calib_fits_path: FULL path of the 'calibration goodfits.dat' file
% returned by the 'Z_calibrate_astig.m' code earlier.

% defocusing_param_path: FULL path of the defocusing curve parameters
% returned by the 'Z_calibrate_astig.m' code earlier.

% max_allowed_D: The maximum allowable separation of an input width pair
% and the 2 calibration curves. Fits falling too far from the calibration
% curve will be rejected and have NaN for the z-center value. It is highly
% recommended to use the same number here and in cell imaging experiments
% where the defocusing curve is searched for z-positions of real signals.

% pxsize: Pixel size (nm/px).

% z_uncertainty_curve_interval: (unit: nm) The desired z-sampling-inverval
% of the output "uncertainty VS z-position" profile (lookup table). The
% uncertainties associated with z-positions not experimentally sampled are
% interpolated linearly using exising uncertainty values at experimentally
% sampled z-positions. Set to NaN to skip interpolation.

% OUTPUTS:

% z_std_LUT: A lookup table for localization uncertainties (defined as the
% standard deviation of fitted z-positions at a known z-position, unit: nm)
% at various z-positions. It is a 2-by-n matrix, where the 1st row stores
% the sampled z-positions and the 2nd row stores the (interpolated)
% z-uncertainties. The file name is '[save path]/Z-Uncertainty LUT.mat',
% where the [save path] is the same path under which the calibration
% goodfits file (and other 3D calibration files) are saved.

% A figure plotting z-uncertainty VS z-position is also saved under the
% same path and named as '[save path]/Z-Uncertainty Profile.fig'.

% The 'calib_fits' variable contains the following fields: .fitparam;
% .focal_depth; .center_w; defocusing_param.central_z_readouts and
% .z_range. Lengths are in unit of nm.
calib_fits=savedgoods;

%% ------------------------------------------------------------------------
%  Calculate the z-positions (and stdev) of fits from each known z-readout
%  ------------------------------------------------------------------------

z_readouts=unique(calib_fits(:,2));

% Pre-allocate a 2-by-n matrix to store the stdev of z-positions, where n
% is the number of unique z-readouts given by the piezo-stage, i.e., the
% z-positions where calibration is done. The 1st row stores the z-readouts,
% and 2nd row stores the standard deviation of z-positions obtained from
% curve-searching.
z_std_nm=zeros(2,length(z_readouts));

% Unique piezo-stage readouts.
z_std_nm(1,:)=unique(calib_fits(:,2))*1000;
z_std_nm_col_ind=1;

% figure, % FOR DEBUGGING db1
% Loop through each piezo-stage z-readout
for curr_z_readout=z_readouts'
    % Calibration data of fits from the current z-readout.
    curr_calib_data=calib_fits(calib_fits(:,2)==curr_z_readout,:);
        
    curr_Wx_nm=curr_calib_data(:,7)*pxsize;
    curr_Wy_nm=curr_calib_data(:,16)*pxsize;
    
    % Search the calibration curve to find z-positions
    z_nm=find_astigZ_position2(curr_Wx_nm,curr_Wy_nm,defocusing_param,...
        max_allowed_D,1);
    
    % plot(repmat(curr_z_readout, [1, length(z_nm)]), z_nm, 'x'), hold all %
    % % FOR DEBUGGING db1
    
    z_std_nm(2,z_std_nm_col_ind)=nanstd(z_nm);
    z_std_nm_col_ind=z_std_nm_col_ind+1;
end

%% ------------------------------------------------------------------------
%  Interpolate z-uncertainties at specified z-positions
%  ------------------------------------------------------------------------
z_rel_range=z_std_nm(1,:)-defocusing_param.central_z_readouts;
% The z-position values w.r.t. the central focal plane. (i.e., not the
% original piezo-readouts)

if ~isnan(z_uncertainty_curve_interval)
    z_std_LUT=defocusing_param.z_range(1):z_uncertainty_curve_interval:...
        defocusing_param.z_range(2);
    
    z_std_LUT(2,:)=interp1(z_rel_range,z_std_nm(2,:),z_std_LUT(1,:));
    
else % Skip interpolation
    z_std_LUT(1,:)=z_rel_range;
    z_std_LUT(2,:)=z_std_nm(2,:);
end
%% ------------------------------------------------------------------------
%  Generate graphical output (z-uncertainty VS z-position)
%  ------------------------------------------------------------------------
z_uncertainty_fig=figure;
plot(z_rel_range,z_std_nm(2,:),'ro','MarkerSize',5,'MarkerFaceColor',...
    'r','MarkerEdgeColor','none','LineWidth',2);

% Plot the uncertainty VS z-position curve, containing only data points
% from experimentally sampled z-positions (i.e., not interpolation here).
xlabel('z-position (nm)','Fontsize',26);
ylabel('Localization Uncertainty (nm)','FontSize',26);
set(gca,'FontSize',22,'LineWidth',2)

% Plot the uncertainty VS z-position profile for interpolated z-positions.
hold all, plot (z_std_LUT(1,:), z_std_LUT(2,:), 'kx');
legend('Sampled','Interpolated')

%  Save outputs
zuncLUT=z_std_LUT;
saveas(z_uncertainty_fig,'z-uncertainty Profile.fig');
end