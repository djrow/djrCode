function z_nm=find_astigZ_position2(Wx,Wy,defocusing_param,max_allowed_D,indRefr_corr)

% This is a part of the 'Astig3D_master_program.m' code and the
% 'zUncertaintyLUT.m' code. Based on the input x-width and y-width (nm)
% and the input astigmatism defocusing parameters, it returns the
% corresponding z-locations for the fits in nm.

% INPUTS:

% Wx: x-width (nm, sigma of Gaussian, not FWHM) of fits. It can be an
% n-by-1 vector where n is the number of fits.

% Wy: y-width (nm) of fits. Same dimension as Wx.

% defocusing_param: Defocusing curve parameters returned by the
% 'Z_calibrate_astig.m' code.

% max_allowed_D: The maximum allowable separation of an input width pair
% and the 2 calibration curves. Fits falling too far from the calibration
% curve will be rejected and have NaN for the z-center value.

% indRefr_corr: Correction factor for index of refraction mismatch. Usually
% less than 1.

% OUTPUTS:

% z_nm: An n-by-1 vector storing the Z-distance (nm) from the central focal
% plane for the fits. The Z values are interpolated/extrapolated from the
% defocusing curve linearly.

%-------------------------------------------------------------------------%
if isempty(defocusing_param)
    z_nm=nan;
    return
end

xparam=defocusing_param.fitparam(1,:);
yparam=defocusing_param.fitparam(2,:);
% d=defocusing_param.d;
w0=defocusing_param.center_w;

% Make a defocusing curve for each direction (x or y) with the input
% parameters that has a z-sampling interval of 1 nm. This curve will be
% searched subsequently to find the z-location of fits based on the input
% widths.
z_sample=defocusing_param.z_range(1):defocusing_param.z_range(2);

astigmatism_defocusing_eq=@(p,z,w)...
    w*sqrt(1+((z-p(1))/p(2)).^2+p(3)*((z-p(1))/p(2)).^3+p(4)*((z-p(1))/p(2)).^4);

wx_calib=astigmatism_defocusing_eq(xparam,z_sample,w0(1));
wy_calib=astigmatism_defocusing_eq(yparam,z_sample,w0(2));

rato=wx_calib./wy_calib;
if mean(double(diff(rato)>0))>.5
    ind=logical([diff(rato)>0,0]);
elseif mean(double(diff(rato)<0))>.5
    ind=logical([diff(rato)<0,0]);
end
z_sample=z_sample(ind);
rato=rato(ind);

z_nm=interp1(rato,z_sample,Wx./Wy)*indRefr_corr;
% z_nm(z_nm>max_allowed_D)=nan;

% % Minimize the distance between the input widths and the calibration
% % curve in the sqrt(Wx)-sqrt(Wy) space.
% [sqrt_wx_calib_grid,sqrt_Wx_grid]=meshgrid(sqrt(wx_calib),sqrt(Wx));
% [sqrt_wy_calib_grid,sqrt_Wy_grid]=meshgrid(sqrt(wy_calib),sqrt(Wy));
% 
% sqrt_Dx=sqrt_wx_calib_grid-sqrt_Wx_grid;
% sqrt_Dy=sqrt_wy_calib_grid-sqrt_Wy_grid;
% 
% D=sqrt((sqrt_Dx).^2+(sqrt_Dy).^2);
% 
% [min_D,min_D_ind]=min(D,[],2);
% 
% z_nm=z_sample(min_D_ind)*indRefr_corr;
% 
% % Input fit with widths lying too far from the calibration curve will have
% % z-position as NaN.
% z_nm(min_D>max_allowed_D) = nan;
end