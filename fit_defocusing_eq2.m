function defocusing_param=...
    fit_defocusing_eq2(cali_data,guess_param,d,fit_z_range)
warning off %#ok<WNOFF>
%-------------------------------------------------------------------------%
% Fit the calibration data (x- and y-widths versus z-distance from the
% central plane) to an analytical defocusing equation which will be used to
% determine the z-location of a single-molecule emitter in an astigmatism
% 3D imaging experiment.

% INPUTS:

% calibration_data_path: The FULL path where the z-calibration data are
% stored. The data should be a 5-row matrix having the form of [z-readouts;
% sigma_x; sigma_y; stderr_x_width; stderr_y_width], all in unit of nm. The
% matrix has the variable 'cali_data'. Note the z-readouts input here is
% the original z-distance given by the piezo stage, and will be converted
% to z-distance (nm) from the central focal plane by the code
% automatically for defocusing function fitting.

% guess_param: A 4-member vector storing the guessed parameters for the
% x-plane offset (nm), the y-plane offset (nm), the 3rd order and the 4th
% order correction term (unitless).

% d: The focal depth of the microscope in nm, which is a fixed parameter in
% the defocusing curve.

% fit_z_range: A 2-member vector specifying the lower and upper range for
% z-distance to central plane in unit of nm, data obtained with z-distance
% to the central plane out of this range will not be fit. If you want to
% use all available z-distance values, specify this parameter as string 'all'.

% OUTPUTS:

% defocusing_param: A structure having the following fields:
% defocusing_param.fitparam: 2-by-3 matrix storing the fit parameters of
% the defocusing equation in the following form: [focal-plane-offset (nm),
% the 3-order correction (unitless), the 4th order correction (unitless)],
% the 1st row is for the x-defocusing curve and the 2nd row is for the
% y-defocusing curve.
% defocusing_param.focal_depth: Focal depth (nm) used in the fitting.
% defocusing_param.center_w: Width (nm) at the central focal plane.
% defocusing_param.central_z_readouts: The original piezo stage z-readouts
% (nm) for the central plane.
% defocusing_param.z_range: The z-sampling range stored in a 3-element
% vector [minimal z(nm), maximal z(nm), sampling interval (nm)];

%-------------------------------------------------------------------------%

%% ------------------------------------------------------------------------
%  Extract calibration data and fit it with a defocusing equation
%  ------------------------------------------------------------------------
if size(cali_data, 1) ~= 5
    display('The input calibration data does not have the correct format.')
    return
end

z=cali_data(1,:);
sigma_x=cali_data(2,:);
sigma_y=cali_data(3,:);

% Plot the average X- and Y-widths versus z-readout and ask the user to
% left- and right-click to specify the z-distance range over which the
% central plane location will be found and the fitting will be done. I
% avoid making the code to automatically find the z-location associated
% with the minimal aspect ratio since the minimal aspect ratio can also
% occur if the emitter is too out of focus. Thus I ask the user to specify
% a range here.

figure
% errorbar(z, sigma_x, stderr_x_width, 'b.');
% I commented the above line since the standard error is very small - no
% need to display it.

plot(z,sigma_x,'bo--','MarkerSize',5,'MarkerFaceColor','b'); hold all,
plot(z,sigma_y,'ro--','MarkerSize',5,'MarkerFaceColor','r');
set(gca,'xlim',[z(1),z(end)],'ylim',...
    [min(min(sigma_x),min(sigma_y)),max(max(sigma_x),max(sigma_y))]);
ylabel('Width (nm)','FontSize',16)
xlabel('Piezo Stage Z-Readouts (nm)','FontSize',16)
legend('Wx','Wy'); set(gca,'FontSize', 6)
grid on

title('Left- and Right-Click to Specify the Fitting Domain', 'FontSize', 16);
[z_readout_range,~,~]=jsbginput(2);

% Remove values that fall outside the chosen domain
sigma_x=sigma_x(z>z_readout_range(1)&z<z_readout_range(2));
sigma_y=sigma_y(z>z_readout_range(1)&z<z_readout_range(2));
z=z(z>z_readout_range(1)&z<z_readout_range(2));

% Determine the z-location that correspond to the central focusing
% plane, that is, the z-location where x-width and y-width are essentially
% equal.
rato=sigma_x./sigma_y;
[~,equal_width_loc_ind]=min(abs(rato-1));
subrange=equal_width_loc_ind-2:equal_width_loc_ind+2;
% central_plane_z=z(equal_width_loc_ind);
central_plane_z=interp1(rato(subrange),z(subrange),1);
% central_plane_wx=interp1(rato(subrange),sigma_x(subrange),1);
% central_plane_wy=interp1(rato(subrange),sigma_y(subrange),1);

% Subtract 'central_plane_z' from the original z-readout values to
% get the z-distance from the central plane, which will be used for
% defocusing function fitting
fit_z=z-central_plane_z;

if ~ischar(fit_z_range)&&length(fit_z_range)==2
    display('Defocusing function fitting will use specified z-readouts values only.')
    fit_z(fit_z<fit_z_range(1)|fit_z>fit_z_range(2))=nan;
else
    display('Defocusing functiong fitting will use all available z-readouts among the chosen range.')
end

options=statset('MaxFunEvals',50000,'MaxIter',50000);
astigmatism_defocusing_eq=@(p,z,w)...
    w*sqrt(1+((z-p(1))/p(2)).^2+p(3)*((z-p(1))/p(2)).^3+p(4)*((z-p(1))/p(2)).^4);
lb=[-1e3,0,0,0];
ub=[1e3,inf,inf,inf];

w0=[sigma_x(sigma_x==min(sigma_x)),sigma_y(sigma_y==min(sigma_y))];

% Fit x-width to defocusing equation
[x_fitparam,~,~]=lsqcurvefit(@(param,z)astigmatism_defocusing_eq(param,z,w0(1)),...
    [guess_param(1),guess_param(3:end)],fit_z,sigma_x,lb,ub,options);

% Fit y-width to defocusing equation
[y_fitparam,~,~]=lsqcurvefit(@(param,z)astigmatism_defocusing_eq(param,z,w0(2)),...
    [guess_param(2),guess_param(3:end)],fit_z,sigma_y,[],[],options);

% Define output variable
defocusing_param.fitparam=[x_fitparam;y_fitparam];
% defocusing_param.w0=w0;
defocusing_param.center_w=w0;
defocusing_param.central_z_readouts=central_plane_z;
defocusing_param.z_range=[min(fit_z),max(fit_z),abs(fit_z(2)-fit_z(1))];
%% ------------------------------------------------------------------------
%  Plot the Fitting Curve
%  ------------------------------------------------------------------------

% Plot the defocusing data
figure
plot(fit_z,sigma_x,'bo','MarkerFaceColor','b','MarkerEdgeColor','b'),
hold all,
plot(fit_z,sigma_y,'ro','MarkerFaceColor','r','MarkerEdgeColor','r'),

% Plot the defocusing fit curve
hold all
sample_z=linspace(-1e3,1e3); %fit_z(1):5:fit_z(end); % Sample the width every 5 nm of z value

plot(sample_z,astigmatism_defocusing_eq(x_fitparam,sample_z,w0(1)),...
    'b--','LineWidth',2); % x-defocusing curve
hold all
plot(sample_z, astigmatism_defocusing_eq(y_fitparam,sample_z,w0(2)), ...
    'r--','LineWidth',2); % x-defocusing curve
axis tight
% % Lower and upper bounds for plotting.
% z_plot_lb=fit_z(find(~isnan(fit_z),1,'first'));
% z_plot_ub=fit_z(find(~isnan(fit_z),1,'last'));
% 
% set(gca,'xlim',[z_plot_lb,z_plot_ub],'ylim',[0,max(max(sigma_x),max(sigma_y))]);
ylabel('Width (nm)','FontSize',16)
xlabel('Z-Distance from Central Plane (nm)','FontSize',16)
set(gca,'FontSize',16)
grid on

% Save the figure under the same directory where the calibration data are
saveas(gcf,'Calibration Curve.fig');
end % End of the function 'fit_defocusing_eq'.