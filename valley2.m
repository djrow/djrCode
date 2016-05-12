function PhaseMask=valley2(f,params)
% clear all, close all;clc
% Update History:

% 6/27/2013 YL: Added the "autofill" function. Now it is OK to use lower
% values for "higher_thresh" and/or higher values for "lower_thresh"
% without worrying too much about over-segmentation (i.e., cell masks
% having holes inside). Also, the connectivity criterion for the binary
% labelling has been changed from 8(default) to 4 in the autofill section.

% 6/172013 YL: Updated Comments.

% 5/3/2016 DR: parameters are now fields of a structure

%% ------------------------------------------------------------------------
%  User-Defined Parameters:
%  ------------------------------------------------------------------------

dilate_factor=params.nDilation; % 0;
% Change the to zero only if you still have under-segmentation regardless
% of threshold.

lower_thresh=params.lowThresh;% 0;
% Can be negative, but usually use [-1 0.5].

higher_thresh=params.highThresh;% .001;
% Usually use positive values [0 3] , should be larger than
% lower_threshold. *Decrease* this value first when cells are
% UNDER-segmented (e.g., 2 or more cells grouped as 1 region).

autofill=params.autofill;% 1;

% Note: If cells are under-segmented, increase the lower_thresh and/or
% decrease the higher_thresh (recommended); If cells are over-segmented,
% decrease the lower_thresh and increase the higher_thresh.
min_area=params.minArea;% 500;
% minimum allowable area (in px) for a single segmented region. Regions
% with area smaller than this will be discarded.
max_area=params.maxArea; % 10000;
% maximum allowable area (in px) for a single segmented region. Regions
% with area larger than this will be discarded.
%%

% [wlimg_name, wlimg_path, ~] = uigetfile({'*.tif','*.tiff'}, ...
%     'Select a white-light image');
% if wlimg_path == 0
%     display('No white light image selected. Aborting program.')
%     return
% end
% 
% frame=imread([wlimg_path wlimg_name]);
% f = frame;
% figure, imshow(f, []); title('orginal')
f=uint16(f);
f=imcomplement(f); % Invert intensity

% Laplacian of Gaussian filtering
[g,~]=edge(f,'log', 0);
% figure,imshow(~g,[])

% h=fspecial('sobel');
% fd=double(f);
% g=sqrt(imfilter(fd,h,'replicate').^2+imfilter(fd,h','replicate').^2);
% figure, imshow(g,[]); title('sobel');

f = imfilter(f, fspecial('average', 5), 'replicate');
% figure, imshow(f, []); title('average- or gaussian- filtered image');


% figure, imshow(f, []); title('inverted')

%-------------------------------------------------------------------------%
% Define Valley filters
V = zeros(3,3); V(2,2) = -1;
A1 = V; A1(3,1) = 1;
A2 = V; A2(2,1) = 1;
A3 = V; A3(1,1) = 1;
B1 = V; B1(3,2) = 1;
B3 = V; B3(1,2) = 1;
C1 = V; C1(3,3) = 1;
C2 = V; C2(2,3) = 1;
C3 = V; C3(1,3) = 1;
%-------------------------------------------------------------------------%

A1 = imfilter(f, A1, 'corr', 'replicate', 'same');
C3 = imfilter(f, C3, 'corr', 'replicate', 'same');
A2 = imfilter(f, A2, 'corr', 'replicate', 'same');
C2 = imfilter(f, C2, 'corr', 'replicate', 'same');
A3 = imfilter(f, A3, 'corr', 'replicate', 'same');
C1 = imfilter(f, C1, 'corr', 'replicate', 'same');
B3 = imfilter(f, B3, 'corr', 'replicate', 'same');
B1 = imfilter(f, B1, 'corr', 'replicate', 'same');

V1 = min(A1, C3); V2 = min(A2, C2); V3 = min(A3, C1); V4 = min(B3, B1);
V1 = max(V1, V2); V2 = max(V3, V4); V = max(V1, V2);

% figure, imshow(V,[]), title('valley values')
% figure, imshow(V~=0,[]), title('non-zero valley values')

% 8-connectivity correlation
mean_V = mean(V(:)); std_V = std2(V);
T = [mean_V+lower_thresh*std_V, mean_V + higher_thresh*std_V];
V2 = uint8(V>= T(2))*10; V2(V2~=10)=1; % Strong threshold
V1 = uint8(V>= T(1)); % Weak threshold
V1 = imfilter(V1.*V2, [1 1 1; 1 1 1; 1 1 1], 'corr', 'replicate', 'same')>=10;

% figure, imshow(V2,[]); title('strong threshold')
% figure, imshow(~V1,[]); title('strong and weak threshold')

% Thresholding
f2 = im2bw(f, graythresh(f));
% figure, imshow(f2, []); title('thresholding')

f3 = f2 & ~V1 & ~g;
% figure, imshow(f3, []); title('threshold and edge combined')

f4 = bwmorph(f3,'close',inf);
% figure, imshow(f4,[]); title('close')

f5 = imfill(f4,'holes');
% figure, imshow(f5,[]); title('imfill')

f6 = bwmorph(f5, 'majority', inf');
% figure,imshow(f6,[]); title('majority')

f7 = bwlabel(f6, 4);
% figure, imshow(f7, []); title('bwlabel')

seg_area = regionprops(f7, 'area');

% for i = 1:length(seg_area)
%     if seg_area(i,1).Area < min_area || seg_area(i,1).Area > max_area
%         f7(f7 == i) = 0; % Get rid of small regions
%     end
% end
f7(ismember(f7,find([seg_area.Area]<min_area|[seg_area.Area]>max_area)))=0;

f7 = bwmorph(f7, 'dilate', dilate_factor);
% figure, imshow(f7, []); title('final dilation by 1')
% f7 = bwlabel(f7 ~= 0, 4);
% figure, imshow(f7, []); title('bwlabel without small regions and dilated by 1')
% % Reassign region label after getting rid of small regions

f7 = bwlabel(~bwmorph(~f7,'diag',5));
% figure, imshow(f7, []); title('Unfilled Phase Mask')

conv_hull = regionprops(f7, 'Convexhull');

% for i = 1:length(conv_hull)
%     hold all,
%     plot(conv_hull(i,1).ConvexHull(:,1),conv_hull(i,1).ConvexHull(:,2), 'r-' ,'linewidth', 2)
%     axis ij equal
% end

% dah= bwconvhull(f7~=0,'objects',4);figure,imshow(dah)

PhaseMask = f7;

if autofill == 1
    PhaseMask = bwconvhull(PhaseMask~=0,'objects',4);
    PhaseMask = bwlabel(PhaseMask, 4);
    conv_hull = regionprops(PhaseMask, 'Convexhull');
%     close all
%     figure, imshow(frame, []); title('orginal');
%     figure,imshow(PhaseMask, []); title('Autofilled Phase Mask');
%     for i = 1:length(conv_hull)
%         hold all,
%         plot(conv_hull(i,1).ConvexHull(:,1),conv_hull(i,1).ConvexHull(:,2), 'c-' ,'linewidth', 2)
%         axis ij equal
%     end
end
end

% %-------------------------------------------------------------------------%
% % Plot contour of convex hull and scale it up by x10
% conv_contour = zeros(size(f7).*10);
% figure,imshow(conv_contour,'border','tight');
% for i = 1:length(conv_hull)
%     hold all,
%     plot(conv_hull(i,1).ConvexHull(:,1).*10,conv_hull(i,1).ConvexHull(:,2).*10, 'c-' ,'linewidth', 2)
%     axis ij equal
% end
% %-------------------------------------------------------------------------%