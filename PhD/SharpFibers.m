function SharpFibers
Sites={'F042'};%{'F064','F065','F015','F039','F056','F074'};
cells={[25,30,29,34,26],[7 9 11 6 3],[11 15 20 26 21 27 32 40 34 28 25 17],[3 6 8 10 15 16 12 9 13 20 19 21],[23 18 14 13 5 8 6 27 16],[11,5,10,4,16]};
files=dir('.');
files=extractfield(files,'name');
numberOfunits=[];
k=1;
for site=1:length(Sites)
    idx = all(ismember(files,Sites{site}),2);
    
basename=find(~cellfun(@isempty,regexp(files,sprintf('%s\\w*C03.jpg',Sites{site}))));
basename2=find(~cellfun(@isempty,regexp(files,sprintf('%s\\w*C01\\w*Cells',Sites{site}))));
basename3=find(~cellfun(@isempty,regexp(files,sprintf('%s\\w*C01\\w*Nuclei',Sites{site}))));
 
I = double(imread(files{basename}))/2^16;
%mask=imread((files{basename2}));
%Nuclei=imread((files{basename3}));
 
if ~exist(sprintf('angles%s.mat',Sites{site}))    
I=repelem(I,3,3);    
I=I-0.0007;
I(I<0)=0;
L=[0 -1 0;-1 4 -1;0 -1 0];
%img=zeros(size(I,1)+4,size(I,2)+4,60);
i=1;
img=[];
%as suggested in https://link.springer.com/article/10.1007%2Fs10237-015-0706-9#Sec2
c=linspace(0,180,30);
for theta=c   
    gauss=customgauss([30 30],theta,sqrt(2));
    LOG=conv2(L,gauss);
    imgC=conv2(LOG,I);
    img(:,:,i)=imgC;
    i=i+1;
end
[maxI,indI]=max(img,[],3);
 
for i=1:30
    indI(indI==i)=c(i);
end
maxI(maxI<0.001)=0;
level=graythresh(maxI);
maxI=logical(maxI);
% figure;
%  imagesc(indI)
%  colormap jet
%  axis off
% axis equal
%maxI=maxI(1080:1157,2256:2353);
%indI=indI(1080:1157,2256:2353);
EnchImg=CoherenceFilter(maxI); %coherence-enhancing diffusion filtering (CEDF)
level = graythresh(EnchImg);
BW = imbinarize(EnchImg,level);
[finalFib,indI]=localOrientFilt(BW,indI);
finalFib=bwareaopen(finalFib, 50);
indI=indI.*finalFib;
% figure;
% imagesc(indI);
% colormap jet
% axis off
% forHist=reshape(indI,[],1);
%     forHist(forHist==0)=[];
%     binranges=min(forHist):6:max(forHist);
%     bincounts=histc(forHist,binranges)/length(forHist);
%     bar( binranges,bincounts);
save(sprintf('angles%s.mat',Sites{site}),'finalFib','indI')
rgbImage = ind2rgb(uint16(indI), jet(180));
imwrite(rgbImage,sprintf('%s_map.png',Sites{site}),'BitDepth',16)
else
    load(sprintf('angles%s.mat',Sites{site}))
end
%figure;
FA=[];
i=1;
for cell=16%cells{site}
    maskCell=repelem(mask,3,3);
    maskNuc=Nuclei;
    maskCell(maskCell~=cell)=0;
%     im=imfuse(I,maskCell,'ColorChannels',[1 2 0]);
%     BW = roipoly(im);
%     maskCell=BW;
    maskNuc(maskNuc~=cell)=0;
    maskNuc=logical(maskNuc);
    maskCell=logical(maskCell);
    maskNuc=repelem(maskNuc,3,3);
    propsMN=regionprops(maskNuc,'ConvexArea','Area');
    Solidity=propsMN.Area/propsMN.ConvexArea;    
    maskCell=logical(maskCell);
    cellIm=indI(16:end-16,16:end-16).*maskCell;
    cellMask=indI(16:end-16,16:end-16).*maskCell;
    
    cent=regionprops(maskNuc,'Centroid');
    sectors={};
    r=size( maskCell,2);
    i=1;
    currsector=1;
    unitsn=1;
if ~exist(sprintf('%s_sectors_cell%d_sectorsRaw.mat',Sites{site},cell))
    for theta=2*pi/100:2*pi/100:2*pi
     alpha=theta-2*pi/100:0.0001:theta;
     xr = r * cos(alpha) + round(cent.Centroid(1));
     yr = r * sin(alpha) + round(cent.Centroid(2));
     x = [round(cent.Centroid(1)), xr, cent.Centroid(1)];
     y = [round(cent.Centroid(2)), yr, round(cent.Centroid(2))];
     sector = poly2mask(x,y, size( maskCell, 1), size( maskCell, 2));
     sectorIm=cellIm.*sector;
     sectorMask=cellMask.*sector;
     meanAng=sum(sectorIm(:))/length(find(sectorMask));
     if  theta~=2*pi/100
         if  abs(meanAngOld-meanAng)<12%~kstest2(sectorImOld(find(sectorImOld)),sectorIm(find(sectorIm)),'Alpha',0.01)%
             xrnew=[xrOld xr];
             yrnew=[yrOld yr];
             x = [round(cent.Centroid(1)), xr, round(cent.Centroid(1))];
             y = [round(cent.Centroid(2)), yr, round(cent.Centroid(2))];
         sector = poly2mask(x,y, size(cellIm, 1), size(cellIm, 2));
        sectorImN=cellIm.*sector;
        sectorMaskN=cellMask.*sector;
        meanAngnew=sum(sectorImN(:))/length(find(sectorMaskN));
        if std(sectorImN(find(sectorMaskN)))<36
            xr=xrnew;
            yr=yrnew;
            meanAng=meanAngnew;
            sectorIm=sectorImN;
            unitsn(i)=unitsn(i)+1;
            else
            currsector=currsector+1; 
            i=i+1;
        end
         else
            currsector=currsector+1; 
            i=i+1;
         end
     end
     %imshow(sectorIm)     
     sectors{currsector}=[xr;yr];
     meanAngOld=meanAng;
     xrOld=xr;
     yrOld=yr;
     if length(unitsn)<i
        unitsn(i)=1;
     end
    end
numberOfunits{k}=[Solidity unitsn];
save(sprintf('%s_sectors_cell%d_sectorsRaw.mat',Sites{site},cell),'sectors','numberOfunits','unitsn')
else
    load(sprintf('%s_sectors_cell%d_sectorsRaw.mat',Sites{site},cell))
end
 
if ~exist(sprintf('%s_Dirsectors_cell%d_sectors*',Sites{site},cell))
sector=zeros(size(maskCell,1),size(maskCell,2),length(sectors));
theta_start=0;
numaligned=0;
alignedmatrix=zeros(1,length(sectors));
cellIm(1:round(cent.Centroid(2)),:)=cellIm(1:round(cent.Centroid(2)),:)+180;
for i=1:length(sectors)
    theta_stop=theta_start+unitsn(i)*2*pi/100;
    mtheta=180*(theta_stop+theta_start)/2/pi;    
    x = [cent.Centroid(1), sectors{i}(1,:), cent.Centroid(1)];
     y = [cent.Centroid(2), sectors{i}(2,:), cent.Centroid(2)];
     sector(:,:,i) = double(poly2mask(x,y, size(maskCell, 1), size(maskCell, 2))); 
     cellsector=cellIm.* sector(:,:,i);
     cellsectormask=cellMask.* sector(:,:,i);
     
     cellsectorm=cellsector(find(cellsectormask));
     idx=length(cellsectorm(cellsectorm>180*theta_start/pi & cellsectorm<180*theta_stop/pi))/length(cellsectorm);
     alignedmatrix(i)=idx;
     sectororig(:,:,i)=sector(:,:,i)*i;
     sector(:,:,i)=sector(:,:,i)*idx;
     
     theta_start=theta_stop;
end
numberOfunits{k}=[numberOfunits{k};Solidity alignedmatrix];
k=k+1;
c=distinguishable_colors(length(sectors)+1);
c(1,1:3)=1;
allsec=max(sector,[],3);
sectororig2=max(sectororig,[],3);
allsec=allsec.*maskCell;
sectororig2=sectororig2.*maskCell;
rgbImage2=ind2rgb(uint8(sectororig2), c);
rgbImage = uint8(allsec*255);%ind2rgb(uint8(allsec), c);
imwrite(rgbImage,sprintf('%s_Dirsectors_cell%d_sectors%d.png',Sites{site},cell,numaligned),'BitDepth',8)
imwrite(rgbImage2,sprintf('%s_sectors_cell%d_sectors%d.png',Sites{site},cell,length(sectors)),'BitDepth',8)
end
end
end
%save('Sol_numunits.mat','numberOfunits');
%axis equal
end
 
function [EnchImg,indI]=localOrientFilt(EnchImg,indI)
%angles=linspace(0,360,60);
newImg=zeros(size(EnchImg,1),size(EnchImg,2));
 
while sum(sum(abs(EnchImg-newImg)))>0
    fprintf('%d\n',sum(sum((EnchImg-newImg))))
    newImg=EnchImg;
    indI=EnchImg.*indI;
    bs=5;
for i= ceil(bs/2): size(EnchImg,1)- ceil(bs/2)
    for j= ceil(bs/2): size(EnchImg,2)- ceil(bs/2)
         if EnchImg(i,j)~=0
             box2=indI(i-floor(bs/2):i+floor(bs/2),j-floor(bs/2):j+floor(bs/2));
            tt=cos(pi*(box2(ceil(bs/2),ceil(bs/2))-(sum(box2(:))-box2(ceil(bs/2),ceil(bs/2)))/(length(find(box2))-1))/180);
            if abs(tt)<0.97 || length(find(box2))<5
                EnchImg(i,j)=0;
            end
        end
    end
end
fprintf('%d\n',sum(sum((EnchImg-newImg))))
end
indI=EnchImg.*indI;
end
 
function ret = customgauss(gsize, theta,sigma)
sigmax=sigma;
sigmay=sigma*10;
offset=0;
factor=1; 
center=[0 0];
ret     = zeros(gsize);
rbegin  = -round(gsize(1) / 2);
cbegin  = -round(gsize(2) / 2);
for r=1:gsize(1)
    for c=1:gsize(2)
        ret(r,c) = rotgauss(rbegin+r,cbegin+c, theta, sigmax, sigmay, offset, factor, center);
    end
end
%imagesc(ret)
end
function val = rotgauss(x, y, theta, sigmax, sigmay, offset, factor, center)
xc      = center(1);
yc      = center(2);
theta   = (theta/180)*pi;
xm      = (x-xc)*cos(theta) - (y-yc)*sin(theta);
ym      = (x-xc)*sin(theta) + (y-yc)*cos(theta);
u       = (xm/sigmax)^2 + (ym/sigmay)^2;
val     = offset + factor*exp(-u/2);
end
 
function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);
% Copyright 2010-2011 by Timothy E. Holy
  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
    bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end
  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end
function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end
function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1
    c = rgbspec(k,:);
  elseif length(c)>2
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end
