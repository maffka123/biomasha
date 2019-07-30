function skeletANDcross
list=ls('images/*_F04_*C03.png*');
for kk=1%:size(list,1)
%figure;
img=list(kk,:);
basenameN=strfind(img,'01C03')+4;
basename=img(1:basenameN);
name=strcat('images/',basename,'_eq_Probabilities.tiff');
P=im2double(imread(name));
img=P(:,:,1);
img(img<0.6)=0;img(img>=0.6)=1;
 
imgClean=bwareaopen(img,20);
 
skelimg = bwmorph(imgClean,'skel', Inf);
skel = bwmorph(skelimg,'thin');
cross = bwmorph(skel, 'branchpoints');
for i=1:6
ep= bwmorph(skel, 'endpoints');
D = bwdistgeodesic(skel,find(cross),'quasi'); D(isnan(D))=0;
epD=ep.*D;
epD(epD<6 & epD>0)=NaN; epD(~isnan(epD))=1; epD(isnan(epD))=0;
skel=logical(skel.*epD);
cross = bwmorph(skel, 'branchpoints');
end
skel=imclose(skel,1);
skel=bwmorph(skel,'thin');
skel=imclose(skel,[1 1;1 1]);
skel=bwmorph(skel,'thin');
skel=imclose(skel,[1 1 1;1 1 1]);
skel=bwmorph(skel,'thin');
skel=imclose(skel,[1 1 1 1;1 1 1 1]);
skel=bwmorph(skel,'thin');
skel=bwareaopen(skel,5);
cross = bwmorph(skel, 'branchpoints');
cross=filterCross(cross);
 fin=imoverlay(img,skel,[0 1 0]);
% 
 fin2=imoverlay(fin,cross,[0.9 0.2 0.2]);
 
 
p3=P(:,:,2);p3(p3<0.6)=0;p3(p3>=0.6)=1;
 fin3=imoverlay(fin2,p3,[0.9 0.9 0.2]);
 
 dapiname=strcat('SEGMENTATION\',strrep(basename,'A03Z01C03','A01Z01C01'),'_SegmentedKDNuclei.png');
 dapi=logical(im2double(imread(dapiname)));
% 
Segmentation=im2double(imread(strcat('SEGMENTATION\',strrep(basename,'A03Z01C03','A01Z01C01'),'_SegmentedKDCytoplasm.png')));
B= bwboundaries(Segmentation);
% 
Cell=im2double(imread(strcat('SEGMENTATION\',strrep(basename,'A03Z01C03','A01Z01C01'),'_SegmentedCells.png')));
Cell(Cell~=13)=0;
 
 fin4=imoverlay(fin3,dapi,[0.2 0.2 0.9]);
 
cross=cross-dapi;
cross(cross<0)=0;
cross=cross-p3;
cross(cross<0)=0;
cross=cross-Cell;
cross(cross<0)=0;
 
figure; imshow(fin4)
 
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'color',[0.87 0.49 0], 'LineWidth', 2)
end
%rgbImage = ind2rgb(uint16(fin4),colormap);
imwrite(fin4,sprintf('Segm/%s_map.png',basename),'BitDepth',8)
cross(:,1:10)=0; cross(1:10,:)=0;cross(:,end-10:end)=0;cross(end-10:end,:)=0;
[crosidx1,crosidx2]=find(cross);
numNot3wJ=0; num3wJ=0;
for i=1:sum(cross(:))
%    flag=0;
    padd=skel(crosidx1(i)-4:crosidx1(i)+4,crosidx2(i)-4:crosidx2(i)+4);
%     image(padd, 'XData', [1 1], 'YData', [1 1],'CDataMapping','scaled');hold on
    pad2=zeros(9,9); pad2(4:6,4:6)=1; pad2=pad2.*padd;
    gg=[32 31 40 49 50 51 42 33];
    ggpos=find(padd(gg)==1);
    [p1,p2]=ind2sub([9,9],gg(ggpos));
    p1=p1';p2=p2';%c=(p1==5 & p2==5); p1(c)=[];p2(c)=[];
    if length(p1)~=3
        numNot3wJ=numNot3wJ+1;
        continue;
    end
    num3wJ=num3wJ+1;
    points=cell(length(p1),2);
    Fit=cell(length(p1),1);
    points(:)={5};
    for k=1:3
        points{k,1}(2)=p1(k); points{k,2}(2)=p2(k);
        for j=2:4
            pad2=zeros(9,9); pad2(points{k,1}(j)-1:points{k,1}(j)+1,points{k,2}(j)-1:points{k,2}(j)+1)=1; pad2=pad2.*padd;
        [pp1,pp2]=find(pad2);c=(pp1==points{k,1}(j) & pp2==points{k,2}(j)); pp1(c)=[];pp2(c)=[];
        c=(pp1==points{k,1}(j-1) & pp2==points{k,2}(j-1)); pp1(c)=[];pp2(c)=[];
         try
       [~,c]=intersect([pp1 pp2],[p1 p2],'rows'); pp1(c)=[];pp2(c)=[];
         catch
            break
         end
       if length(pp1)~=1
%            num3wJ=num3wJ-1;
%            numNot3wJ=numNot3wJ+1;
%           flag=1;
           break            
       end
        points{k,1}(j+1)=pp1(1);points{k,2}(j+1)=pp2(1);
        end
%         if flag==1
%             flag=0;
%             break
%         end
warning off
        p=polyfit(points{k,2},10-points{k,1},1);
        x=[points{k,2}(1) points{k,2}(end)];
        
        if polyval(p,x(1))==polyval(p,x(2))
           Fit{k}=[x' [10-points{k,1}(1) 10-points{k,1}(end)]']; 
        else
        Fit{k}=[x' polyval(p,x)'];
        end
%         hold on
%         plot(Fit{k}(:,1),10-Fit{k}(:,2),'LineWidth',5)
    end
  
    cc=[1 2; 2 3; 3 1];
    Fit2=cellfun(@diff,Fit,'UniformOutput' , false);
 warning on 
    angle=zeros(3,1);
    for j=1:size(cc,1)
        angle(j)=atan2d(det([Fit2{cc(j,1)}; Fit2{cc(j,2)}]),dot(Fit2{cc(j,1)},Fit2{cc(j,2)}));%(Fit(cc(j,1))-Fit(cc(j,2)))*180/pi;
        if angle(j)>0
            angle(j)=360-angle(j);
        end
    end
    Ag1(:,i)=angle;
    %hold off
end
Ag{kk}=Ag1;
end
ep= bwmorph(skel, 'endpoints');
D = bwdistgeodesic(skel,find(cross),'quasi');
dd=D(ep); dd(dd==Inf)=[];
mDist=mean(dd);
Segmentation=Segmentation-p3;
Segmentation(Segmentation<0)=0;
Segmentation(Segmentation~=0)=1;
img=imgClean.*Segmentation;
num3wJ=num3wJ/sum(img(:));
ep=sum(ep(:))/sum(img(:));
 [i,j] = find(bwmorph(skel,'endpoints'));
D = bwdistgeodesic(skel,find(cross),'quasi');
 
for n = 1:numel(i)
    text(j(n),i(n),[num2str(D(i(n),j(n)))],'color','g');
end
Ag=cellfun(@abs,Ag,'UniformOutput' , false);
 
t=cell2mat(Ag);
ff=find(t(1,:)==0);
t(:,ff)=[];
three_way_angles(t(1,:),t(2,:),t(3,:))
 
t2=sort(t);
binranges=0:5:360;
bincounts1=histc(t2(1,:),binranges)/size(t2,2);
bincounts2=histc(t2(2,:),binranges)/size(t2,2);
bincounts3=histc(t2(3,:),binranges)/size(t2,2);
subplot(1,3,1)
bar( binranges,bincounts1);
xlim([0 200]); ylim([0 0.3])
subplot(1,3,2)
bar( binranges,bincounts2);
xlim([0 200]); ylim([0 0.3])
subplot(1,3,3)
bar( binranges,bincounts3);
xlim([30 200]); ylim([0 0.3])
end
 
function crossO=filterCross(crossI)
filt(:,:,1)=[-1 0 0;0 1 0;0 0 -1];
filt(:,:,2)=[0 0 -1;0 1 0;-1 0 0];
filt(:,:,3)=[0 -1 0;0 1 0;0 -1 0];
filt(:,:,4)=[0 0 0;-1 1 -1;0 0 0];
cross=crossI;
for i=1:4
   cross= bwhitmiss(cross,filt(:,:,i));
end
crossO=cross;
end
function BW_Thinned=ZhangSuen(BW_Original)
changing = 1;
[rows, columns] = size(BW_Original);
BW_Thinned = BW_Original;
BW_Del = ones(rows, columns);
while changing
    % BW_Del = ones(rows, columns);
    changing = 0;
    % Setp 1
    for i=2:rows-1
        for j = 2:columns-1
            P = [BW_Thinned(i,j) BW_Thinned(i-1,j) BW_Thinned(i-1,j+1) BW_Thinned(i,j+1) BW_Thinned(i+1,j+1) BW_Thinned(i+1,j) BW_Thinned(i+1,j-1) BW_Thinned(i,j-1) BW_Thinned(i-1,j-1) BW_Thinned(i-1,j)]; % P1, P2, P3, ... , P8, P9, P2
            if (BW_Thinned(i,j) == 1 &&  sum(P(2:end-1))<=6 && sum(P(2:end-1)) >=2 && P(2)*P(4)*P(6)==0 && P(4)*P(6)*P(8)==0)   % conditions
                % No. of 0,1 patterns (transitions from 0 to 1) in the ordered sequence
                A = 0;
                for k = 2:size(P,2)-1
                    if P(k) == 0 && P(k+1)==1
                        A = A+1;
                    end
                end
                if (A==1)
                    BW_Del(i,j)=0;
                    changing = 1;
                end
            end
        end
    end
    BW_Thinned = BW_Thinned.*BW_Del;  % the deletion must after all the pixels have been visited
    % Step 2 
    for i=2:rows-1
        for j = 2:columns-1
            P = [BW_Thinned(i,j) BW_Thinned(i-1,j) BW_Thinned(i-1,j+1) BW_Thinned(i,j+1) BW_Thinned(i+1,j+1) BW_Thinned(i+1,j) BW_Thinned(i+1,j-1) BW_Thinned(i,j-1) BW_Thinned(i-1,j-1) BW_Thinned(i-1,j)];
            if (BW_Thinned(i,j) == 1 && sum(P(2:end-1))<=6 && sum(P(2:end-1)) >=2 && P(2)*P(4)*P(8)==0 && P(2)*P(6)*P(8)==0)   % conditions
                A = 0;
                for k = 2:size(P,2)-1
                    if P(k) == 0 && P(k+1)==1
                        A = A+1;
                    end
                end
                if (A==1)
                    BW_Del(i,j)=0;
                    changing = 1;
                end
            end
        end
    end
    BW_Thinned = BW_Thinned.*BW_Del;
end%while
end
