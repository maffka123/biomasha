function probModel(probfusion, probcollapse, probpickup, pool, saveimg)
warning('off','all')
if nargin==0
    probfusion=70;
    probcollapse=10;
    probpickup=60;
    saveimg=0;
    pool=1000; %limited number of subunits
end
ss=200;
w3j=zeros(ss,ss);
bkg=zeros(ss,ss);
bkg(ss/2,ss/2)=1;
bkgER=generatebkgER(ss,10);
mts={};
mtsER=[];
conn2mt=[];
AcMT=[];
ERs={};
ERpool=5000;
for t=1:10000
    fprintf('time: %d \n',t)
    if pool>0
        newx=ss/2+randi([2,4],1)-3;
        newy=ss/2+randi([2,4],1)-3;  
        if ~(newx==ss/2 && newy==ss/2)
            pool=pool-1;
            mts{1+length(mts)}=[ss/2 ss/2 ;newx,newy];
            if randi([0,100],1)<probpickup && ERpool>0
                mtsER(length(mts))=1;
                ERs{1+length(ERs)}=[newx,newy];
                bkgER(newx,newy)=1;
                conn2mt(length(ERs))=length(mts);
                ERpool=ERpool-1;
            else
                mtsER(length(mts))=0;   
            end
            AcMT(length(mts))=0;
        end
    end
for mt=1:length(mts)
    if length(mts{mt})>30 && mtsER(mt)==1 && ~isempty(find(conn2mt==mt)) && randi([0,100],1)<1 % long MT attached to the ER can get acetylated
        AcMT(mt)=1;
    end
    canfuse=0;
    if ((randi([0,100],1)<1 || (randi([0,100],1)<probcollapse && mtsER(mt)==1)) && AcMT(mt)==0) ...% 1% for MT to collapse or collapse if attached to ER
       ||  (randi([0,500],1)<1 && AcMT(mt)==1)
        pool=pool+length(mts{mt})-1;
         mts{mt}=[];
        if mtsER(mt)==1 && ~isempty(find(conn2mt==mt))
                if ~isempty(intersect(sub2ind(size(bkgER),ERs{find(conn2mt==mt)}(:,1),ERs{find(conn2mt==mt)}(:,2)),find(w3j),'stable'))
                    line=ERs{find(conn2mt==mt)};
                    points=intersect(sub2ind(size(bkgER),line(:,1),line(:,2)),find(w3j));
                    [p1,p2]=ind2sub(size(bkgER),points);
                    B=[p1,p2];
                    distances = sqrt(sum(bsxfun(@minus, line(end,end), B).^2,2));
                    closest = B(find(distances==min(distances)),:);
                    pos=find(ismember(line,closest,'rows'));
                    bkgER(sub2ind(size(bkgER),ERs{find(conn2mt==mt)}(pos:end,1),ERs{find(conn2mt==mt)}(pos:end,2)))=0;
                    ERpool=ERpool+length(ERs{find(conn2mt==mt)}(pos:end,:));
                    ERs{find(conn2mt==mt)}(pos:end,:)=[]; 
                else
            ERpool=ERpool+length(ERs{find(conn2mt==mt)});
            bkgER(sub2ind(size(bkgER),ERs{find(conn2mt==mt)}(:,1),ERs{find(conn2mt==mt)}(:,2)))=0;
            ERs{find(conn2mt==mt)}=[];
                end
                mtsER(mt)=0; conn2mt(conn2mt==mt)=0; 
        end     
%        cleanANDdraw(mtsER,mts,conn2mt,ERs,ss,1)
    elseif mts{mt}(end,1)+1<ss && mts{mt}(end,1)-1>0 && mts{mt}(end,2)+1<ss && mts{mt}(end,2)-1>0 && pool>0 && length(mts{mt})+1<80  % if not collapsed         
        newx=randi([2,4],1)-3; newy=randi([2,4],1)-3;
        if newx==0 && newy==0; continue;  end;
        border=DetectBorders(bkgER); % can only pool ER from the edges
            pos2=cellfun(@(x) ismember(mts{mt},x,'rows'),ERs,'UniformOutput',false,'ErrorHandler',@errorfun);
            pos3=cellfun(@(x) length(find(x(2:end)~=0)),pos2);
        if ((randi([0,100],1)<probpickup && border(mts{mt}(end,1),mts{mt}(end,2))~=0 && w3j(mts{mt}(end,1),mts{mt}(end,2))~=1 && isempty(find(pos3>5)))... %Probability to pick up ER tubule but not from 3wj, not from long ER tubule
                || mtsER(mt)==1) && ERpool>0 % or existing ER-MT connection
            mtsER(mt)=1; 
            if size(mts{mt},1)>1
                [~,newx,newy]=checkconn(mts{mt},newx,newy,1); %grow straight if picked up ER
            end
        else % if did not pick up ER
            if size(mts{mt},1)>1
                [~,newx,newy]=checkconn(mts{mt},newx,newy,0);
                if length(newx)>1 && length(newy)>1
                    newxo=newx;newx=0;
                    newyo=newy;newy=0;
                    while newx==0 && newy==0
                        newx=randi(newxo,1);
                        newy=randi(newyo,1);
                    end
                else
                if length(newx)>1; newx=randi(newx,1);end
                if length(newy)>1; newy=randi(newy,1); end
                end
            end
        end 
        
        if mts{mt}(end,1)+newx==mts{mt}(end-1,1) && mts{mt}(end,2)+newy==mts{mt}(end-1,2) && isempty(mts{mt})==0 %shrink
            mts{mt}(end,:)=[];
            if mtsER(mt)==1 && ~isempty(find(conn2mt==mt))
               bkgER(mts{mt}(end,1),mts{mt}(end,2))=0; 
               ERs{find(conn2mt==mt)}(end,:)=[];
               ERpool=ERpool+1;
            elseif mtsER(mt)==1
                bkgER(mts{mt}(end,1),mts{mt}(end,2))=1;
                ERs{length(ERs)+1}=[mts{mt}(end,1),mts{mt}(end,2)];
                conn2mt(length(ERs))=mt;
                ERpool=ERpool-1;
            end
            pool=pool+1;
        else
            mts{mt}(end+1,:)=[mts{mt}(end,1)+newx,mts{mt}(end,2)+newy];
            pool=pool-1;
            pos2=cellfun(@(x) ismember(mts{mt},x,'rows'),mts,'UniformOutput',false,'ErrorHandler',@errorfun);
            pos3=cellfun(@(x) length(find(x(2:end)~=0)),pos2);
            pos3(mt)=[];
            if bkgER(mts{mt}(end,1),mts{mt}(end,2))==1 && isempty(find(pos3>5)) && w3j(mts{mt}(end,1),mts{mt}(end,2))~=1 % can fuse if new MT part touches ER edge and not 3wj and they are not parallel
                canfuse=1; end;
            
            if ~isempty(find(conn2mt==mt))
                ERs{find(conn2mt==mt)}(end+1,:)=[mts{mt}(end,1),mts{mt}(end,2)];
                bkgER(mts{mt}(end,1),mts{mt}(end,2))=1;
                ERpool=ERpool-1;
            elseif mtsER(mt)==1 && isempty(find(conn2mt==mt))
                bkgER(mts{mt}(end,1),mts{mt}(end,2))=1;
                ERs{length(ERs)+1}=[mts{mt}(end,1),mts{mt}(end,2)];
                conn2mt(length(ERs))=mt;
                ERpool=ERpool-1;
            end            
            if canfuse==1 && mtsER(mt)==1 %if meets ER looses the connection
                if randi([0,100],1)<probfusion % probability to fuse
                    w3j(mts{mt}(end,1),mts{mt}(end,2))=1;
                    mtsER(mt)=0;
                    conn2mt(conn2mt==mt)=0;
                    l=cellfun(@length, ERs);
                end  
            end
        end
    end 
  %  cleanANDdraw(mtsER,mts,conn2mt,ERs,ss,1);
end
[mtsER,mts,conn2mt,ERs,bkgER,bkg,AcMT]=cleanANDdraw(mtsER,mts,conn2mt,ERs,ss,AcMT,0);
%pause(0.3)
end
im2=imfuse(bkg,bkgER,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
if saveimg==0
    figure('name',sprintf('prob fusion: %d prob collapse: %d', probfusion, probcollapse))
    imagesc(im2)
    colormap jet;
else
    imwrite(im2,sprintf('pFusion_%d_pCollapse_%d_pPickUp_%d.png',probfusion, probcollapse,probpickup),'BitDepth',8)
end
end

function [out,newx,newy]=checkconn(line,newx,newy, tocken)
out=0;
if newx==0 && newy==0
    out=0;
else
if size(line,1)<5
        newx=line(end,1)-line(end-1,1);
        newy=line(end,2)-line(end-1,2);
else
    if tocken==0 && ((abs(line(end,1)-line(end-1,1))==1 && abs(line(end-1,1)-line(end-2,1))==1 && abs(line(end,2)-line(end-1,2))==1 && abs(line(end-1,2)-line(end-2,2))==1) ...
            || (line(end,1)-line(end-1,1)==0 && line(end-1,1)-line(end-2,1)==0) || (line(end,2)-line(end-1,2)==0 && line(end-1,2)-line(end-2,2)==0))
           if line(end,1)-line(end-1,1)==1 && line(end,2)-line(end-1,2)==0
               newx=1; newy=[-1 1];
  
           elseif line(end,1)-line(end-1,1)==-1 && line(end,2)-line(end-1,2)==0
               newx=-1; newy=[-1 1];

           elseif line(end,2)-line(end-1,2)==1 && line(end,1)-line(end-1,1)==0
               newy=1; newx=[-1 1];
  
           elseif line(end,2)-line(end-1,2)==-1 && line(end,1)-line(end-1,1)==0
               newy=-1; newx=[-1 1];
     
           elseif line(end,2)-line(end-1,2)==1 && line(end,1)-line(end-1,1)==1
               newy=[0 1]; newx=[0 1];
      
           elseif line(end,2)-line(end-1,2)==-1 && line(end,1)-line(end-1,1)==1
               newy=[-1 0]; newx=[0 1];

           elseif line(end,2)-line(end-1,2)==1 && line(end,1)-line(end-1,1)==-1
               newy=[0 1]; newx=[-1 0];
       
           elseif line(end,2)-line(end-1,2)==-1 && line(end,1)-line(end-1,1)==-1
               newy=[-1 0]; newx=[-1 0];
            
           end  
    else
        newx=line(end,1)-line(end-1,1);
        newy=line(end,2)-line(end-1,2); 
    end
end
end
end

function border=DetectBorders(bkgER)
% bkgERout=imfill(bkgER,'holes');
% B1=bwboundaries(bkgERout);
% bkgERinv=imclose(imcomplement(bkgER),[1 1 ;1 1]);
B=bwboundaries(bkgER);
border=zeros(size(bkgER,1),size(bkgER,2));
% border(sub2ind(size(border),B1{1}(:,1),B1{1}(:,2)))=1;
for i=1:length(B)
    border(sub2ind(size(border),B{i}(:,1),B{i}(:,2)))=1;
end
end

function [mtsER,mts,conn2mt,ERs,bkgER,bkg,AcMT]=cleanANDdraw(mtsER,mts,conn2mt,ERs,ss,AcMT,draw)
    bkg=zeros(ss,ss);
    bkgER=generatebkgER(ss,10);
    %---------------MTs
    idxbefore=find(mtsER);
    idx=cellfun(@isempty,mts);
    mts(idx)=[];
    mtsER(idx)=[];
    AcMT(idx)=[];
    for mt=1:length(mts)
        idx=sub2ind(size(bkg),mts{mt}(:,1),mts{mt}(:,2));
        bkg(idx)=1;
    end
    %---------------ERs
    idxafter=find(mtsER);
    for i=1:length(idxbefore)
        conn2mt(conn2mt==idxbefore(i))=idxafter(i);
    end
    idx=cellfun(@isempty,ERs);
    ERs(idx)=[];
    conn2mt(idx)=[];
    if ~isempty(ERs)
        for ertub=1:length(ERs)
            idx=sub2ind(size(bkgER),ERs{ertub}(:,1),ERs{ertub}(:,2));
            bkgER(idx)=1;
        end
    end
    
    %---------------Drawing
    if draw==1
        im2=imfuse(bkg,bkgER,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        imagesc(im2)
        colormap jet;
    end
    drawnow;
end

function bkgER=generatebkgER(ss,r)
bkgER=zeros(ss,ss);
 if r==0     
     bkgER(ss/2,ss/2)=1;
  else
     x = 1:ss;
     y = 1:ss;
     [xx, yy] = meshgrid(x,y);
     bkgER(((xx-ss/2-r*(2/3)).^2+(yy-ss/2-r*(2/3)).^2)<r^2)=1;
     bkgER(((xx-ss/2-r*(2/3)).^2+(yy-ss/2-r*(2/3)).^2)<(r-1)^2)=0;
 end
end
function result = errorfun(S, varargin)
   warning(S.identifier, S.message);
   result = 0;
end
