SPMaddres='D:\MATLAB\spm8';
addpath(SPMaddres)
addpath(strcat(SPMaddres,'\matlabbatch'))
addpath(strcat(SPMaddres,'\matlabbatch\cfg_basicio'))
addpath(strcat(SPMaddres,'\config'))
addpath(strcat(SPMaddres,'\toolbox\DARTEL'))
addpath(strcat(SPMaddres,'\toolbox\FieldMap'))
addpath(strcat(SPMaddres,'\toolbox\Seg'))
%% dicom import
clear matlabbatch
load('dicom import.mat')

pp = sprintf('D:\\gaba_data\\20141210_PengBY_svs\\scan2\\MPR');

x=dir('*.IMA');
y={x.name};
for i=1:176
   nn=sprintf('%s\\%s',pp,y{i});
   tempcell{i,1}=nn;
end

matlabbatch{1,1}.spm.util.dicom.data=tempcell;
matlabbatch{1,1}.spm.util.dicom.outdir={[pp]};
matlabbatch{1,1}.spm.util.dicom.convopts.format = 'img';
matlabbatch{1,1}.spm.util.dicom.convopts.icedims = 0;
spm_jobman('run',matlabbatch); 

%% segment
clear matlabbatch
load('segment1.mat')
x=dir('*-000176-01*.img');
x=x.name;
y=sprintf('%s\\%s,1',pp,x);
matlabbatch{1,1}.spm.spatial.preproc.data={[y]};
spm_jobman('run',matlabbatch); 
%% normalise
clear matlabbatch
load('normalise1.mat')
% specify the transformation you want to apply

gg =  sprintf('%s\\SMA_voxel.nii',pp);

x=dir('*inv*.mat');
x=x.name;
y=sprintf('%s\\%s',pp,x); %自行修改路徑
matlabbatch{1}.spm.spatial.normalise.write.subj.matname={[y]};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample ={[gg]};
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-255 -255 -255;255 255 255];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%% load 對位後voxel 取得position
roi=MRIread('wvc_voxel.nii'); 
a = roi.vol;
a(isnan(a))=0;
b=logical(a);
stats = regionprops(b,'basic');
cen = stats.Centroid;
cennew = [cen(2)-128.5 cen(1)-128.5 cen(3)-128.5] * 2;
fprintf(' %5.1f  %5.1f  %5.1f\n', cen(3),cen(2),cen(1))
fprintf('L %5.1f P %5.1f H %5.1f\n',-cennew(3),cennew(2),cennew(1))
%%
cenf = floor(cen);
S_z=(flipud(squeeze(b(:,:,cenf(3)))));
IM1 = bwmorph(S_z,'remove');
%IM1=edge(S_z,'canny');
figure;imshow(IM1);
%%
IM1=double(IM1);
line1=zeros(256,256);
line2=zeros(256,256);

u=0;
 for j= 1:256
   for i= 1:256
     if (IM1(i,j) == 1)
        y=i;
        x=j;
        u=1;
        break
     end
   end
   if(u==1)
       break
   end
end
x1=x;
y1=y;
for j = 1:256
   for i =1:256
       if(IM1(i,j)==1)    
        if(j>=x1 && i>=y1)
           y1=i;
           x1=j;
           %a=polyfit(y1,x1,1);
           line1(y1,x1)=IM1(y1,x1);
        end
       end
   end
end  
figure(1);imshow(line1);

m1 = zeros(256);
p1 = [x,y];
p2 = [x1,y1];
x = p1(1):p2(1);
y =ceil((x - p1(1)) * (p2(2) - p1(2)) / (p2(1) - p1(1)) + p1(2));
m1(sub2ind(size(m1),y,x)) = 1;
figure;imshow(m1);

A=[x(1)-1,256-y(1),cenf(3)];
B=[x(5)-1.5,256-y(5),cenf(3)];
% CC=[x(5)-1.5,256-y(5),cenf(3)+5];
%%
% VAB=B-A;
% VCCB=B-CC;
% D=cross(VAB,VCCB);
% 
% Vx=[1,0,0];
% Vy=[0,1,0];
% Vz=[0,0,1];
% 
% Dx = dot(D,Vx);
% Dy = dot(D,Vy);
% Dz = dot(D,Vz);
% 
% deg_raw = atan(Dz/Dy)*180/pi;
% deg_pit = atan(Dz/Dx)*180/pi;
% deg_yaw = atan(Dy/Dx)*180/pi;
% 
% fprintf('L %f P %f H %f\n',-cennew(3),cennew(2),cennew(1))
% %fprintf('T > C %f\n',90 + deg_yaw);
% 
%  fprintf('MPFC\n');
% 
%  if (90+deg_yaw) > 45
%     fprintf('C > T %f\n',-deg_yaw);
%     fprintf('T > C %f\n',-(90-deg_yaw));
%  else
%     fprintf('C > T %f\n',90 + deg_yaw);  
%  end  
 
%%
S_y=((squeeze(b(256-y(5),:,:))));
IM2 = edge(S_y,'canny',0.01,0.1);
% IM3 = bwmorph(S_y,'remove');
figure;imshow(IM2);
% figure;imshow(IM3);

xxx = 256-y(5);

line3=zeros(256,256);
line4=zeros(256,256);

u=0;
for i= 1:256
   for j= 1:256
     if (IM2(i,j) == 1)
        yy=i;
        xx=j;
        u=1;
        break
     end
   end
   if(u==1)
       break
   end
end
xx1=xx;
yy1=yy;
for i = 1:256
   for j =1:256
       if(IM2(i,j)==1)    
        if(j>xx1 && i>=yy1)
           yy1=i;
           xx1=j;
           %a=polyfit(y1,x1,1);
           line3(yy1,xx1)=IM2(yy1,xx1);
        end
       end
   end
end  
figure(1);imshow(line3);

m2 = zeros(256);
p1 = [xx,yy];
p2 = [xx1,yy1];
xx = p1(1):p2(1);
yy =ceil((xx - p1(1)) * (p2(2) - p1(2)) / (p2(1) - p1(1)) + p1(2));
m2(sub2ind(size(m2),yy,xx)) = 1;
figure;imshow(m2);

if  yy(2) ~= yy(5);
    yy(1) = yy(1);
else
    yy(2) = B(1);
end   

C=[yy(2),xxx,xx(2)];
%%
% degg =  deg_yaw * pi/180;
% deg_sin = sin(degg);
% deg_cos = cos(degg);
%%
VAB=B-A;
VCB=B-C;
D = cross(VAB,VCB);

Vx=[1,0,0];
Vy=[0,1,0];
Vz=[0,0,1];

Dx = dot(D,Vx);
Dy = dot(D,Vy); 
Dz = dot(D,Vz);

raw = atan(Dz/Dy)*180/pi;
yaw = atan(Dz/Dx)*180/pi;
pit = atan(Dy/Dx)*180/pi;

fprintf('VC\n',pit);
fprintf('T > C %f\n',pit);

 if (pit > 45)
    fprintf('T > C %f\n',pit);
 else
    fprintf('C > T %f\n',pit);
    fprintf('T > C %f\n',-pit);
 end  
    fprintf('C > S %f\n',yaw);
%     fprintf('C > S %f\n',yaw+90);
