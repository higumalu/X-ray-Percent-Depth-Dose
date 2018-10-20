clear all;
close all;

img_dir='C:\Users\HIGUMA_LU\Desktop\放射治療實驗\tif';
data_list=dir([img_dir '\*.tif']);

for n=1:length(data_list) 

im=imread([img_dir '\' data_list(n).name]);%6MV X-ray 照野10x10

imall = double(im);
imr = imall(:,:,1);
[a1 b1] = size (imr); % Y, X

xc = b1/2; % Central axis position
lfb = 0.1/0.02; % 0.1 cm / 0.02 cm; 0.02 cm is the resolution  % 1pix=0.2mm=0.02cm  1mm=0.1cm
udb = 0.1/0.02; %up down
step = 0.2; %%% Step is 0.2 cm 每2mm 取1點
for i = 1:25.3/step % The last column 25.4 cannot be used, since no 25.5 for average
pv_mean_ud(i,:) = mean(imr(int32(i*step/0.02-udb+1):int32(i*step/0.02+udb),:));% Combine up and down for central axis
    for j = -9/step : 9/step  % Left Right, Choose -9 cm to 9 cm
    pv_mean(i,j+9/step+1) =  mean(pv_mean_ud(i,int32(xc+j*step/0.02-lfb+1):int32(xc+j*step/0.02+lfb)));% % This is optical density now; Combine left and right for the certain depth
    end
end



%%% Data Output for Curve Fitting %%%
load dose300cGy;
pv = dose300cGy(:,2);   %pixel value
ds = dose300cGy(:,1);   %dose
[p, s] = polyfit(pv,ds,3);  % Fitting
%%%%%%%%%%% Transfer to Pixel Value to Dose %%%%%%%%%%
imd = p(1)*pv_mean.^3+p(2)*pv_mean.^2+p(3)*pv_mean+p(4);

%%%%%%%%% Get PDD at Center %%%%%%%%%%%%%%%%%%
[h2, w2] = size(imd);

w_cen60=uint8(w2*0.2:w2*0.8);
w_cen50=uint8(w2*0.25:w2*0.75);
imd_60=imd(:,w_cen60);
imd_50=imd(:,w_cen50);

[h_60, w_60]=size(imd_60); %    L&R-20%
[h_50, w_50]=size(imd_50);% L&R-25%

depth_dose = imd(:, int16(w2/2));%find dose-pic center
depth_dose_all(:,n)=depth_dose;
file_name=data_list(n).name

d_max = max(depth_dose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fitx = []; fity = []; %fitcy = [];
for x = pv(1)-1000:1:pv(length(pv))+1000
  fitx = [fitx;x];
  fity = [fity; p(1)*x^3+p(2)*x^2+p(3)*x+p(4)];
end



%%%make PDD %%%
depth = 0.2:0.2:25.2;
%figure; plot(depth, depth_dose, '*'); 
figure; plot(depth,  depth_dose/max(depth_dose)*100);
title (['PDD--' data_list(n).name]);
PDD = depth_dose/max(depth_dose)*100;
%t = imd(10/0.2,:);
%figure; plot(t)
depth_N = 1:1:25; 

depth_dose_N(:,n)= imd(depth_N/0.2, int16(w2/2));

depth = 0.2:0.2:25.2;

%plot(depth, depth_dose, '*');


%%%savefi
%{
%saveas(PDD1,['PDD1' int2str(n) '.png'])
PDD=figure;
plot(pv, ds, 'rx', fitx, fity, 'g');
%saveas(PDD,['PDD' int2str(n) '.png'])
PDD2=figure; plot(depth,  depth_dose/max(depth_dose)*100);
title ('PDD');
%saveas(PDD2,['PDD2' int2str(n) '.png'])
PDD = depth_dose/max(depth_dose)*100;
t = imd(15,:);
FWD=figure; plot(t)
%saveas(FWD,['FWD' int2str(n) '.png'])
%}


%%%get dose at 3,10,20cm%%%
t1 = imd(15,:);
t2 = imd(50,:);
t2_60=imd_60(50,:);
t2_RD=uint8(t2*100/max(t2));
t3= imd(100,:);
rand_d=randi([3 w_60/2.5],1,1);


%%%Calc field size%%%
FD=max(t2)/2;
%{
for i=1:1:uint8((w2-1)/2)
    if t2_RD==50
        %[i uint8((w2-1)/2)+i t2_RD(i) t2_RD(uint8((w2-1)/2)+i)];
        FS=i
        
    end
end
%}


%%%Calc flatness%%%
flatness=(max(t2_60)-min(t2_60))/(max(t2_60)+min(t2_60))*100

%%%Calc symmetry%%%
D1=imd_60(10,int16(w_60/2)+rand_d);
D2=imd_60(10,int16(w_60/2)-rand_d);
if D1-D2 < 0
    symmetry=(D2-D1)/D1*100

else
    symmetry=(D1-D2)/D1*100
end
    
%%%FWD ALL%%%
figure;
x=-9:18/(w2-1):9;
xFD=zeros(w2,1);
xFD(:)=FD;
hold on;
plot(x,t1)
plot(x,t2)
plot(x,xFD,'--')
plot(x,t3)
hold off;
set(gca,'xtick',(-10:1:10)) 
title(['FWD--' data_list(n).name])
xlabel('width(cm)'), ylabel('dose(cGy)');
legend('3cm','10cm','10cm 50% dose level','20cm');
[x y]=ginput(2);
fieldsize=((x(2)-x(1))^2+(y(2)-y(1))^2)^0.5
end


%%%plot pixvel value & dose%%%
figure; plot(pv, ds, 'rx', fitx, fity, 'g');

%%%PDD ALL%%%
depth = 0.2:0.2:25.2;
PDD1=figure; 
hold on;
title('P D D');
xlabel('depth(cm)'), ylabel('dose(cGy)');

plot(depth, depth_dose_all(:,1), '.');
plot(depth, depth_dose_all(:,2), '.');
plot(depth, depth_dose_all(:,3), '.');
plot(depth, depth_dose_all(:,4), '.');
hold off;
legend(data_list(1).name,data_list(2).name,data_list(3).name,data_list(4).name);

save('depth_dose_N.mat','depth_dose_N')



%b11128

