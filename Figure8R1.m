% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% Figure 8
%
% Victor M. Vergara, PhD
% Created 2018-August-3
% 
% 2018-October-10: 
%	- Corrected the membership functions displayed. We are comparing tuned 
%	  SWC with tuned ASWC. SWC and ASWC results were run trough k-means and
%	  the membership functions displayed correspond to these two cases.
%	- We matched the dFNC states with those in (Vergara 2018). The k-means
%	  membership functions were modified to match the indexes of mentioned
%	  dFNC states.
%	- The correct tuning thus is:
%		f = 0.01 Hz;
%		TR= 2 sec;
%		1/(2f) = 50 sec
%		g_aswc = 1/(2f) = 50 sec
%		h_aswc = 44 sec
%		g_aswc + h_aswc = 94 sec
%
% References:
%
% Vergara, V.M., Weiland, B.J., Hutchison, K.E., Calhoun, V.D. (2018).
% The Impact of Combinations of Alcohol, Nicotine, and Cannabis on Dynamic
% Brain Connectivity. Neuropsychopharmacology, 43:877-890.
% 
%
%% Initialize
clear;close all;
%% Some constants
TR = 2; % sec
cmp = parula(8);
SWCLen = 100; % tuned  SWC (UsedLen = h = 100 sec) : membership plots
ASWCLen = 96; % tuned ASWC (UsedLen =g+h=  96 sec) : membership plots
SWCMarker = '-.';
ASWCMarker= '-';
Orange = [0.9 0.35 0];
% These are exclusive for ASWC
half_ave_len_aswc = 25;	% Averaging Length is 50 sec (tune 0.01Hz)
half_win_len_aswc = 22;  % Window length is 44 sec
nabla = floor(half_ave_len_aswc/TR);
% These are exlusive for SWC
half_win_len_swc = 50;  % Window length is 100 sec
% FontSizes
axisFS = 16;
%% Initialize the Figure
figure('rend','painters','pos',[10 100 900 600]);

%% Plot Constant Connectivity. Example Subject 37
subplot(3,3,1);
% ----------- SWC --------------
load Figure8_membership_sub037_SWC;
t = 2*(0:1:length(m)-1) + floor(SWCLen/2); % temporal offset
m(m==6)=2;	% Match the State Index
plot(t,m,SWCMarker,'Color',Orange,'LineWidth',4);
% -- Configure Subplot--
A = [0 TR*length(m)+SWCLen 0 7];
axis(A);
set(gca,'yticklabels',{'','1','2','3','4','5','6',''},'FontSize',axisFS);
set(gca,'ytick',[0 1 2 3 4 5 6 7])
title('Static Case','FontSize',18);
xlabel('time [sec]','FontSize',axisFS);
ylabel('dFNC state','FontSize',axisFS);
grid on;
hold on;
% ----------- ASWC --------------
load Figure8_membership_sub037_ASWC;
t = 2*(0:1:length(m)-1) + floor(ASWCLen/2); % temporal offset
m(m==1)=2;	% Match the State Index
plot(t,m,ASWCMarker,'Color',cmp(2,:),'LineWidth',2);

%% Sharp Transition Subject 248
subplot(3,3,4);
% Sharp transitions must be synchornized with the window length
load Figure8_membership_sub248_SWC;
t = 2*(0:1:length(m)-1) + floor(SWCLen/2); % temporal offset
m(m==1) = 1;% Match the State Index
m(m==3) = 2;% Match the State Index
plot(t,m,SWCMarker,'Color',Orange,'LineWidth',3);
% -- Configure Subplot--
A = [0 TR*length(m)+SWCLen 0 7];
axis(A);
set(gca,'yticklabels',{'','4','6','','','4','6',''},'FontSize',axisFS);
set(gca,'ytick',[0 1 2 3 4 5 6 7])
title('Sharp Edge Transition','FontSize',18);
xlabel('time [sec]','FontSize',axisFS);
ylabel('dFNC state','FontSize',axisFS);
grid on
hold on;
% ----------- ASWC --------------
load Figure8_membership_sub248_ASWC;
t = 2*(0:1:length(m)-1) + floor(ASWCLen/2); % temporal offset
m(m==2) = 5;% Match the State Index and the displayed y-axis
m(m==6) = 6;% Match the State Index and the displayed y-axis
plot(t,m,ASWCMarker,'Color',cmp(2,:),'LineWidth',3);
%% 100 sec period Subject 41
subplot(3,3,7);
% ----------- SWC --------------
load Figure8_membership_sub041_SWC;
t = 2*(0:1:length(m)-1) + floor(SWCLen/2); % temporal offset
m(m==4) = 1;
m(m==2) = 2;
plot(t,m,SWCMarker,'Color',Orange,'LineWidth',3);
% -- Configure Subplot--
A = [0 TR*length(m)+SWCLen 0 7];
axis(A);
title('Fluctuating at approx. 0.012 Hz','FontSize',18);
xlabel('time [sec]','FontSize',axisFS);
ylabel('dFNC state','FontSize',axisFS);
set(gca,'yticklabels',{'','3','5','','','3','5',''},'FontSize',axisFS);
set(gca,'ytick',[0 1 2 3 4 5 6 7])
grid on;
hold on;
% ----------- ASWC --------------
load Figure8_membership_sub041_ASWC;
t = 2*(0:1:length(m)-1) + floor(ASWCLen/2); % temporal offset
A = [0 TR*length(m)+ASWCLen 0 7];
m(m==5) = 6;
m(m==3) = 5;
plot(t,m,ASWCMarker,'Color',cmp(2,:),'LineWidth',3);

%% Constant Connectivity : temporal window estimation
load Figure8_WL44sec_gica_dfnc_sub_037;
selectedRSNpair = [53 82];
subplot(3,3,2);
dfnc = FNCdyn(:,495); % ica comp 53 82
% --------------- ASWC ---------------
p=0;
for kk=nabla+1:size(dfnc,1)-nabla
	p=p+1;
	% notice that length(kk-nabla+1:kk+nabla) = 25 samples: 25*TR = 50 sec
	% this corresponds to the setting g = 50 sec 
	avef(p) = mean(dfnc(kk-nabla:kk+nabla));
end
t = TR*(nabla:size(dfnc,1)-nabla-1) + half_win_len_aswc;
plot(t,avef,ASWCMarker,'Color',cmp(2,:),'LineWidth',3);
set(gca,'FontSize',axisFS);
A(3)=0.6;A(4)= 1.4;axis(A);
xlabel('time [sec]','FontSize',axisFS);
ylabel('atanh (corr)','FontSize',axisFS);
title('h_{ASWC}=44 sec  g_{ASWC}=50 sec','FontSize',14);
% other details
grid on;
legend(['\sigma = ',sprintf('%.2f',std(avef))]);
% --------------- SWC ---------------
load Figure8_WL100sec_gica_dfnc_sub_037;
subplot(3,3,3);
dfnc = FNCdyn(:,495); % ica comp 53 82
t = TR*(1:1:length(dfnc)) + half_win_len_swc; % WL 100 we adda on offset
plot(t,dfnc,SWCMarker,'Color',Orange,'LineWidth',3);
set(gca,'FontSize',axisFS);
A(3)=0.6;A(4)= 1.4;axis(A);
xlabel('time [sec]','FontSize',axisFS);
ylabel('atanh (corr)','FontSize',axisFS);
title('h_{SWC}= 100 sec','FontSize',14);
% other details
grid on;
legend(['\sigma = ',sprintf('%.2f',std(dfnc))]);

%% Sharp Edige Transition Connectivity : temporal window estimation
load Figure8_WL44sec_gica_dfnc_sub_248;
subplot(3,3,5);
dfnc = FNCdyn(:,381); % ica comp 16 51
% --------------- ASWC ---------------
p=0;
clear avef;
for kk=nabla+1:size(dfnc,1)-nabla
	p=p+1;
	% notice that length(kk-nabla+1:kk+nabla) = 26 samples: 26*TR = 52 sec
	% this corresponds to the setting g = 50 + TR = 52	
	avef(p) = mean(dfnc(kk-nabla:kk+nabla));
end
t = TR*(nabla:size(dfnc,1)-nabla-1) + half_win_len_aswc;
plot(t,avef,ASWCMarker,'Color',cmp(2,:),'LineWidth',3);
set(gca,'FontSize',axisFS);
A(3)= -0.20;A(4)= 0.8;axis(A);
xlabel('time [sec]','FontSize',axisFS);
ylabel('atanh (corr)','FontSize',axisFS);
title('h_{ASWC}=44 sec  g_{ASWC}=50 sec','FontSize',14);
% other details
grid on;
% --------------- SWC ---------------
load Figure8_WL100sec_gica_dfnc_sub_248;
selectedRSNpair = [53 82];
subplot(3,3,6);
dfnc = FNCdyn(:,381); % ica comp 16 51
t = TR*(1:1:length(dfnc)) + half_win_len_swc; % WL 100 we adda on offset
plot(t,dfnc,SWCMarker,'Color',Orange,'LineWidth',3);
set(gca,'FontSize',axisFS);
A(3)= -0.20;A(4)= 0.8;axis(A);
xlabel('time [sec]','FontSize',axisFS);
ylabel('atanh (corr)','FontSize',axisFS);
title('h_{SWC}= 100 sec','FontSize',14);
% other details
grid on;

%% High Freq. Varying Connectivity : temporal window estimation
load Figure8_WL44sec_gica_dfnc_sub_041;
subplot(3,3,8);
dfnc = FNCdyn(:,521); % ica comp 6 25
% --------------- ASWC ---------------
p=0;
clear avef;
for kk=nabla+1:size(dfnc,1)-nabla
	p=p+1;
	% notice that length(kk-nabla+1:kk+nabla) = 25 samples: 25*TR = 50 sec
	% this corresponds to the setting g = 50 sec	
	avef(p) = mean(dfnc(kk-nabla:kk+nabla));
end
t = TR*(nabla:size(dfnc,1)-nabla-1) + half_win_len_aswc;
plot(t,avef,ASWCMarker,'Color',cmp(2,:),'LineWidth',3);
set(gca,'FontSize',axisFS);
A(3)= -0.6;A(4)= 0.3;axis(A);
xlabel('time [sec]','FontSize',axisFS);
ylabel('atanh (corr)','FontSize',axisFS);
title('h_{ASWC}=44 sec  g_{ASWC}=50 sec','FontSize',14);
% other details
grid on;
% --------------- SWC ---------------
load Figure8_WL100sec_gica_dfnc_sub_041;
subplot(3,3,9);
dfnc = FNCdyn(:,521); % ica comp 6 25
t = TR*(1:1:length(dfnc)) + half_win_len_swc; % WL 100 we adda on offset
plot(t,dfnc,SWCMarker,'Color',Orange,'LineWidth',3);
set(gca,'FontSize',axisFS);
A(3)= -0.60;A(4)= 0.3;axis(A);
xlabel('time [sec]','FontSize',axisFS);
ylabel('atanh (corr)','FontSize',axisFS);
title('h_{SWC}= 100 sec','FontSize',14);
% other details
grid on;

