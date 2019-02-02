% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% Average Sliding Window Correlation (ASWC)
% Figure 5
%
% This code plots the frequency responses of SWC and ASWC methods.
%
% Victor M. Vergara, PhD
% Created on 2018-July-15
%
clear;close all;clc;

TR = 1;
N = 1000;

f = 0:1/N:0.05;
FZ = 18;

hswc = 100 ; % in sec

haswc=  44;  % in sec
gaswc=  50;  % in sec

FH = figure('rend','painters','pos',[128 128 900 300]);

subplot(1,2,1);
x = pi*TR*f;
y = (1/hswc)*abs(sin(x*hswc)./sin(x));
plot(f,y,'LineWidth',2);
title('SWC Frequency Response','FontSize',FZ);
xlabel('Frequency Hz','FontSize',FZ);
legend('H_h: h_{SWC} = 100 sec');
set(gca,'FontSize',14);
axis([0 0.05 0 1]);
grid on;


subplot(1,2,2);
x = pi*TR*f;
y1 = (1/haswc)*abs(sin(x*haswc)./sin(x));
y2 = (1/gaswc)*abs(sin(x*gaswc)./sin(x));
y = y1.*y2;
plot(f,y1,'.',f,y2,'--',f,y,'LineWidth',2)
title('ASWC Frequency Response','FontSize',FZ);
xlabel('Frequency Hz','FontSize',FZ);
legend('H_h: h_{ASWC} = 44 sec',['H_g: g_{ASWC} = ',	...
	num2str(gaswc),' sec'],'ASWC: H_h H_g');
set(gca,'FontSize',14);
axis([0 0.05 0 1]);
grid on;

%% Save the figure

saveas(FH,'Figure5R1','png');
saveas(FH,'Figure5R1','fig');