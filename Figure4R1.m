% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% Average Sliding Window Correlation (ASWC)
% Figure 4
%
% Victor M. Vergara, PhD
% Created on 2018-June-20
%
%
clear;close all;
%% Setup the simulation configuration
TR =   1;	% [sec]
halfSimLenSec = 100;% [sec]

F = 0.025;		% Cosine Frequency
A=sqrt(2);		% Cosine Amplitude

% IN this code we dont want to change the window, so set one size
winLenSec = 40;		% Window Length [sec]
winLen = ceil(winLenSec/TR);

aveLenSec = 20;
aveLen = floor(aveLenSec/TR);

halfnTR = ceil(halfSimLenSec/TR);
nTR = 1 + 2*halfnTR;	% How many points
FSZ = 18; %Figure font size

% For this simulation we have 4 cosines to simulate. 
% So the problem of defining the phases between them is left to the user.
corr1 =  0.5;	% you can use these two correlations and check the phases
corr2 = -0.5;	% 
Dphase1 = acos(corr1);
Dphase2 = acos(corr2);
phase1 = pi/4; %	[rads]
phase2 = phase1 - Dphase1;
phase3 = -pi/4; %	[rads]
phase4 = phase3 - Dphase2;

%% Create the cosines to be estimated
TCOPC = (halfnTR + 1); % timeCourseOfPhaseChange
t = ((1:nTR) - TCOPC) * TR;
x=A*[cos(2*pi*F*t(1:TCOPC-1) + phase1),cos(2*pi*F*t(TCOPC:end) + phase3)];
y=A*[cos(2*pi*F*t(1:TCOPC-1) + phase2),cos(2*pi*F*t(TCOPC:end) + phase4)];

%% Display the Simulation Configuration
fprintf('Sliding Window Simulation (Cosine Model)\n');
fprintf('----------------------------------------\n');
fprintf('\n\t     TR [sec] = %.1f\n',TR);
fprintf('\tNumber of TRs = %d\n',nTR);
fprintf('\tDuration[sec] = %.1f\n',TR*nTR);
fprintf('\tWindow Duration[sec] = %.1f\n',winLenSec);
fprintf('\n');
fprintf('Parameters for the 2 Simulated Cosines\n');
fprintf('\t [Hz] \t Amp \t\n');
for kk=1:length(F)
	fprintf('\t%.3f\t%.3f\n',F(kk),A(kk));
end

cc1 = corr1;
cc2 = corr2;
fprintf('Interval 1 Correlation = %.3f\n',cc1);
fprintf('Interval 2 Correlation = %.3f\n',cc2);

%% Plot a) The Cosines
FH = figure('rend','painters','pos',[10 100 900 600]);
subplot(2,2,1);
plot(t,x,'-o',t,y,'--x','LineWidth',1);legend('x','y');
xlabel('Time [sec]','FontSize',FSZ)
ylabel('Cosine Amplitude','FontSize',FSZ);
title('a) Sharp Phase Transition','FontSize',FSZ);
%% Plot b) 
% Covariance by Sliding Window
for tc = 1:nTR-winLen
	cm = cov(x(tc:tc+winLen-1)',y(tc:tc+winLen-1)');
	cc(tc) = cm(1,2);
end
% Create a covariance averaging
tt = TR*((1:length(cc)) - floor(length(cc)/2));
for tc = 1:length(cc)-aveLen
	at(tc) = tt(tc+floor(aveLen/2));
	ac(tc) = mean(cc(tc:tc+aveLen-1));
end
% Plot the results
subplot(2,2,2);
plot(tt,cc,at,ac,[-100 0],cc1*[1 1],	...
	'--',[0 100],cc2*[1 1],'--','LineWidth',2);
grid on;
A = axis();A(1) = min(tt)-TR;A(2) = max(tt)+TR;axis(A);
xlabel('Time [sec]','FontSize',FSZ);
ylabel('Covariance','FontSize',FSZ);
legend('SWC','ASWC',['corr 1 =',num2str(cc1,1)],	...
	['corr 2 =',num2str(cc2,1)]);
title(['b) h=',num2str(winLen),'   g=',num2str(aveLen)],'FontSize',FSZ);

%% Plot c) Double averaging 
aveLenOriginal =aveLen;
aveLen = 2*aveLen;
clear cc; clear cm;clear at;clear ac;
% Covariance by Sliding Window
for tc = 1:nTR-winLen
	cm = cov(x(tc:tc+winLen-1)',y(tc:tc+winLen-1)');
	cc(tc) = cm(1,2);
end
% Create a covariance averaging
tt = TR*((1:length(cc)) - floor(length(cc)/2));
for tc = 1:length(cc)-aveLen
	at(tc) = tt(tc+floor(aveLen/2));
	ac(tc) = mean(cc(tc:tc+aveLen-1));
end
% Plot the results
subplot(2,2,3);
plot(tt,cc,at,ac,[min(tt) 0],cc1*[1 1],	...
	'--',[0 max(tt)],cc2*[1 1],'--','LineWidth',2);
grid on;
A = axis();A(1) = min(tt)-TR;A(2) = max(tt)+TR;axis(A);
xlabel('Time [sec]','FontSize',FSZ);
ylabel('Covariance','FontSize',FSZ);
legend('SWC','ASWC',['corr 1 =',num2str(cc1,1)],	...
	['corr 2 =',num2str(cc2,1)]);
title(['c) h=',num2str(winLen),'   g=',num2str(aveLen)],'FontSize',FSZ);

%% Plot d) Double windowing
aveLen = aveLenOriginal;
winLen = 2*winLen;
clear cc; clear cm;clear at;clear ac;
% Covariance by Sliding Window
for tc = 1:nTR-winLen
	cm = cov(x(tc:tc+winLen-1)',y(tc:tc+winLen-1)');
	cc(tc) = cm(1,2);
end
% Create a covariance averaging
tt = TR*((1:length(cc)) - floor(length(cc)/2));
for tc = 1:length(cc)-aveLen
	at(tc) = tt(tc+floor(aveLen/2));
	ac(tc) = mean(cc(tc:tc+aveLen-1));
end
% Plot the results
subplot(2,2,4);
plot(tt,cc,at,ac,[min(tt) 0],cc1*[1 1],	...
	'--',[0 max(tt)],cc2*[1 1],'--','LineWidth',2);
grid on;
% A = axis();A(1) = min(tt)-TR;A(2) = max(tt)+TR;
axis(A); % Use the same x-axis range
xlabel('Time [sec]','FontSize',FSZ);
ylabel('Covariance','FontSize',FSZ);
legend('SWC','ASWC',['corr 1 =',num2str(cc1,1)],	...
	['corr 2 =',num2str(cc2,1)]);
title(['d) h=',num2str(winLen),'   g=',num2str(aveLen)],'FontSize',FSZ);
%% Save the figure
return
saveas(FH,'Figure4R1','png');
saveas(FH,'Figure4R1','fig');