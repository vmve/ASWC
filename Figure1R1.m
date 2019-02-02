% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% Average Sliding Window Correlation (ASWC)
% Figure 1
%
% Victor M. Vergara, PhD
% Created on 2018-February-26
%
% This code does 2 things:
%	a) Replicates the cosine-covariance model for 
%	b) Show the new estimate for the average of the model
%
% Modified on 2018-10-22
%
%  - Changes on the title to ASWC.
%  - The first lobe occurs at h = 48 and Figure 1c reflect this change.
%    Previously it was set to 50, but this is incorrect.
%
clear;close all;clc;
%% Setup the simulation configuration
TR =  1;	% [sec]
timeL = 120;% [sec]
sCorr = 0.2;% Simulated correlation

F = 1/40;		% Cosine Frequency
A=sqrt(2);		% Cosine Amplitude

minW = 1;		% Minimum Window Size [sec]
maxW = 120;		% Maximum Window Size [sec]

plotWLength = 48;% Selection for Figure 1c 

FSZ = 14; % Figure font size

% The next values are calculated from configured ones
nTR = ceil(timeL/TR);	% How many points
P=acos(sCorr);	% Phase between Cosines

%% Check failsafe
if maxW > ceil(nTR*TR);
	maxW = ceil(nTR*TR);
	warning('Cannot Simulate Window Sizes larger than Number of Points');
end;
%% Create the cosines to be simulated
w = 2*pi*F;	% Frequency in rads/sec
t = TR*(0:nTR-1);	% time variable
x = A*cos(w*t);
y = A*cos(w*t	+	P);

%% Display the Simulation Configuration
fprintf('Sliding Window Simulation (Cosine Model)\n');
fprintf('----------------------------------------\n');
fprintf('\n\t     TR [sec] = %.1f\n',TR);
fprintf('\tNumber of TRs = %d\n',nTR);
fprintf('\tDuration[sec] = %.1f\n',TR*nTR);
fprintf('\n');
fprintf('Parameters for the 2 Simulated Cosines\n');
fprintf('\t [Hz] \t Amp \tPhase [rad]\n');
for kk=1:length(F)
	fprintf('\t%.3f\t%.3f\t%.3f\n',F(kk),A(kk),P(kk));
end
fprintf('\ncorrelation :: cos(%.3f) = %.3f\n',P,corr(x',y'));
fprintf('\nwindow sizes from %d to %d\n',minW,maxW);

%% Perform the Simulation for Window Length
for wSzSec = minW:maxW
	wSz = ceil(wSzSec/TR) - 1;
	clear cc;
	for tc = 1:nTR-wSz
		cm = cov(x(tc:tc+wSz)',y(tc:tc+wSz)');
		cc(tc) = cm(1,2);
	end
	if ~exist('cc');continue;end
	windowShifts{wSzSec}  = (1:1:tc)*TR;
	covEstimation{wSzSec} = cc;
	windowLengths(wSzSec) = wSzSec;
end
%% Calculate the average over all the covariances
meanEstimation = zeros(1,length(covEstimation));
for kk = 1:length(covEstimation)
	meanEstimation(kk) = mean(covEstimation{kk});
end
%% Calculate the theoretical average (Equation)
D = ((windowLengths/TR))/2;
eqEstimation = sCorr*(1 - (sinc(F.*windowLengths)).^2);

%% Simulate for different averaging lengths
% find the data first (closest to the simulation we obtained)
[minV,minIndex] = min(abs(windowLengths - plotWLength));
wLenCovData = covEstimation{minIndex};
timeCovData = windowShifts {minIndex};
nWin = length(wLenCovData);
for aveSz = 1:nWin
	clear aveW;
	for tc = 1:nWin - aveSz
		aveW(tc) = mean(wLenCovData(tc:tc+aveSz-1));
	end
	if ~exist('aveW');continue;end
	aveShifts{aveSz}  = (1:1:tc)*TR;
	aveEstimation{aveSz} = aveW;
	aveLengths(aveSz) = aveSz*TR;
end
%% Obtain data using the Empirical (need ot analyze the equations) avelen
g = round(1/(2*F));
aveLen = ceil(g/TR);
nWL = length(windowLengths);
for kk=1:nWL
	cc = covEstimation{kk};
	if (aveLen > length(cc));continue;end;
	aveSWCval(kk) = mean(cc(1:aveLen));
	aveSWClen(kk) = windowLengths(kk);
end

%% Create the figure, first by opening the figure handle
FH = figure('rend','painters','pos',[10 100 900 600]);

%% Plot The Cosines
subplot(2,2,1);
plot(t,x,'-o',t,y,'--x','LineWidth',1);legend('x','y');
xlabel('Time [sec]','FontSize',FSZ)
ylabel('Cosine Amplitude','FontSize',FSZ);
title(['a) Cosines of Frequency = ',num2str(F),' Hz'],'FontSize',FSZ);
%% Plot the Simulation results for Window Length
subplot(2,2,2);
% plot the theoretical configured correlation
plot([0 max(windowLengths)],sCorr*[ 1 1],'b--','LineWidth',2);
hold on;grid on;
% plot each of the estimated covariances
for kk = 1:length(covEstimation)
	cc = covEstimation{kk};
	wL = windowLengths(kk);
	plot(wL*ones(1,length(cc)),cc,'.k','MarkerSize',1.5);
end
xlabel('Window Length h [sec]','FontSize',FSZ)
ylabel('Covariance','FontSize',FSZ)
hold off;
if (sCorr>0)
	legend(	'Ground Truth Correlation','Shifted SWCs',	...
		'Location','SouthEast');
else
	legend(	'Ground Truth Correlation','Shifted SWCs',	...
		'Location','NorthEast');
end
title('b) All Possible Shifts of SWCs','FontSize',FSZ);
% save the axis for the other plots
AX = axis();
%% Plot the Simulation results for Average Length
subplot(2,2,3);
% plot the theoretical configured correlation
plot([0 max(windowLengths)],sCorr*[ 1 1],'b--','LineWidth',2);
hold on;grid on;
% plot each of the estimated covariances
for kk = 1:length(aveEstimation)
	cc = aveEstimation{kk};
	aL = aveLengths(kk);
	plot(aL*ones(1,length(cc)),cc,'.k','MarkerSize',1.5);
end
xlabel('Averaging Length g [sec]','FontSize',FSZ)
ylabel('Covariance','FontSize',FSZ)
hold off;
axis(AX);
if (sCorr>0)
	legend(	'Ground Truth Correlation','Shifted Averages',	...
		'Location','SouthEast');
else
	legend(	'Ground Truth Correlation','Shifted Averages',	...
		'Location','NorthEast');
end
title(['c) ASWC for h= ',num2str(plotWLength),' sec'],'FontSize',FSZ);
%% Plot the Simulation results for Average Length
subplot(2,2,4);
% plot the theoretical configured correlation
plot([0 max(windowLengths)],sCorr*[ 1 1],'b--','LineWidth',2);
hold on;grid on;
% plot the theoretical average of the correlation
plot(windowLengths,eqEstimation,'r.-','LineWidth',2);
% plot the simulated mean of the sliding windows
plot(aveSWClen,aveSWCval,'k','LineWidth',2);
% plot each of the estimated covariances
xlabel('Window Length h [sec]','FontSize',FSZ)
ylabel('Covariance','FontSize',FSZ)
hold off;
axis(AX);
if (sCorr>0)
	legend(	'Ground Truth Correlation','Equation (5)',	...
		'Average of SWCs','Location','SouthEast');
else
	legend(	'Ground Truth Correlation','Equation (5)',	...
		'Average of SWCs','Location','NorthEast');
end	
title(['d) ASWC (any shift) g= ',num2str(g),' sec'],'FontSize',FSZ);
%% Save the figure
return
saveas(FH,'Figure1R1','png');
saveas(FH,'Figure1R1','fig');