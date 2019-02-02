% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% Average Sliding Window Correlation (ASWC)
% Figure 2
%
% Victor M. Vergara, PhD
% 2018-March-01
%
% Modified:
% 2018-April-16	: updated the cuttoff and filter equations
%	
%
clear;close all;clc;
%% Setup the simulation configuration
TR =   0.1;	%	TR [sec]
timeL = 880;%	Time Length of the Cosines [sec]
sCorr = 1.0;%	Simulated correlation

F = 1/1000 : 1/1000 : 1/20;		% Frequency range in Hz
A=sqrt(2);		% Cosine Amplitude
P=acos(sCorr);	% Phase between Cosines

winSizeSec = 44;	% Window Length [sec]
winSizeTR = ceil(winSizeSec/TR);
winSizeSecRange = 2:100;	% SEcond subplot

aveSizeSecRange = 2:40;	% SEcond subplot 

FSZ = 12; %Figure font size

ButterWF = [0.01 0.025];% Hz (add these Botterwort responses (compare)
FilterOrd = 2;	% filter order = 2

%% Calculate some setup values based on the setup

nTR = ceil(timeL/TR);	% How many TRs in the simulation (time courses)

% Obtain the cutoff frequency
a = 1/2;
r = 1 - sqrt(1-a);
frqHfac = (1/pi)*sqrt( 10*(1 - sqrt (1-(6/5)*r)));
frq = (frqHfac/winSizeSec);
%frq = 0.4441/winSizeSec; % only for a = 1/2;

%% Display the Simulation Configuration
fprintf('Sliding Window Simulation (Cosine Model)\n');
fprintf('----------------------------------------\n');
fprintf('Configuration:\n');
fprintf('\n\t     TR [sec] = %.1f\n',TR);
fprintf('\tNumber of TRs = %d\n',nTR);
fprintf('\tDuration[sec] = %.1f\n',TR*nTR);
fprintf('\twin Len [sec] = %d\n',winSizeSec);
fprintf('\nEstimated Simulation :\n');
fprintf('\n\tcutoff  [Hz ] = %.3f\n',frq);
fprintf('\n');

%% Perform the Simulation
asympCov = zeros(1,length(F));
t = TR*(0:nTR-1);	% time variable
%% We are going to get the asymptotic behavior for different frequencies
for f = 1:length(F)
	%% Get the frequency to use in this iteration  [ rads/sec ]
	w = 2*pi*F(f);
	%% Create the cosines to be simulated
	x = A*cos(w*t);
	y = A*cos(w*t	+	P);
	clear cc;
	for tc = 1:nTR-winSizeTR
		cm = cov(x(tc:tc+winSizeTR-1)',y(tc:tc+winSizeTR-1)');
		cc(tc) = cm(1,2);
	end
	asympCov(f) = mean(cc);
end

%% get the equation
E = sCorr*(1-(sinc(F.*winSizeSec)).^2);

%% -------  Plot Simulation and Equation -----
% figure('rend','painters','pos',[10 100 1200 400]);
FH = figure('rend','painters','pos',[10 100 900 600]);
%% First Plot x = F*h
% this first plot shows the normalized version for x=F*h
x = 0:0.01:2;
y = 1 - sinc(x).^2 ;
subplot(2,2,1);
plot(x,y,'LineWidth',2);
xlabel('Frequency (f) X Window Length (h)','FontSize',FSZ);
ylabel('Normalized Covariance','FontSize',FSZ);
title('a) Normalized Asymptotic Plot x=h x f','FontSize',FSZ);
grid on;
%% Second plot: Window Length by frequency in the range of interested
for kk=1:length(winSizeSecRange)
	for jj=1:length(F)
		C(kk,jj) = sCorr*(1-(sinc(F(jj).*winSizeSecRange(kk))).^2);
	end
end
subplot(2,2,2);
ax = gca;
imagesc(C);

colorbar;
% Set the x axis
xticks = 1:10:length(F);  %adjust as appropriate, positive integers only
if(xticks(end) ~= F(end));xticks=[xticks,length(F)];end;
for kk=1:length(xticks);xlabels{kk} = num2str( F(xticks(kk)),'%.2f' );end
  set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xlabel('Frequency [Hz]','FontSize',FSZ);
% Set the y axis
yticks = (10:10:length(winSizeSecRange))-1;  %adjust as appropriate, positive integers only
if(yticks(end)~=winSizeSecRange(end));yticks=[yticks,length(winSizeSecRange)];end;
for kk=1:length(yticks);
	ylabels{kk} = num2str( winSizeSecRange(yticks(kk)),'%d' );
end
set(ax, 'YTick', yticks, 'YTickLabel', ylabels);
ylabel('Window Length h [sec]','FontSize',FSZ);
set(gca,'fontsize',10)
colormap(hot);
% Set grid lined
ax.GridColor = [0.0 0.0 1.0];
grid on;
title('b) Normalized Asymptotic h vs. f','FontSize',FSZ);

%% Third Plot For the window length of interest 
subplot(2,2,3);
plot(frq,0.5*sCorr,'ko',F,asympCov,'b','LineWidth',2,'MarkerSize',10);
LegendLabel = {'Cutoff (C_x_y/2)','Averaged C_x_y'};
%% create Butterworth responses
mark = {'--','-.',':'};
ColOrd = get(gca,'ColorOrder');
for kk=1:length(ButterWF)
	if(kk>3) continue;end;	% max of 2 filter plots
	hold on;
	cutoff = ButterWF(kk);
	BF = ((F./cutoff).^(2*FilterOrd))./(1+(F./cutoff).^(2*FilterOrd));
	plot(F,BF,mark{kk},'Color',ColOrd(kk+1,:),'LineWidth',2);
	LegendLabel{kk+2} = ['Butterworth Cutoff ',num2str(ButterWF(kk))];
end
hold off
set(gca,'fontsize',FSZ)
lgd = legend(LegendLabel,'Location','SouthEast');
lgd.FontSize = 8;
xlabel('Frequency [Hz]','FontSize',FSZ);
ylabel('Covariance','FontSize',FSZ);
title(['c) Asymptotic Averaging (C_x_y=1) for h=',	...
	num2str(winSizeSec),' sec'],'FontSize',FSZ);
grid on;
A = axis();A(4)=1.1;axis(A);

%% Fourth Plot: for h= winSizeSec : averaging and frequency 
clear mse;
mse = zeros(	length(aveSizeSecRange),	length(F)	);
for kk=1:length(aveSizeSecRange)
	for jj=1:length(F)
		f = F(jj);
		h = winSizeSec;
		h2 = h/2;	% in the equations delta = h/2
		g = aveSizeSecRange(kk);
		g2 = g/2;	% in the equations nabla = g/2
		w = 2*pi*f;
		
		K = (1-(sinc(f*h).^2));
		G = (2/(g*h*(w^2))) * sin(2*w*g2*TR) * sin(2*w*h2*TR) * ...
			(	cos(w*h2*TR) - (sin(w*h2*TR)/(h*w))	);
		
		Cxy = sCorr*K + cos(2*w*TR*t + acos(sCorr))*G;
		
		mse(kk,jj) = mean((sCorr - Cxy).^2); %mean square error
	end
end
subplot(2,2,4);
ax = gca;
imagesc(sqrt(mse),[0 0.5]);
h = colorbar;
set(get(h,'title'),'string','RMSE');
% Set the x axis
xticks = 1:10:length(F);  %adjust as appropriate, positive integers only
if(xticks(end) ~= F(end));xticks=[xticks,length(F)];end;
for kk=1:length(xticks);xlabels{kk} = num2str( F(xticks(kk)),'%.2f' );end
  set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xlabel('Frequency [Hz]','FontSize',FSZ);
% Set the y axis
yticks = (10:10:length(aveSizeSecRange))-1;  %adjust as appropriate, positive integers only
if(yticks(end)~=aveSizeSecRange(end));yticks=[yticks,length(aveSizeSecRange)];end;
for kk=1:length(yticks);
	ylabels{kk} = num2str( aveSizeSecRange(yticks(kk)),'%d' );
end
set(ax, 'YTick', yticks, 'YTickLabel', ylabels);
ylabel('Averaging Length g [sec]','FontSize',FSZ);
set(gca,'fontsize',10)
title('d) g vs. f Response for C_x_y=1','FontSize',FSZ);
colormap(hot);
% Set grid lined
ax.GridColor = [0.0 0.0 1.0];
grid on;
%% Save the figure
return
saveas(FH,'Figure2R1','png');
saveas(FH,'Figure2R1','fig');