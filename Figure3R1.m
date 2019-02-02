% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% Average Sliding Window Correlation (ASWC)
% Figure 3
%
% Victor M. Vergara, PhD
% 2018-April-26
%
%
clear;close all;clc;
%% Setup the simulation configuration
TR =  0.1;% [sec]
timeL = 300;% [sec]

F  = 0.001:0.001:0.1;%:0.01:0.1;%:0.01:0.1;		% Cosine Frequencies
nF = length(F);
A = sqrt(2);	% Cosine Amplitudes sqrt(2/N of cosines)
PhaseM = 0;
% we will randomize the phase...
%P=acos(sCorr)*[0.5 1];	% Phase between Cosines

wSzSecASWC = 44;	% sec:	 h= window length
wSzASWC = ceil(wSzSecASWC/TR);
gSzSecASWC = 50;
gSzASWC = floor(gSzSecASWC/TR);

wSzSecSWC = 100;	% sec:	 h= window length
wSzSWC = ceil(wSzSecSWC/TR);

nTR = ceil(timeL/TR);	% How many points

%% Iterate as many times as requested
t = TR*(0:nTR-1);	% time variable
FSZ = 16;
%% Create the cosines to be simulated

fc = 0.01;
[b,a] = butter(5,2*TR*fc,'high');

x = zeros(1,nTR);
y = zeros(1,nTR);
for f=1:nF	% Add each frequency
	w = 2*pi*F(f);	% Frequency in rads/sec
	x = A*cos(w*t	);
	y = A*cos(w*t+	PhaseM);
	
	% Filter the data
	xf = filter(b,a,x);
	yf = filter(b,a,y);	
	
	% Calculate the Leonardi
	cc = zeros(1,nTR-wSzSWC);
	cr = zeros(1,nTR-wSzSWC);
	for tc = 1:nTR-wSzSWC
 		cm = cov(xf(tc:tc+wSzSWC-1)',yf(tc:tc+wSzSWC-1)');
 		cc(tc) = cm(1,2);
	end
	L{f} = cc;
	
	% Calculate the Averaged SWC (Vergara)
	cc = zeros(1,nTR-wSzASWC);
	for tc = 1:nTR-wSzASWC
 		cm = cov(x(tc:tc+wSzASWC-1)',y(tc:tc+wSzASWC-1)');
 		cc(tc) = cm(1,2);
	end
	nCC = length(cc);
	ac = zeros(1,nCC - gSzASWC);
	for tc = 1:nCC-gSzASWC
		ac(tc) = mean(cc(tc:tc+gSzASWC));
	end
	V{f} = ac;

end

FH = figure('rend','painters','pos',[10 100 1100 300]);
stdSWCerror = zeros(1,nF);
stdASWCerror = zeros(1,nF);
for f=1:nF	% Add each frequency

	subplot(1,3,1);
	cc = L{f};
	nc = length(cc);
	plot(F(f)*ones(1,nc),cc,'k.','MarkerSize',3);
	hold on;
	stdSWCerror(f) = std(cc);
	
	subplot(1,3,2);
	cc = V{f};
	nc = length(cc);
	plot(F(f)*ones(1,nc),cc,'k.','MarkerSize',3);
	hold on;
	stdASWCerror(f) = std(cc);

end

%fix the axis
subplot(1,3,1);Ax = axis();grid on;
xlabel('Frequency Hz','FontSize',FSZ);
ylabel('Covariance','FontSize',FSZ);
title(['a) SWC with h = ',num2str(wSzSecSWC),' sec'],'FontSize',FSZ);

subplot(1,3,2);axis(Ax);grid on;
xlabel('Frequency Hz','FontSize',FSZ);
ylabel('Covariance','FontSize',FSZ);
title(['b) ASWC with h = ',num2str(wSzSecASWC),	...
	' sec, g = ',num2str(num2str(gSzSecASWC)),' sec'],'FontSize',FSZ);

subplot(1,3,3);
plot(F,log2(stdSWCerror./stdASWCerror));
grid on;
xlabel('Frequency Hz','FontSize',FSZ);
ylabel('log_2( \sigma_{SWC} / \sigma_{ASWC} )','FontSize',FSZ);
title('c) Artifact Comparison','FontSize',FSZ);


%% Save the Figures
saveas(FH,'Figure3R1','fig');
saveas(FH,'Figure3R1','png');


