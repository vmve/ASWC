% An average sliding window correlation method for dynamic functional 
% connectivity.
% Vergara VM, Abrol A, Calhoun VD
%
% Hum Brain Mapp. 2019 Jan 19. doi: 10.1002/hbm.24509
%
%
% ASWC: Simulate a wide spectrum signal with cosines
%
% Figure 6 and Figure 7
% Victor M. Vergara, PhD
% 2018-July-27
%
%
clear;clc;close all;
%% CONFIGURATION: Setup the simulation 
% Configuring variables
TR =  1;% [sec]
timeL = 300;% [sec]
nIter = 100;	% number of iterations
highestFreq = 0.1 ;	% in Hz
% variables calculated from configuration
nTR = ceil(timeL/TR);	% How many TR to simulate
t = TR*(0:nTR-1);	% time variable in seconds
Fstep = 1/(TR*nTR);
F  = Fstep:Fstep:highestFreq; % in Hz
nF = length(F);	% Number of Frequencies
%% CONFIGURATION: Frequency Spectrum Amplitudes 
% A = ones(1,nF);
% more BOLD-like
A = ones(1,nF) - (0:1:nF-1)/nF;
[n,m] = min(abs(F-0.025));
A(1:m) = A(m);
%% CONFIGURATION: Set the dynamic functional connectivity
dFC{1} = 0.9*tanh(ones(1,nTR));	% Constant (static connectivitY)
dFC{2} = 0.9*tanh(cos((  pi/nTR)*(1:1:nTR)));	% smooth transition
dFC{3} = 0.9*tanh(cos((2*pi/nTR)*(1:1:nTR)));	% a bit faster transition
Tcorr = 100;	% in seconds
dFC{4} = 0.9*tanh(cos(2*pi*t/Tcorr));	% corr frequency of 0.01 Hz
%%CONFIGURATION: Set the different window and ave sizes
hSWC = 100;	% window length in seconds
hASWC= [10,20,44,50,100];% ASWC window length in seconds
gASWC= 50*ones(1,length(hASWC));%floor(hASWC / 2); % ASWC averaging length in seconds
%% CONFIGURATION: Setup the Butterworth filter
fc = 0.01;	% cut-off 
[b,a] = butter(5,2*TR*fc,'high');

%% Calculate the variables for the subplot
nrows = length(dFC);
ncols = length(hASWC)+1;

%% Iterate and create the plots
for irow = 1:nrows
	for icol = 2:ncols
		% configure this iteration
		CC = dFC{irow};
		% Setup the ASWC method
		winlenASWCsec = hASWC(icol-1);	% sec:	 h= window length
		winlenASWC = ceil(winlenASWCsec/TR);
		avelenASWCsec = gASWC(icol-1);
		avelenASWC = floor(avelenASWCsec/TR);
		% Setup the SWC method
		winlenSWCsec = hSWC;	% sec:	 h= window length
		winlenSWC = ceil(winlenSWCsec/TR);
		%% Simulate
		MSE_SWC  = zeros(1,nIter);
		MSE_ASWC = zeros(1,nIter);

		for ii=1:nIter
			x = zeros(1,nTR);
			y =	zeros(1,nTR);
			for f=1:nF	% Add each frequency
				PP = 2*pi*rand(1,1);
				w = 2*pi*F(f);	% Frequency in rads/sec
				x = x + A(f)*cos(w*t + PP );
				y = y + A(f)*cos(w*t + PP + acos(CC) );
			end

			% Filter the data
			xf = filter(b,a,x);
			yf = filter(b,a,y);	

			yf = y;xf=x;
% 			%% plot the frequency spectrum to check
% 			figure(1);
% 			frq = 0:1:ceil(0.1/Fstep);
% 			Fori = abs(fft(x)).^2;
% 			Ffil = abs(fft(xf)).^2;
% 			subplot(2,1,1);
% 			plot(TR*(0:1:nTR-1),x,TR*(0:1:nTR-1),xf);
% 			legend('x','filtered x');
% 			xlabel('Time sec');
% 			subplot(2,1,2);
% 			plot(Fstep*frq,Fori(frq+1),Fstep*frq,Ffil(frq+1));
% 			xlabel('Frequency Hz');
% 			legend('Unfiltered','Filtered');

			%% Calculate the SWC
			SWCcorr = zeros(1,nTR-winlenSWC);
			for tc = 1:nTR-winlenSWC
				cr = corr(xf(tc:tc+winlenSWC-1)',yf(tc:tc+winlenSWC-1)');
				SWCcorr(tc) =atanh(cr);
			end
			SWCshift = ceil( winlenSWCsec / (2*TR) ) ;
			ccRange = SWCshift	+	(1:1:length(SWCcorr));
			MSE_SWC(1,ii) = mean((SWCcorr - atanh(CC(ccRange))).^2);
			%% Calculate the ASWC (Vergara)
			% Perform initial SWC
			ASWCcorr = zeros(1,nTR-winlenASWC);
			for tc = 1:nTR-winlenASWC
				 cr = corr(xf(tc:tc+winlenASWC-1)',yf(tc:tc+winlenASWC-1)');
				 ASWCcorr(tc) =atanh(cr);
			end
			% Perform Averging of SWC
			ASWCacorr = zeros(1, length(ASWCcorr) - avelenASWC);
			for tc = 1:length(ASWCcorr)-avelenASWC
				ASWCacorr(tc) = mean(ASWCcorr(tc:tc+avelenASWC));
			end
			ASWCshift = ceil( (winlenASWCsec + avelenASWCsec) / (2*TR) ) ;
			ccRange = ASWCshift	+	(1:1:length(ASWCacorr));
			MSE_ASWC(1,ii) = mean((ASWCacorr - atanh(CC(ccRange))).^2);
			fprintf('done with iter# %d out of %d\n',ii,nIter);

%			Just a Test. Comment if not testing.
% 			figure(2); 
% 			SWCshift = floor(winlenSWCsec /2);
% 			ASWCshift= floor((winlenASWCsec + avelenASWCsec)/2);
% 			plot(	t,atanh(CC),	...
% 					SWCshift + TR*(1:length(SWCcorr)),SWCcorr,	...
% 					ASWCshift+ TR*(1:length(ASWCacorr)),ASWCacorr)
% 			ylabel('correlation');
% 			xlabel('time [sec]');
% 			grid on;	
% 			legend('ground truth','SWC','ASWC')


		end
	%% Plot MSE subplots
	H = figure(1);
	subplot(nrows,ncols, ncols*(irow-1)+icol);
	scatter(MSE_SWC,MSE_ASWC,'b.',	...
		'LineWidth',2);
	hold on;
	xlabel('SWC MSE','FontSize',16);
	ylabel('ASWC MSE','FontSize',16);
% 	AX = axis();mx = max(abs(AX));mn=min(abs(AX));AX = [mn mx mn mx];
	mx = max([MSE_SWC,MSE_ASWC]);mn = min([MSE_SWC,MSE_ASWC]);
	axis([mn mx mn mx]);
	hold on;
	plot([mn mx],[mn mx]);
	hold off;
	title(['ASWC h=',num2str(winlenASWCsec), ...
		' g=',num2str(avelenASWCsec)],'FontSize',16);
	
	% plot the dFC
	H = subplot(nrows,ncols, ncols*(irow-1)+1);
	plot(	t,atanh(CC),'LineWidth',2);
	ylabel('atanh(correlation)','FontSize',16);
	xlabel('time [sec]','FontSize',16);
	grid on;
	axis([0 max(t) -1 1]);
	title('Ground Truth','FontSize',16);
	%% Plot Tracking  subplots
	H = figure(2);
	subplot(nrows,ncols, ncols*(irow-1)+icol);
	SWCshift = floor(winlenSWCsec /2);
	ASWCshift= floor((winlenASWCsec + avelenASWCsec)/2);
	plot(	t,atanh(CC),	...
			SWCshift + TR*(1:length(SWCcorr)),SWCcorr,'.-',	...
			ASWCshift+ TR*(1:length(ASWCacorr)),ASWCacorr,'k--','LineWidth',2);
	title(['ASWC h=',num2str(winlenASWCsec), ...
		' g=',num2str(avelenASWCsec)],'FontSize',16);
	if(irow > 1);
		axis([0 (max(t)+1) -1.5 1.5]);
	end
	% plot the dFC
	H = subplot(nrows,ncols, ncols*(irow-1)+1);
	plot(	t,atanh(CC),'LineWidth',2);
	ylabel('atanh(correlation)','FontSize',16);
	xlabel('time [sec]','FontSize',16);
	axis([0 max(t) -1 1]);
	grid on;	
	title('Ground Truth','FontSize',16);

	
	
	end

end
% plot(CC,ASWCartstd,'o',CC,SWCartstd,'x')
% plot(CC,ASWCmean,'o',CC,SWCmean,'x')
return
%% Plot each CC
	winlenASWCsec = hASWC(icol-1);	% sec:	 h= window length
		winlenASWC = ceil(winlenASWCsec/TR);
		avelenASWCsec = gASWC(icol-1);
		avelenASWC = floor(avelenASWCsec/TR);
		% Setup the SWC method
		winlenSWCsec = hSWC;	% sec:	 h= window length
		winlenSWC = ceil(winlenSWCsec/TR);
return
