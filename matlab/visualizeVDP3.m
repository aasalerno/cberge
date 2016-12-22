% ### visualizeVDP3.m ###
% analyze/visualize output from vdModelP#.m
% [original A. Salerno; re-hacked by C. Bergevin]

% [02.11.16] - option to analyze (summed) waveform via plotSOAEwf10.m (Figs.xxx)
% [02.12.16] - new bit to allow for better tiling when there are multiple figs.;
% also incorporated use of plotSOAEwf10.m
% [02.29.16] - added functionality (via "fig.10") to plot heat map of
% steady-state for all oscillators (stems from ASA 2015 slides); also
% incorporating the elements of ./LineVideo.m here so to animate (Fig.111) a line
% heat-map so to visually see the dynamics

clear
% ----------------
file= 'test.mat';
oscR= [30 33];         % specify oscillators to examine analytic signal (Figs.3)
oscRB= 1;               % specify oscillator to examine waveform (real(Z)) and s.s. spectrum (Fig.4)
oscRC= [38:41];      % specify oscillators to compare s.s. spectrum (Fig.5)
lengthS= 8192*2;      % length of buffer (at end of total waveform) when summing oscillators to calculate spectrum [8192*2]

% filter properties for a "ring of fire"-type plot
In.CF= 2.64;
In.BW= 0.1;

% ----------------
baseCoords= [720 350];  % base position for Fig.1
Wsize= [560 520];        % window size; [560 520] is default
offW= 20;
fNUM= 0;    % dummy indexer for fig #
In.nPts= 128;   % waveform length of "ring of fire" fig.

load(file);     % load in data
SR= 1/data.info.params.stepSize;
freq= [0:lengthS/2];    % create a freq. array (for FFT bin labeling)
freq= SR*freq./lengthS;

% =====================================================================
%set(gca,'DefaultFigureColormap',jet)

% =====================================================================
% Figure 1 [vdSpec(data)]
% sum waveforms for displacement (real part of Z) for all oscillators, then
% compute FFT and plot spectrum
if 1==1
    %data.X = real(data.Z);
    %data.wf = sum(data.X(end-lengthS+1:end,:),2)';     % sum rows for 'steady-state' section
    %data.wf = sum(data.X(end-lengthS+1:end,oscR(1)),2)';  % sum only a single oscillator (same as Fig.4)
    %data.FFT = rfft(data.wf);                                % calculate the FFT
    %nPts = length(data.wf);
    data.FFT = rfft(sum(real(data.Z(end-lengthS+1:end,:)),2)'); % do all the above in a single line
    X1 = db(abs(data.FFT));
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    plot(freq,X1,'k'); hold on; grid on;
    xlim([0 2*max(data.info.params.w/(2*pi))]); xlabel('Frequency'); ylabel('Magnitude (dB)');
    title('Spectrum of summed waveform of all oscillators')
end


% =====================================================================
% Figure 2   [vdAnalyze(data)]
% takes each of the individual waveforms from the 'steady-state' (re
% lengthS), then uses two methods to determine the 'dominant' frequency:
% Method 1 - simply grab largest peak from spectrum
% Method 2 - calculate the average (unwrapped) phase velocity
if 1==1
    N= size(data.Z,2);  % # of oscillators (needed for indexing)
    X = real(data.Z(end-lengthS+1:end,:));  % extract 'steady-state' at the end re lengthS, only with the real part
    xOscF = zeros(1,N);     % method 1 storage array
    xArgF = zeros(1,N);     % method 2 storage array (argument for average slope)
    tDom = data.info.params.T(end-lengthS+1:end); % timevalues for extracted s.s. (used for polyfit)
    % ---
    findloc = 1;    % indexer for while loop
    for k = 1:N
        xFFT = rfft(X(:,k)); %The RFFT of the dataset,
        % ---
        % Method 1 - Most Prominent Freqency using FFT (we don't care about the value)
        [~,locs] = findpeaks(abs(xFFT),'SORTSTR','descend'); % MAXIMUM peak comes first, because of 'descend'
        if findloc == 1 && ~isempty(locs)
            xOscF(k) = freq(locs(1)); % If locs found something append it to the list
        else
            findloc = 0; %if not, this prevents any overflow errors from occuring
        end
        % ---
        % Method 2 - steady-state phase slope [CB: not presently sure how legit this 'method' is, but interesting idea...]
        phs = unwrap(angle(data.Z(end-lengthS+1:end,k)))/(2*pi); %Unwraps the arg information, diving by 2*pi b/c of omega
        m = polyfit(tDom,phs',1); % Find the best fit
        xArgF(k) = m(1); %m(1) is m from y = m*x+b
    end
    % ---
    %figure(2); clf;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    plot(xOscF(2:end),'-o','MarkerSize',5) %Data from method 1 being plotted in blue
    hold on; grid on;
    plot(xArgF(2:end),'-+','MarkerSize',5,'Color','r');    plot(data.info.params.w/(2*pi),'k--')
    %plot(mean(abs(real(data.Z)).^2).^.5,'-.k');
    plot(mean(abs((data.Z(:,2:end)))),'k-');
    xlabel('k - Oscillator Number'); ylabel('Strongest Oscillatory Frequency')
    legend('max. spectral peak','phase velocity','CF distribution','mean(abs(Z))','Location','NorthWest')
    plot([oscR(1):oscR(2)],xOscF([oscR(1):oscR(2)]),'gs','LineWidth',2); % also indicate osc. for Fig.3
    title('Frequency plateaus a la Vilfan & Duke (2008) analysis');
    %set(2,'OuterPosition',[baseCoords(1)+offW baseCoords(2)-offW Wsize]);
end
%

% =====================================================================
% Figure 3   [plotHilb(data,[oscMin oscMax])]
% For group of oscillators specified in oscR, calculate the analytic signal
% of the entire waveform and plot associated envelope and inst. phase (the
% latter are referenced to the first oscillator)
if 1==0
    intStp = 1;
    data.X = real(data.Z);
    oscN= oscR;
    if intStp==1;
        oscN = oscN(1):oscN(end);
    end
    osc1 = round(oscN(1));
    for j = 1:numel(oscN)
        data.Hilb(:,j) = hilbert(data.X(:,osc1+j-1));
        data.HMag(:,j) = abs(data.Hilb(:,j));
        %data.HArg(:,j) = unwrap(angle(data.Hilb(:,j)))/(2*pi);
        data.HArg(:,j) = angle(data.Hilb(:,j));
        oscNstr{j} = ['Oscillator ' num2str(oscN(j))];
    end
    HArgRef = data.HArg(:,1);   % set 1st oscillator as phase ref.
    % ---
    % "correct" phase re to the 1st oscillator
    for j = 1:numel(oscN)
        %data.HArg(:,j) = data.HArg(:,j) - HArgRef - fix(data.HArg(:,j) - HArgRef);
        data.HArg(:,j) = data.HArg(:,j) - HArgRef;
    end
    % ---
    %figure(3); clf;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    for nn=1:size(data.HMag,2)
        shd=0.7*(nn-1)/size(data.HMag,2);
        subplot(211); plot(data.info.params.T,data.HMag(:,nn),'Color',[0 0 0]+shd); hold on;
        %subplot(212); plot(data.T,data.HArg(:,nn),'Color',[0 0 0]+shd); hold on;
        subplot(212); plot(data.info.params.T,unwrap(data.HArg(:,nn))/(2*pi),'Color',[0 0 0]+shd); hold on;
    end
    subplot(211); xlabel('Time'); ylabel('Envelope'); grid on
    subplot(212); xlabel('Time'); ylabel('(relative) Phase [cycles]'); grid on
end


% =====================================================================
% Figure 4
% Plot the entire waveform (real(Z)) for the first oscillator specified in
% oscR, then also plot the spectrum for the 'steady-state' (lengths re the end) portion
if 1==0
    %figure(4); clf;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    % plot waveform and spectrum (of last segment portion)
    wf= real(data.Z(:,oscRB));
    wfSS= wf(end-lengthS+1:end);
    % ---
    subplot(211);
    plot(data.info.params.T,wf); hold on; grid on;  plot(data.info.params.T(end-lengthS+1:end),wfSS,'r.')
    xlabel('Time'); ylabel('Position');
    title(['Entire waveform for oscillator #',num2str(oscRB)])
    % ---
    subplot(212);
    %plot(db(rfft(real(data.Z(:,oscR((1)))))))
    plot(freq,db(rfft(wfSS))); xlim([0 10]); hold on; grid on;
    xlabel('Frequency'); ylabel('Magnitude [dB]');
    title(['Steady-state spectrum for oscillator #',num2str(oscRB),' (red waveform portion)'])
    % ---
    % if Fig.3 is plotted, can also include the analytic signal here as a
    % reality check (asumes oscRB= oscR(1)); looks good/consistent
    if (1==0),  subplot(211); plot(data.T,data.HMag(:,1),'g-'); end
end

% =====================================================================
% Figure 5
% compare steady-state spectra
if 1==0
    %figure(5); clf;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    for nn=1:numel(oscRC)
        shd=0.7*(nn-1)/numel(oscRC);    % shading factor for visualizing
        wf= real(data.Z(:,oscRC(nn)));
        wfSS= wf(end-lengthS+1:end);
        plot(freq,db(rfft(wfSS)),'Color',[0 0 0]+shd); hold on;
    end
    xlim([0 10]); grid on; xlabel('Frequency'); ylabel('Magnitude [dB]');
    title(['Steady-state spectrum for oscillator #s',num2str(oscRC)])
end


% ====================================================================
% Figure 6
% True "clustering"
if 1==1
    %figure(6); clf
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    xFFT = rfft(X);
    imagesc((db(abs(xFFT))),'YData',freq); colormap default;
    set(gca,'YDir','normal')
    ylim([0 6]); h= colorbar; caxis([-100 0]); ylabel(h, 'dB')
    xlabel('Oscillator Number (k)'); ylabel('Frequency (kHz)')
end


% ====================================================================
% Figure 7
% Spectrogram of a single oscillator (specified via oscRB)
if 1 == 0
    %figure(7); clf;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    xSpec = real(data.Z(end-lengthS+1:end,oscRB));
    [S,~,T,P] = spectrogram(xSpec,1024,ceil(0.95*1024));
    F = linspace(freq(1),freq(end),size(S,1));
    surf(F,T,10*log10(abs(P))','EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90); colorbar;
    ylabel('Time (s)');
    xlabel('Frequency (kHz)');
    title(['Steady State Spectrogram of Oscillator ' num2str(oscRB)])
    xlim([data.info.params.frange(1)-1 data.info.params.frange(2)+1])
end


% ====================================================================
% Figure 8
% Instantaneous freq. -- check for "bobble"
if 1==0
    %figure(8); clf; hold on;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    oscs = [40 41];
    xHilb = hilbert(real(data.Z(end-8192+1:end,:)));
    xEnv = abs(xHilb);
    xPh = unwrap(angle(xHilb));
    t = data.info.params.T(end-8192+1:end);
    instfreq = SR/(2*pi)*diff(xPh(:,oscs(1)));
    plot(t(1:end-1),instfreq,'k'); hold on;
    plot(t,5 + X(end-8192+1:end,oscs(1)),'b')
    plot(t,0 + X(end-8192+1:end,oscs(2)),'r')
    legend(['f_{inst}^{' num2str(oscs(1)) '}'],['Osc ' num2str(oscs(1))],['Osc ' num2str(oscs(2))])
    ylim([-1.5 6.5])
    ylabel('Position (arb) & f_{inst} [kHz]')
    xlabel('Time')
end

% =====================================================================
% Figure 9
% phase coherence plot for oscillators specified in oscR (not oscRC!)
% (see Bergevin et al. PNAS 2015 eqn.S1)
% (requires Fig.3 to be made such that data.Hilb is generated)
if 1==0
    
    PC= sum(data.Hilb,2)./sum(abs(data.Hilb),2);
    
    % ---
    %figure(9); clf;
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    subplot(211); plot(data.info.params.T,abs(PC));
    subplot(212); plot(data.info.params.T,angle(PC)/(2*pi));
    subplot(211); ylabel('Envelope'); title('Phase coherence'); grid on
    subplot(212); xlabel('Time'); ylabel('Phase [cycles]'); grid on
    
end


% =====================================================================
% Figure 10
% "heat map" plot of the (steady-state) time waveform for all
% oscillators (see Anthony's 2015 ASA slide 16)
if 1==1
    ssAll= data.Z(end-lengthS+1:end,:);    % extract complex steady-state response for all osc.
    ssAllS= abs(ssAll).*sin(angle(ssAll));  % make a "signed" version (i.e., real-valued, but phase incl.)
    tSS= data.info.params.T(end-lengthS+1:end);  % time array for SS portion
    ossNum= linspace(1,size(ssAll,2),size(ssAll,2));
    % ====
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    f10= imagesc('XData',tSS,'YData',ossNum,'CData',ssAllS'); colormap jet; c10= colorbar;
    xlabel('Time (steady-state)'); ylabel('Osc. #'); ylabel(c10,'Amplitude'); ylim([0 size(ssAll,2)])
end


% =====================================================================
% Figure 11
% animated "line video" of the steady-state portion
if 1==0
    SR= 1/data.info.params.stepSize;
    chunk= real(data.Z(end-lengthS+1:end,:));
    chunk= rot90(chunk); % rotate 90 to get proper orientation
    ossNum= linspace(1,size(ssAll,2),size(ssAll,2));
    chunk= flipud(chunk);
    
    % ====
    fNUM= fNUM+1; figure(fNUM); clf;    % next two lines to handle offset tiling
    set(fNUM,'OuterPosition',[baseCoords(1)+offW*(fNUM-1) baseCoords(2)-offW*(fNUM-1) Wsize]);
    colormap('bone')
    for n=1:size(chunk,2)
        % as a 1-D-ish color heatmap
        imagesc('YData',ossNum,'CData',chunk(:,n)); Z= colorbar; caxis([-1 1]); ylim([0 size(ssAll,2)]);
        ylabel('Oscillator #'); ylabel(Z,'Inst. Amplitude');
        pause(0.05);
        
        % as a more simple scatter plot
        %plot(ossNum,chunk(:,n),'ko-'); grid on; axis([0 size(ssAll,2) -1.2 1.2])
        %pause(0.05);
    end
end

% ====================================================================
% % Figure xxff
% "ring of fire" type fig. via plotSOAEwf10.m
if 1==0
    In.waveform= sum(real(data.Z(end-lengthS+1:end,2:end)),2)';     % sum all oscillators (sans papilla)
    
    % -----
    In.SR= SR; In.plotScaleD= 1;
    In.arAR= 0; % no need for AR here....
    In.plotFMAX= 10;  % changes to lower freq. scale due to "SR"
    In.AnalyticHist= 1;
    In.AmplDistSpec= 1;
    In.AmplDistSpec= 1;
    In.xxOFF= 300;
    In.xxHistL= 50;
    In.FiltSpecBoth= 0;
    % -----
    Data= plotSOAEwf10(In); view([0 90]);
end