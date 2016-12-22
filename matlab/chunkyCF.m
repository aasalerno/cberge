% ### chunkyCF.m ###    03.11.16 CB

% fiddling w/ a front end re tonotopic map assumption for the lizard model
% (i.e., developing option for P.CFtype in vdModelP#.m)
% ---> basically the oscillators form groups w/ roughly the smae CF, as
% motivated by the notion that bundles are morphologically (and
% physiologically?) similar within a given longitudinal cross-section

clear

% ---------------------------------------
P.N = 140; % number of active oscillators
P.CFtype= 1;  % tonotopic map type: 0-linear, 1-log (determines w)
P.frange= [1 4.5];    % min/max freqs. for tonotopic map (determines w) [1 5] {[1 4.5]}
P.CFnoise= 1;   % add noise to CF distribution? (determines w) [0] {1}
P.CFnoiseFact= 2*10^(-2);    % factor for noisiness of CF distrib. [10^(-2)] {2*10^(-2)}

P.chunkify= 1; % boolean re making the tonotopic map "chunky" or not 0=no, 1=yes [0]
P.rows= 22; % if P.chunkify=1, approx. how many groups? (since it can be noisy) {30?}
P.rowNoise= 0.2; % adds a bit of variability as to the integer # of oscillators per CF
% ---------------------------------------

% +++
% ensure the number of oscillators is even (for simplicity in further calculations) (why??)
if mod(P.N,2); P.N= P.N+1; end

% ++++ Determine angular frequency (CF) distribution
if (P.CFtype==0),   P.w= linspace(P.frange(1),P.frange(2),P.N)*2*pi;
elseif (P.CFtype==1), P.w= logscale(P.frange(1),P.frange(2),P.N)*2*pi;  end


% ++++ chunkify the tonotopic map
if (P.chunkify==1),
    nn=1; mm= 1;    % indexers
    while nn < P.N
        chunks(mm)= floor(P.N/P.rows+ P.rowNoise*randn(1,1));  % # of oscillators for this CF
        if nn+chunks(mm)>=P.N
            break     % if near the end, break out of while loop
        else
            P.wc(nn:nn+chunks(mm))= mean(P.w(nn:nn+chunks(mm)));  % avg. the CF for that "chunk"
        end
        nn= nn+ chunks(mm); mm= mm+1;
    end
    % take care of any final stragglers at the end
    if (numel(P.wc)<=P.N), P.wc(numel(P.wc)+1:P.N)= mean(P.w(numel(P.wc)+1:P.N));   end
end

% ++++ Add noise to oscillator CFs?
if (P.CFnoise==1),    dw = (8*randn(1,P.N) - 4)*P.CFnoiseFact;     P.w= P.w+dw;  P.wc= P.wc+0.5*dw;  end




% ---------------------------------------
figure(1); clf;
oscN= linspace(1,P.N,P.N);
plot(oscN,P.w/(2*pi),'bx-'); grid on; hold on;
plot(oscN,P.wc/(2*pi),'rs-'); grid on; hold on;
xlabel('Oscillator #'); ylabel('CF [Hz]');
legend('Cont. map (w/ roughness)','Chunky map (w/ rough)','Location','NorthWest');
