
% test code to visualize dynamically the oscillator dynamics

clear
% ----------------
file= '050815test.mat';


% ----------------

load(file);     % load in data
SR= 1/data.info.params.stepSize;




% ===============================================
% make a "heatmap" for all oscillators across a "steady-state" time
figure(1); clf;
N= 1000;    % length of "time" points from end (thereby assuming "steady-state")
chunk= real(data.Z(end-N:end,:));
chunk= rot90(chunk); % rotate 90 to get proper orientation
chunk= flipud(chunk);
imagesc(chunk); Z= colorbar;
xlabel('Time'); ylabel('Oscillator #'); ylabel(Z,'Inst. Amplitude');

% ===============================================
% make a "heatmap" for all oscillators at an instant, but then animate
% across time
%imagesc(real(data.Z(1,:)));
figure(2); clf; colormap('bone')
for n=1:size(chunk,2)
    imagesc(chunk(:,n));
    pause(0.05);
end