%for i = 1:17
  %  figure
 %   plot(S(:,i))
%end
[S HDR] = sload('c30o30_s4_t1.bdf');
Fs = HDR.SampleRate;
% Implements supplementary_code/part_2/extract_epochs.m
% Find length of signal
L = length(S);
% Differentiate trigger channel
delta_trig = S(2:L,17) - S(1:L-1,17);
% Look for rising edges
trig = find(delta_trig>0);
% Debounce if required
% I.e. remove spurious short triggers, if any
% these can occur when using a manual switch-based trigger
% e.g. look at the trigger channel in 'data_files/c30o30_s4_t1.bdf'
T_debounce = 0.5; % (in seconds) - set as desired
dt_trig = trig(2:end) - trig(1:end-1);
delete_short_trigs = find(dt_trig < Fs*T_debounce) + 1;
trig(delete_short_trigs) = [];
% How many triggers were extracted in total, after debouncing?
N_trigs = length(trig)-1;
% Ensure we have an odd number of triggers (ignore final incomplete
% epoch if not).
% Trigger 1 = start of experiment
% Triggers 2, 4, 6, 8... (2k) = eyes closed
% Triggers 3, 5, 7, 9... (2k+1) = eyes open
if mod(N_trigs, 2) == 0
 N_trigs = N_trigs - 1;
end
disp(cell2mat(strcat({'Found '}, int2str(N_trigs), {' triggers.'})));

N_epochs = N_trigs;
% Discard 2 seconds of data before & after each cue
discard = 2*Fs;
% Extract each epoch of data and store in the relevant structure for
% each class
% - eyes_closed (even triggers)
% - eyes_open (odd triggers)
for i=1:N_epochs/2
 eyes_closed{i} = S(trig(2*i)+discard:trig(2*i+1)-discard,1:16);
 eyes_open{i} = S(trig(2*i+1)+discard:trig(2*i+2)-discard,1:16);
end

%11.
% Specify the electrode montage on scalp (the electrode numbers laid out 
% in a matrix, according to their positions on the scalp and as we want 
% them to be displayed
montage = [-1 -1 1 -1 -1; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];

% The 10/20 labels for each electrode channel (again, corresponding to the 
% montage above. In thsi case electrode 1 = Fz, electrode 2 = FC3, 
% electrode 3 = FC1 etc.
electrode_labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4','C3', 'C1',... 
                       'Cz', 'C2', 'C4','CP3', 'CP1', 'CPz', 'CP2', 'CP4'};

% Create a new figure (or select an exisiting one if you prefer)

%plotFFT2(eyes_closed_FFT(:,electrode), Fs)

figure

% Iterate through each row (j) and column(i) of the montage matrix that
% will form your subplot figure
for j=1:4
    for i=1:5
        % Only create plots for electrodes that exist (ignore -1 values)
        if montage(j, i) > 0
            % which electrode is at this location in the montage matrix?
            electrode = montage(j, i);
            % select the correct (row, colum) --> what does the number 5 do
            % here?
            subplot(4, 5, i+(5*(j-1)))
            
            % plot something for the current electrode. You can plot
            % whatever data you want here, but for example:
            
            %Fourier Transfer:
            eyes_closed_FFT(:,electrode) = ssfft(eyes_closed{1}(:,electrode));
            plotFFT2(eyes_closed_FFT(:,electrode), Fs);
            
            % Set axes limits if necessary, e.g.
            ylim([0 2.5])
            xlim([0 50])
            
            % Label each subplot with the corresponding electrode poistion
            title(electrode_labels{electrode})
            
            % Label the axes for the bottom left plot only (to prevent the
            % figure getting too crowded and illegible)
            if i == 1 && j == 4
                xlabel('Frequency (Hz)');
                ylabel('FFT (V)')
            end
        end
    end
end

% Add a title to the entire plot - play around, as required
annotation('textbox', [0.1,0.75, 0.2, 0.2], 'String', 'Eyes Closed',...
                                'FontWeight', 'bold', 'LineStyle', 'none')



%%
%for i = 1:17
  %  figure
 %   plot(S(:,i))
%end
[S HDR] = sload('c30o30_s4_t1.bdf');
Fs = HDR.SampleRate;
% Implements supplementary_code/part_2/extract_epochs.m
% Find length of signal
L = length(S);
% Differentiate trigger channel
delta_trig = S(2:L,17) - S(1:L-1,17);
% Look for rising edges
trig = find(delta_trig>0);
% Debounce if required
% I.e. remove spurious short triggers, if any
% these can occur when using a manual switch-based trigger
% e.g. look at the trigger channel in 'data_files/c30o30_s4_t1.bdf'
T_debounce = 0.5; % (in seconds) - set as desired
dt_trig = trig(2:end) - trig(1:end-1);
delete_short_trigs = find(dt_trig < Fs*T_debounce) + 1;
trig(delete_short_trigs) = [];
% How many triggers were extracted in total, after debouncing?
N_trigs = length(trig)-1;
% Ensure we have an odd number of triggers (ignore final incomplete
% epoch if not).
% Trigger 1 = start of experiment
% Triggers 2, 4, 6, 8... (2k) = eyes closed
% Triggers 3, 5, 7, 9... (2k+1) = eyes open
if mod(N_trigs, 2) == 0
 N_trigs = N_trigs - 1;
end
disp(cell2mat(strcat({'Found '}, int2str(N_trigs), {' triggers.'})));

N_epochs = N_trigs;
% Discard 2 seconds of data before & after each cue
discard = 2*Fs;
% Extract each epoch of data and store in the relevant structure for
% each class
% - eyes_closed (even triggers)
% - eyes_open (odd triggers)
for i=1:N_epochs/2
 eyes_closed{i} = S(trig(2*i)+discard:trig(2*i+1)-discard,1:16);
 eyes_open{i} = S(trig(2*i+1)+discard:trig(2*i+2)-discard,1:16);
end

%11.
% Specify the electrode montage on scalp (the electrode numbers laid out 
% in a matrix, according to their positions on the scalp and as we want 
% them to be displayed
montage = [-1 -1 1 -1 -1; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];

% The 10/20 labels for each electrode channel (again, corresponding to the 
% montage above. In thsi case electrode 1 = Fz, electrode 2 = FC3, 
% electrode 3 = FC1 etc.
electrode_labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4','C3', 'C1',... 
                       'Cz', 'C2', 'C4','CP3', 'CP1', 'CPz', 'CP2', 'CP4'};

% Create a new figure (or select an exisiting one if you prefer)

%plotFFT2(eyes_closed_FFT(:,electrode), Fs)

figure

% Iterate through each row (j) and column(i) of the montage matrix that
% will form your subplot figure
for j=1:4
    for i=1:5
        % Only create plots for electrodes that exist (ignore -1 values)
        if montage(j, i) > 0
            % which electrode is at this location in the montage matrix?
            electrode = montage(j, i);
            % select the correct (row, colum) --> what does the number 5 do
            % here?
            subplot(4, 5, i+(5*(j-1)))
            
            % plot something for the current electrode. You can plot
            % whatever data you want here, but for example:
            
            %Fourier Transfer:
            eyes_open_FFT(:,electrode) = ssfft(eyes_open{1}(:,electrode));
            plotFFT2(eyes_open_FFT(:,electrode), Fs);
            
            % Set axes limits if necessary, e.g.
            ylim([0 2.5])
            xlim([0 50])
            
            % Label each subplot with the corresponding electrode poistion
            title(electrode_labels{electrode})
            
            % Label the axes for the bottom left plot only (to prevent the
            % figure getting too crowded and illegible)
            if i == 1 && j == 4
                xlabel('Frequency (Hz)');
                ylabel('FFT (V)')
            end
        end
    end
end

% Add a title to the entire plot - play around, as required
annotation('textbox', [0.1,0.75, 0.2, 0.2], 'String', 'Eyes Open',...
                                'FontWeight', 'bold', 'LineStyle', 'none')
                            
%%
%13.
figure
PSD_eyes_closed = (eyes_open_FFT(:,electrode)).*conj(eyes_closed_FFT(:, electrode));
PSD_eyes_closed = (abs(eyes_open_FFT(:,electrode)).^2);
plot(Fs,PSD_eyes_closed,'LineWidth',1)
ylabel('PSD(Eyes Closed)')
xlabel('Frequency (Hz)')
title('Spectrum of Eyes Closed')

%%
%13.
figure
PSD_eyes_open = (eyes_open_FFT(:,electrode)).*conj(eyes_open_FFT(:, electrode));
plot(Fs,PSD_eyes_open)
ylabel('PSD(Eyes Open)')
xlabel('Frequency (Hz)')
title('Spectrum of Eyes Open')
%%

%%
%for i = 1:17
  %  figure
 %   plot(S(:,i))
%end
[S HDR] = sload('c30o30_s4_t1.bdf');
Fs = HDR.SampleRate;
% Implements supplementary_code/part_2/extract_epochs.m
% Find length of signal
L = length(S);
% Differentiate trigger channel
delta_trig = S(2:L,17) - S(1:L-1,17);
% Look for rising edges
trig = find(delta_trig>0);
% Debounce if required
% I.e. remove spurious short triggers, if any
% these can occur when using a manual switch-based trigger
% e.g. look at the trigger channel in 'data_files/c30o30_s4_t1.bdf'
T_debounce = 0.5; % (in seconds) - set as desired
dt_trig = trig(2:end) - trig(1:end-1);
delete_short_trigs = find(dt_trig < Fs*T_debounce) + 1;
trig(delete_short_trigs) = [];
% How many triggers were extracted in total, after debouncing?
N_trigs = length(trig)-1;
% Ensure we have an odd number of triggers (ignore final incomplete
% epoch if not).
% Trigger 1 = start of experiment
% Triggers 2, 4, 6, 8... (2k) = eyes closed
% Triggers 3, 5, 7, 9... (2k+1) = eyes open
if mod(N_trigs, 2) == 0
 N_trigs = N_trigs - 1;
end
disp(cell2mat(strcat({'Found '}, int2str(N_trigs), {' triggers.'})));

N_epochs = N_trigs;
% Discard 2 seconds of data before & after each cue
discard = 2*Fs;
% Extract each epoch of data and store in the relevant structure for
% each class
% - eyes_closed (even triggers)
% - eyes_open (odd triggers)
for i=1:N_epochs/2
 eyes_closed{i} = S(trig(2*i)+discard:trig(2*i+1)-discard,1:16);
 eyes_open{i} = S(trig(2*i+1)+discard:trig(2*i+2)-discard,1:16);
end

%11.
% Specify the electrode montage on scalp (the electrode numbers laid out 
% in a matrix, according to their positions on the scalp and as we want 
% them to be displayed
montage = [-1 -1 1 -1 -1; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];

% The 10/20 labels for each electrode channel (again, corresponding to the 
% montage above. In thsi case electrode 1 = Fz, electrode 2 = FC3, 
% electrode 3 = FC1 etc.
electrode_labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4','C3', 'C1',... 
                       'Cz', 'C2', 'C4','CP3', 'CP1', 'CPz', 'CP2', 'CP4'};

% Create a new figure (or select an exisiting one if you prefer)

%plotFFT2(eyes_closed_FFT(:,electrode), Fs)

figure

% Iterate through each row (j) and column(i) of the montage matrix that
% will form your subplot figure
for x=1:5
for j=1:4
    figure(x)
    for i=1:5
        % Only create plots for electrodes that exist (ignore -1 values)
        if montage(j, i) > 0
            % which electrode is at this location in the montage matrix?
            electrode = montage(j, i);
            % select the correct (row, colum) --> what does the number 5 do
            % here?
            subplot(4, 5, i+(5*(j-1)))
            
            % plot something for the current electrode. You can plot
            % whatever data you want here, but for example:
            
            %Fourier Transfer:
            eyes_open_FFT(:,electrode) = ssfft(eyes_open{1}(:,electrode));
            PSD_eyes_open(:,electrode) = eyes_open_FFT(:,electrode).*conj(eyes_open_FFT(:, electrode));
            plotFFT2(PSD_eyes_open, Fs);
            
            % Set axes limits if necessary, e.g.
            ylim([0 6])
            xlim([0 50])
            
            % Label each subplot with the corresponding electrode poistion
            title(electrode_labels{electrode})
            
            % Label the axes for the bottom left plot only (to prevent the
            % figure getting too crowded and illegible)
            if i == 1 && j == 4
                xlabel('Frequency (Hz)');
                ylabel('||FFT(Eyes Open)||')
            end
        end
    end
end

% Add a title to the entire plot - play around, as required
annotation('textbox', [0.1,0.75, 0.2, 0.2], 'String', 'PSD Eyes Open',...
                                'FontWeight', 'bold', 'LineStyle', 'none')
                            end
%%
%%
%for i = 1:17
  %  figure
 %   plot(S(:,i))
%end
[S HDR] = sload('c30o30_s4_t1.bdf');
Fs = HDR.SampleRate;
% Implements supplementary_code/part_2/extract_epochs.m
% Find length of signal
L = length(S);
% Differentiate trigger channel
delta_trig = S(2:L,17) - S(1:L-1,17);
% Look for rising edges
trig = find(delta_trig>0);
% Debounce if required
% I.e. remove spurious short triggers, if any
% these can occur when using a manual switch-based trigger
% e.g. look at the trigger channel in 'data_files/c30o30_s4_t1.bdf'
T_debounce = 0.5; % (in seconds) - set as desired
dt_trig = trig(2:end) - trig(1:end-1);
delete_short_trigs = find(dt_trig < Fs*T_debounce) + 1;
trig(delete_short_trigs) = [];
% How many triggers were extracted in total, after debouncing?
N_trigs = length(trig)-1;
% Ensure we have an odd number of triggers (ignore final incomplete
% epoch if not).
% Trigger 1 = start of experiment
% Triggers 2, 4, 6, 8... (2k) = eyes closed
% Triggers 3, 5, 7, 9... (2k+1) = eyes open
if mod(N_trigs, 2) == 0
 N_trigs = N_trigs - 1;
end
disp(cell2mat(strcat({'Found '}, int2str(N_trigs), {' triggers.'})));

N_epochs = N_trigs;
% Discard 2 seconds of data before & after each cue
discard = 2*Fs;
% Extract each epoch of data and store in the relevant structure for
% each class
% - eyes_closed (even triggers)
% - eyes_open (odd triggers)
for i=1:N_epochs/2
 eyes_closed{i} = S(trig(2*i)+discard:trig(2*i+1)-discard,1:16);
 eyes_open{i} = S(trig(2*i+1)+discard:trig(2*i+2)-discard,1:16);
end

%11.
% Specify the electrode montage on scalp (the electrode numbers laid out 
% in a matrix, according to their positions on the scalp and as we want 
% them to be displayed
montage = [-1 -1 1 -1 -1; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];

% The 10/20 labels for each electrode channel (again, corresponding to the 
% montage above. In thsi case electrode 1 = Fz, electrode 2 = FC3, 
% electrode 3 = FC1 etc.
electrode_labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4','C3', 'C1',... 
                       'Cz', 'C2', 'C4','CP3', 'CP1', 'CPz', 'CP2', 'CP4'};

% Create a new figure (or select an exisiting one if you prefer)

%plotFFT2(eyes_closed_FFT(:,electrode), Fs)



% Iterate through each row (j) and column(i) of the montage matrix that
% will form your subplot figure
for x=1:5%this loops through the epochs - time variations
    figure(x)
for j=1:4 
    for i=1:5
        % Only create plots for electrodes that exist (ignore -1 values)
        if montage(j, i) > 0 %this loops through electrodes
            % which electrode is at this location in the montage matrix?
            electrode = montage(j, i);
            % select the correct (row, colum) --> what does the number 5 do
            % here?
            subplot(4, 5, i+(5*(j-1)))
            
            % plot something for the current electrode. You can plot
            % whatever data you want here, but for example:
            
            %Fourier Transfer:
            eyes_closed_FFT(:,electrode) = ssfft(eyes_closed{1}(:,electrode));
            PSD_eyes_closed (:,electrode)= eyes_closed_FFT(:,electrode).*conj(eyes_closed_FFT(:, electrode));
            plotFFT2(PSD_eyes_closed, Fs);
            
            % Set axes limits if necessary, e.g.
            ylim([0 6])
            xlim([0 50])
            
            % Label each subplot with the corresponding electrode poistion
            title(electrode_labels{electrode})
            
            % Label the axes for the bottom left plot only (to prevent the
            % figure getting too crowded and illegible)
            if i == 1 && j == 4
                xlabel('Frequency (Hz)');
                ylabel('||FFT(Eyes Closed)||')
            end
        end
    end
end

% Add a title to the entire plot - play around, as required
annotation('textbox', [0.1,0.75, 0.2, 0.2], 'String', 'PSD Eyes Closed',...
                                'FontWeight', 'bold', 'LineStyle', 'none')
end

%%
%%
%for i = 1:17
  %  figure
 %   plot(S(:,i))
%end
[S HDR] = sload('c30o30_s4_t1.bdf');
Fs = HDR.SampleRate;
% Implements supplementary_code/part_2/extract_epochs.m
% Find length of signal
L = length(S);
% Differentiate trigger channel
delta_trig = S(2:L,17) - S(1:L-1,17);
% Look for rising edges
trig = find(delta_trig>0);
% Debounce if required
% I.e. remove spurious short triggers, if any
% these can occur when using a manual switch-based trigger
% e.g. look at the trigger channel in 'data_files/c30o30_s4_t1.bdf'
T_debounce = 0.5; % (in seconds) - set as desired
dt_trig = trig(2:end) - trig(1:end-1);
delete_short_trigs = find(dt_trig < Fs*T_debounce) + 1;
trig(delete_short_trigs) = [];
% How many triggers were extracted in total, after debouncing?
N_trigs = length(trig)-1;
% Ensure we have an odd number of triggers (ignore final incomplete
% epoch if not).
% Trigger 1 = start of experiment
% Triggers 2, 4, 6, 8... (2k) = eyes closed
% Triggers 3, 5, 7, 9... (2k+1) = eyes open
if mod(N_trigs, 2) == 0
 N_trigs = N_trigs - 1;
end
disp(cell2mat(strcat({'Found '}, int2str(N_trigs), {' triggers.'})));

N_epochs = N_trigs;
% Discard 2 seconds of data before & after each cue
discard = 2*Fs;
% Extract each epoch of data and store in the relevant structure for
% each class
% - eyes_closed (even triggers)
% - eyes_open (odd triggers)
for i=1:N_epochs/2
 eyes_closed{i} = S(trig(2*i)+discard:trig(2*i+1)-discard,1:16);
 eyes_open{i} = S(trig(2*i+1)+discard:trig(2*i+2)-discard,1:16);
end

%11.
% Specify the electrode montage on scalp (the electrode numbers laid out 
% in a matrix, according to their positions on the scalp and as we want 
% them to be displayed
montage = [-1 -1 1 -1 -1; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];

% The 10/20 labels for each electrode channel (again, corresponding to the 
% montage above. In thsi case electrode 1 = Fz, electrode 2 = FC3, 
% electrode 3 = FC1 etc.
electrode_labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4','C3', 'C1',... 
                       'Cz', 'C2', 'C4','CP3', 'CP1', 'CPz', 'CP2', 'CP4'};

% Create a new figure (or select an exisiting one if you prefer)

%plotFFT2(eyes_closed_FFT(:,electrode), Fs)



% Iterate through each row (j) and column(i) of the montage matrix that
% will form your subplot figure

for x = 1:5
    figure(x)
for j=1:4
    for i=1:5
        % Only create plots for electrodes that exist (ignore -1 values)
        if montage(j, i) > 0
            % which electrode is at this location in the montage matrix?
            electrode = montage(j, i);
            % select the correct (row, colum) --> what does the number 5 do
            % here?
            subplot(4, 5, i+(5*(j-1)))
            
            % plot something for the current electrode. You can plot
            % whatever data you want here, but for example:
            
            %Fourier Transfer:
            eyes_closed_FFT(:,electrode) = ssfft(eyes_closed{x}(:,electrode));
            %PSD_eyes_closed = eyes_closed_FFT(:,electrode).*conj(eyes_closed_FFT(:, electrode));
            PSD_eyes_closed(:,electrode)= (abs(eyes_open_FFT(:,electrode)).^2);
            avg_PSD_eyes_closed(:,electrode) = mean(PSD_eyes_closed(:, electrode), 2);
            plotFFT2(PSD_eyes_closed, Fs);
            
            % Set axes limits if necessary, e.g.
            ylim([0 6])
            xlim([0 50])
            
            % Label each subplot with the corresponding electrode poistion
            title(electrode_labels{electrode})
            
            % Label the axes for the bottom left plot only (to prevent the
            % figure getting too crowded and illegible)
            if i == 1 && j == 4
                xlabel('Frequency (Hz)');
                ylabel('||FFT(Eyes Closed)||')
            end
        end
    end
end

% Add a title to the entire plot - play around, as required
annotation('textbox', [0.1,0.75, 0.2, 0.2], 'String', 'PSD Eyes Closed',...
                                'FontWeight', 'bold', 'LineStyle', 'none')
end
%M = mean(PSD_eyes_closed,vecdim)
%mean(PSD_eyes_closed, [0 2])

