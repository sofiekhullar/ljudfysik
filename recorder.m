%% sparar ner till wav file

fs = 44100
bits = 16
recObj = audiorecorder(fs, bits, 1);
%# get(recObj)

%# Collect a sample of your speech with a microphone, and plot the signal data:
%# Record your voice for 5 seconds.
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');

%# Play back the recording.
play(recObj);

%# Store data in double-precision array.
myRecording = getaudiodata(recObj);
%disp(size(myRecording));

%# Plot the waveform.
plot(myRecording);

wavwrite(myRecording, fs, bits,'sample01_6k');
%#wavplay(myRecording,fs);