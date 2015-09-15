%Read in sound
[y,Fs] = audioread('father.mp3');

whos y

%Play sound
sound(y,Fs);

