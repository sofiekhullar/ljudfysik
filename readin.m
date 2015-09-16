%% Read in sound mp3 file
[y,Fs] = audioread('./sound/father.mp3');
%Play sound
sound(y,Fs);

%% Read in sound wav file
[mantle, fs] = wavread('./sound/Clock_mantle.wav');
sound(mantle, fs);


%% Spelar ljud baklänges
[mantle, fs] = wavread('./sound/Clock_mantle.wav');
d_reverse = flipud(mantle);
sound(d_reverse, fs);

%% Making echo sound

[beautiful, fs] = wavread('./sound/BeautifulLife.wav');
nSec = 90;

b = beautiful(1: fs*nSec);

b_echo = b;
N = fs/2;
for n = N+1 : length(b)
    % adding N off the phase sound to the original input.
    b_echo(n) = b(n) + b(n-N);
end

time = (1/fs)*length(b);
t = linspace(0, time, length(b));

plot(t,b,'k',t,b_echo,'r');  % hard to tell because too many data
xlabel('time(sec)');
ylabel('signal strength');
title('BeautifulLife');

