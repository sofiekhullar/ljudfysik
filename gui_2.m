function varargout = gui_2(varargin)
% GUI_2 MATLAB code for gui_2.fig
%      GUI_2, by itself, creates a new GUI_2 or raises the existing
%      singleton*.
%
%      H = GUI_2 returns the handle to a new GUI_2 or the handle to
%      the existing singleton*.
%
%      GUI_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_2.M with the given input arguments.
%
%      GUI_2('Property','Value',...) creates a new GUI_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_2

% Last Modified by GUIDE v2.5 28-Sep-2015 14:42:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_2_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_2 is made visible.
function gui_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_2 (see VARARGIN)

% Choose default command line output for gui_2
handles.output = hObject;

ah = axes('unit', 'normalized', 'position', [0 0 1 1]); 
bg = imread('bakgrund_1.jpg'); imagesc(bg);
set(ah,'handlevisibility','off','visible','off')
uistack(ah, 'bottom');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[a,map]=imread('robot.jpg');
[r,c,d]=size(a); 
x=ceil(r/260); % Ändrar storleken?
y=ceil(c/360); 
g=a(1:x:end,1:y:end,:);
g(g==255)=5.5*255;
set(handles.robot,'CData',g); 

[a,map]=imread('bebis.jpg');
[r,c,d]=size(a); 
x=ceil(r/260); % Ändrar storleken?
y=ceil(c/360); 
g=a(1:x:end,1:y:end,:);
g(g==255)=5.5*255;
set(handles.pushbutton2,'CData',g); 

[a,map]=imread('monster.jpg');
[r,c,d]=size(a); 
x=ceil(r/260); % Ändrar storleken?
y=ceil(c/360); 
g=a(1:x:end,1:y:end,:);
g(g==255)=5.5*255;
set(handles.monster,'CData',g); 

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Robot 
% --- Executes on button press in bebis.
function robot_Callback(hObject, eventdata, handles)
% hObject    handle to bebis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



%----- record sound -----
fs = 44100;
bits = 16;
recObj = audiorecorder(fs, bits, 1);
%get(recObj)

% Collect a sample of your speech with a microphone, and plot the signal data:
% Record your voice for 5 seconds.
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');

% Store data in double-precision array.
myRecording = getaudiodata(recObj);

% Plot the waveform.
% plot(myRecording);

% wavwrite(myRecording, fs, bits,'rec_robot');
audiowrite('rec_robot.wav', myRecording, fs);

%#wavplay(myRecording,fs);
%----- user data -----
s_win        = 1024;   % analysis window length [samples]
n1           = 441;    % analysis step [samples]
n2           = n1;     % synthesis step [samples]
[DAFx_in,FS] = audioread('rec_robot.wav');

%----- initialize windows, arrays, etc -----
w1       = hanning(s_win, 'periodic'); % analysis window
w2       = w1;    % synthesis window
L        = length(DAFx_in);
DAFx_in  = [zeros(s_win, 1); DAFx_in; ...
  zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_in));
DAFx_out = zeros(length(DAFx_in),1);

tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin  = 0;
pout = 0;
pend = length(DAFx_in)-s_win;
while pin<pend
  grain = DAFx_in(pin+1:pin+s_win).* w1;
%===========================================
  f     = fft(grain);
  r     = abs(f);
  grain = fftshift(real(ifft(r))).*w2;
% ===========================================
  DAFx_out(pout+1:pout+s_win) = ...
    DAFx_out(pout+1:pout+s_win) + grain;
  pin  = pin + n1;
  pout = pout + n2;
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc

%----- listening and saving the output -----
% DAFx_in = DAFx_in(s_win+1:s_win+L);
DAFx_out = DAFx_out(s_win+1:s_win+L) / max(abs(DAFx_out));
soundsc(DAFx_out, FS);
audiowrite('output/robot.wav', DAFx_out, FS);

% --- Bebis
% --- Executes on button press in pushbutton2.
function bebis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----- record sound -----
fs = 44100;
bits = 16;
recObj = audiorecorder(fs, bits, 1);
%get(recObj)

% Collect a sample of your speech with a microphone, and plot the signal data:
% Record your voice for 5 seconds.
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');

% Store data in double-precision array.
myRecording = getaudiodata(recObj);

% Plot the waveform.
% plot(myRecording);

audiowrite('rec_pitch.wav', myRecording, fs);

%#wavplay(myRecording,fs);
%----- user data -----
s_win        = 2048;   % analysis window length [samples]
n2           = 512;    % synthesis step [samples]
pit_ratio    = 1.2     % pitch-shifting ratio
[DAFx_in,FS] = audioread('rec_pitch.wav');

%----- initialize windows, arrays, etc -----
n1       = round(n2 / pit_ratio);      % analysis step [samples]
tstretch_ratio = n2/n1;
w1       = hanning(s_win, 'periodic'); % analysis window
w2       = w1;    % synthesis window
L        = length(DAFx_in);
DAFx_in  = [zeros(s_win, 1); DAFx_in; ...
   zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_in));
DAFx_out = zeros(length(DAFx_in),1);
omega    = 2*pi*n1*[0:s_win-1]'/s_win;
phi0     = zeros(s_win,1);
psi      = zeros(s_win,1);

%----- for linear interpolation of a grain of length s_win -----
lx   = floor(s_win*n1/n2);
x    = 1 + (0:lx-1)'*s_win/lx;
ix   = floor(x);
ix1  = ix + 1;
dx   = x - ix;
dx1  = 1 - dx;


tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin  = 0;
pout = 0;
pend = length(DAFx_in)-s_win;

while pin<pend
  grain = DAFx_in(pin+1:pin+s_win).* w1;
%===========================================
  f     = fft(fftshift(grain));
  r     = abs(f);
  phi   = angle(f);
  %---- computing phase increment ----
  delta_phi = omega + princarg(phi-phi0-omega);
  phi0  = phi;
  psi   = princarg(psi+delta_phi*tstretch_ratio);
  %---- synthesizing time scaled grain ----
  ft    = (r.* exp(i*psi));
  grain = fftshift(real(ifft(ft))).*w2;
  %----- interpolating grain -----
  grain2 = [grain;0];
  grain3 = grain2(ix).*dx1+grain2(ix1).*dx;
  %plot(grain);drawnow;
% ===========================================
  DAFx_out(pout+1:pout+lx) = DAFx_out(pout+1:pout+lx) + grain3;
  pin    = pin + n1;
  pout   = pout + n1;
  end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc

%----- listening and saving the output -----
% DAFx_in  = DAFx_in(s_win+1:s_win+L);
DAFx_out = DAFx_out(s_win+1:s_win+L) / max(abs(DAFx_out));
soundsc(DAFx_out, FS);
audiowrite('output/pitch_pv.wav', DAFx_out, FS);



% --- Bakis johanna
% --- Executes on button press in pushbutton3.
function monster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----- record sound -----
fs = 44100;
bits = 16;
recObj = audiorecorder(fs, bits, 1);
%get(recObj)

% Collect a sample of your speech with a microphone, and plot the signal data:
% Record your voice for 5 seconds.
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');

% Store data in double-precision array.
myRecording = getaudiodata(recObj);

% Plot the waveform.
% plot(myRecording);

audiowrite('rec_whisper.wav', myRecording, fs);

%#wavplay(myRecording,fs);
%----- user data -----
s_win        = 512;     % analysis window length [samples]
n1           = s_win/8; % analysis step [samples]
n2           = n1;      % synthesis step [samples]
[DAFx_in,FS] = audioread('rec_whisper.wav');

%----- initialize windows, arrays, etc -----
w1       = hanning(s_win, 'periodic'); % analysis window
w2       = w1;    % synthesis window
L        = length(DAFx_in);
DAFx_in  = [zeros(s_win, 1); DAFx_in; ...
  zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_in));
DAFx_out = zeros(length(DAFx_in),1);

tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin  = 0;
pout = 0;
pend = length(DAFx_in) - s_win;
while pin<pend
  grain = DAFx_in(pin+1:pin+s_win).* w1;
%===========================================
  f     = fft(fftshift(grain));
  r     = abs(f);
  phi   = 2*pi*rand(s_win,1);
  ft    = (r.* exp(i*phi));
  grain = fftshift(real(ifft(ft))).*w2;
% ===========================================
  DAFx_out(pout+1:pout+s_win) = ...
    DAFx_out(pout+1:pout+s_win) + grain;
  pin   = pin + n1;
  pout  = pout + n2;
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc

%----- listening and saving the output -----
% DAFx_in = DAFx_in(s_win+1:s_win+L);
DAFx_out = DAFx_out(s_win+1:s_win+L) / max(abs(DAFx_out));
soundsc(DAFx_out, FS);
audiowrite('output/whisper2.wav', DAFx_out, FS);

