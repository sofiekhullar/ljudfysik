     % LOAD DATA
	  load mtlb;
	  x = mtlb;

      [DAFx_in,FS] = audioread('rec_robot.wav');
      
	  figure(1), clf
	  plot(0:220499,DAFx_in)
	  xlabel('n')
	  ylabel('x(n)')

	  % SET PARAMETERS
	  R = 256;               % R: block length
	  window = hamming(R);   % window function of length R
	  N = 512;               % N: frequency discretization
	  L = 35;                % L: time lapse between blocks
	  fs = 7418;             % fs: sampling frequency
	  overlap = R - L;

	  % COMPUTE SPECTROGRAM
	  [B,f,t] = specgram(DAFx_in,N,fs,window,overlap);

	  % MAKE PLOT
	  figure(2), clf
	  imagesc(t,f,log10(abs(B)));
	  colormap('jet')
	  axis xy 
	  xlabel('time')
	  ylabel('frequency')
	  title('SPECTROGRAM, R = 256')