function X = DFT(x)
    N = length(x);
    X = zeros(N, 1);
    for k = 0:N-1
        for n = 0:N-1
            X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
        end
    end
end
function x = InverseDFT(X)
    N = length(X);
    x = zeros(N, 1);
    for n = 0:N-1
        for k = 0:N-1
            x(n+1) = x(n+1) + X(k+1) * exp(1j * 2 * pi * k * n / N);
        end
        x(n+1) = x(n+1) / N;
    end
end
Fs = 8000;
t = 0:1/Fs:1-1/Fs;
f1 = 500;
f2 = 1000;
f3 = 1500;
noise_amp = 0.2;
signal = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t) + 0.3*sin(2*pi*f3*t);
noise = noise_amp * randn(size(t));
audio = signal + noise;
plot(t, audio);
audiowrite('originalaudio.wav', audio, Fs);
sound(audio, Fs);
X = DFT(audio);
f = (0:length(audio)-1) * (Fs / length(audio));
figure;
plot(f, abs(X));
title('Magnitude Spectrum of the Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
cutoff_low = 800;
cutoff_high = 1000;
cutoff_band_low = 600;
cutoff_band_high = 1400;
H_low = double(f <= cutoff_low)';
H_high = double(f >= cutoff_high)';
H_band = double(f >= cutoff_band_low & f <= cutoff_band_high)';
X_low = X .* H_low;
X_high = X .* H_high;
X_band = X .* H_band;
audio_low = InverseDFT(X_low);
audio_high = InverseDFT(X_high);
audio_band = InverseDFT(X_band);
audio_low = real(audio_low) / max(abs(audio_low));
audio_high = real(audio_high) / max(abs(audio_high));
audio_band = real(audio_band) / max(abs(audio_band));
audiowrite('audiolow.wav', audio_low, Fs);
audiowrite('audiohigh.wav', audio_high, Fs);
audiowrite('audioband.wav', audio_band, Fs);
sound(audio_low, Fs);
pause(length(audio)/Fs + 2);
sound(audio_high, Fs);
pause(length(audio)/Fs + 2);
sound(audio_band, Fs);
figure;
subplot(3,1,1);
plot(f, abs(DFT(audio_low)));
title('Low-pass Filtered Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
subplot(3,1,2);
plot(f, abs(DFT(audio_high)));
title('High-pass Filtered Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
subplot(3,1,3);
plot(f, abs(DFT(audio_band)));
title('Band-pass Filtered Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
figure;
subplot(3,1,1);
plot(t, audio_low);
title('Low-pass Filtered Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(3,1,2);
plot(t, audio_high);
title('High-pass Filtered Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(3,1,3);
plot(t, audio_band);
title('Band-pass Filtered Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');