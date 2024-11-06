# Audio Processing with MATLAB

This project implements audio signal processing using Discrete Fourier Transform (DFT) in MATLAB. The script processes an audio signal, adds noise, and applies low-pass, high-pass, and band-pass filters to isolate specific frequency ranges.

## Features

- **DFT & IDFT**: Implementation of custom Discrete Fourier Transform (DFT) and Inverse DFT (IDFT) functions.
- **Signal Generation**: Generates a composite signal made up of multiple sine waves of different frequencies (500 Hz, 1000 Hz, and 1500 Hz) with added noise.
- **Filtering**: Applies low-pass, high-pass, and band-pass filters to the noisy signal.
- **Visualization**: Plots the magnitude spectrum of the original and filtered signals.
- **Audio Playback**: Plays the original and filtered audio signals.
- **File Saving**: Saves the original and filtered audio signals to WAV files.

## Files

- `originalaudio.wav`: The original noisy audio signal.
- `audiolow.wav`: Audio after applying a low-pass filter.
- `audiohigh.wav`: Audio after applying a high-pass filter.
- `audioband.wav`: Audio after applying a band-pass filter.

## Filtering Details

- **Low-pass filter**: Allows frequencies below 800 Hz.
- **High-pass filter**: Allows frequencies above 1000 Hz.
- **Band-pass filter**: Allows frequencies between 600 Hz and 1400 Hz.

## How It Works

1. **Signal Generation**:
   - The signal consists of three sine waves (500 Hz, 1000 Hz, 1500 Hz) combined with noise.
   - The noise amplitude can be adjusted via the `noise_amp` variable.
2. **DFT**:
   - The `DFT(x)` function computes the Discrete Fourier Transform of the signal `x`.
   - The frequency spectrum of the signal is analyzed by plotting the magnitude of its DFT.
3. **Filtering**:
   - Three types of filters (low-pass, high-pass, and band-pass) are applied to the DFT-transformed signal.
   - These filters are implemented by zeroing out the undesired frequency components in the Fourier domain.
4. **InverseDFT**:

   - The filtered signals are transformed back to the time domain using the `InverseDFT(X)` function.
   - The filtered signals are normalized and played back.

5. **Visualization**:
   - The script generates time-domain plots for the filtered audio signals as well as magnitude spectra.
