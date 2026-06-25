# Linear Predictive Coding (LPC) Speech Coding and Synthesis

**Course:** Digital Audio Analysis and Processing

**Programme:** MSc in Music and Acoustic Engineering, Politecnico di Milano

**Work Team**: Di Lorenzo Giuliano, Longhi Filippo

## Overview

This project implements a speech coding and synthesis system based on Linear Predictive Coding (LPC), inspired by the speech synthesis approach used in the Texas Instruments *Speak & Spell* educational toy.

The goal of the project is to analyze speech signals, extract a compact LPC representation, and reconstruct the original audio using the estimated vocal tract model and excitation signal. LPC allows efficient speech representation by modeling the speech production mechanism through a prediction filter and an excitation source.

The implementation was developed in MATLAB2024b and includes both encoding and decoding stages.

## Methodology

The LPC system is divided into two main stages:

1. **Encoding**
2. **Decoding**

The speech signal is processed frame-by-frame, assuming local stationarity of speech. Each frame is analyzed to estimate the LPC parameters that describe the spectral envelope of the signal.

## LPC Analysis

The speech frame is modeled as:

$$ s(n) \approx \sum_{k=1}^{p} a_k s(n-k) $$

where $s(n)$ is the speech signal, $a_k$ are the LPC coefficients, and $p$ is the prediction order.

The LPC coefficients are obtained from the autocorrelation of each frame using the Levinson-Durbin recursion.

Different prediction orders are used depending on the speech characteristics: $p = 10$ for voiced frames, $p = 4$ for unvoiced frames.

Voiced and unvoiced classification allows the system to select an appropriate excitation model during synthesis.

## Encoding Process

The encoding stage performs the following operations:

1. **Frame extraction**: the signal is divided into short frames using a Hamming window, which allows independent analysis of each speech segment.

2. **LPC coefficient computation**: for every frame, the autocorrelation function is computed, LPC coefficients are estimated using the Levinson-Durbin algorithm, and the prediction error is obtained through the whitening filter.
The prediction error represents the excitation information required for reconstruction.

3. **Pitch estimation**: for voiced frames, the pitch period is estimated from the prediction error using autocorrelation analysis.
The pitch information is then used to generate the excitation signal during decoding.

## Decoding Process

The reconstruction stage recreates the speech signal from the stored LPC parameters.

Two excitation models are used:

1. **Voiced frames**: the excitation signal is generated as an impulse train

    - Impulse spacing corresponds to the estimated pitch period.
    - The amplitude is controlled by the estimated gain.

2. **Unvoiced frames**: the excitation signal is modeled as white noise

    - Zero-mean Gaussian noise
    - Variance controlled by the LPC gain estimation

## Speech Reconstruction

The excitation signal is passed through the LPC synthesis filter:

$$ H(z)=\frac{1}{1-\sum_{k=1}^{p}a_kz^{-k}} $$

The reconstructed frames are then combined to obtain the final synthesized speech signal.

## Analysis Results

The project includes analysis plots comparing:

- Speech spectrum and LPC shaping filter response
- Prediction error in the time domain
- Prediction error spectrum

These plots allow evaluation of how accurately the LPC model represents voiced and unvoiced speech components.

## Synthesis Results

The LPC model successfully reproduces the general characteristics of the original speech signal.

Observations:

- Voiced sounds maintain their harmonic structure thanks to the pitch-based excitation model.
- Unvoiced sounds are reconstructed using noise excitation and maintain their spectral behavior.
- The reconstructed signal preserves the temporal evolution of the original speech.

Some differences remain, particularly in spectral detail and consonant clarity, due to the simplified excitation model and LPC approximation.
