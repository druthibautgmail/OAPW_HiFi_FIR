#ifndef FIREQ_H
#define FIREQ_H

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>

class FirEq {
public:
    // fftSize ist die Anzahl der Real-Samples (z.B. 32768)
    FirEq(int sampleRate, int fftSize);

    // Berechnet die Magnitude-Kurve für ein Peaking-Band (analog zu Biquads, aber phasenlinear)
    void addPeakingBand(float centerFreq, float gainDb, float q);
    
    // Setzt die Kurve wieder auf komplett flat (0 dB)
    void reset();

    // Wird in RACEDspEngine direkt nach der Forward-FFT aufgerufen
    void applyToFFT(std::vector<std::complex<float>>& complexData);

    // Wendet FFT statt Biquad Filtern an! Phasenkorrekt, U. Thibaut, 25.10.2026
    void applyToFFTW(fftwf_complex* fftwData);

    // Getter Funktion für den Equalizer auf FFT-Basis, U. Thibaut, 25.10.2026
    const float* getCurveData() const { return eqCurve.data(); }
private:
    int sr;
    int nFft;
    int numBins;
    std::vector<float> eqCurve; // Speichert die linearen Gain-Faktoren für jeden FFT-Bin
};

#endif
