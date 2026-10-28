#ifndef RACE_DSP_ENGINE_H
#define RACE_DSP_ENGINE_H

#include <vector>
#include <mutex>
#include <array>
#include <complex>
#include <cstdint>
#include <fftw3.h>
#include "FirEq.h" // NEU: FIR Equalizer einbinden, phasenlinearer Equalizer! U. Thibaut, 26.10.2026

class FirCrossover {
private:
    int M, L, N;
    int sampleRate;

    float *in_audio, *final_audio_out_low, *final_audio_out_mid, *final_audio_out_high;
    fftwf_complex *out_audio_spec;
    fftwf_complex *multiplied_spec_low, *multiplied_spec_mid, *multiplied_spec_high;

    float *filter_ir_low, *filter_ir_mid, *filter_ir_high;
    fftwf_complex *filter_spec_low, *filter_spec_mid, *filter_spec_high;

    fftwf_plan plan_forward;
    fftwf_plan plan_inverse_low, plan_inverse_mid, plan_inverse_high;
    fftwf_plan plan_filter_low, plan_filter_mid, plan_filter_high;

    std::vector<float> overlap_buffer;

    void generateWindowedSinc(float* ir_buffer, double cutoff1, double cutoff2, bool isBandpass);

public:
    FirCrossover(int blockLength = 32768, int filterLength = 32767, int sr = 192000);
    ~FirCrossover();

    void updateFilters(double freqLow, double freqHigh);
    // NEU: eqCurve Pointer als optionales Argument übergeben
    void processBlock(const float* incoming_audio, float* out_low, float* out_mid, float* out_high, const float* eqCurve = nullptr);
};

// NEU: Hilfsstruktur für die EQ-Parameter
struct EqBand { 
    double f, q, g; 
};

class RACEDspEngine {
private:
    int sampleRate;
    double delayUs, attenuationDb, attenuationLinear, centerP;
    int dn;
    bool freqLimitRACE;
    double currentFreqLow, currentFreqHigh;
    double volume;
    double makeupGainDb;
    double makeupGainLinear;
    bool raceEnabled, filtersEnabled, eqEnabled;

    std::vector<float> spectrumBuffer;
    int spectrumIndex;

    std::vector<float> delayBufferL;
    std::vector<float> delayBufferR;
    int writeIndex;

    FirCrossover firCrossoverL;
    FirCrossover firCrossoverR;

    std::vector<float> inputFifoL;
    std::vector<float> inputFifoR;
    std::vector<float> outputFifoL;
    std::vector<float> outputFifoR;

    std::vector<float> lowL, midL, highL;
    std::vector<float> lowR, midR, highR;

    // NEU: Phasenlinearer EQ anstelle der Biquads, U. Thibaut, 25.10.2026
    FirEq firEqualizer;
    std::array<EqBand, 3> eqBands;

    std::mutex dspMutex;

public:
    RACEDspEngine(double initialDelayUs, double initialAttenuationDb, double initialCenterP, bool initialFreqLimit, int sampleRate = 192000);
    ~RACEDspEngine();

    void setParameters(double newDelayUs, double newAttenuationDb, double newCenterP, bool newFreqLimit);
    void setCrossoverFreqLow(double freq);
    double getCrossoverFreqLow();
    void setCrossoverFreqHigh(double freq);
    double getCrossoverFreqHigh();
    void setVolume(double newVolume);
    void setMakeupGain(double gainDb);
    double getMakeupGain();    
    void setRaceEnabled(bool enabled);
    void setFiltersEnabled(bool enabled);
    void setEqEnabled(bool enabled);
    void setEqBand(int bandIndex, double freq, double q, double gainDb);

    double getDelayUs();
    double getAttenuationDb();
    double getCenterP();
    bool getFreqLimit();
    double getVolume();
    bool getRaceEnabled();
    bool getFiltersEnabled();
    bool getEqEnabled();
    double getEqFreq(int bandIndex);
    double getEqQ(int bandIndex);
    double getEqGain(int bandIndex);

    std::vector<double> getSpectrumBands();
    void processSamples(float* interleavedSamples, size_t frames);
};

#endif // RACE_DSP_ENGINE_H
