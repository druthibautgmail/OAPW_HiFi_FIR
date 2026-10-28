#include "RACEDspEngine.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------
// FFT Helper für den Web-Analyzer
// ---------------------------------------------------------
void computeFFT(std::vector<std::complex<float>>& data) {
    size_t n = data.size();
    if (n <= 1) return;
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        float angle = -2.0f * static_cast<float>(M_PI) / static_cast<float>(len);
        std::complex<float> wlen(std::cos(angle), std::sin(angle));
        for (size_t i = 0; i < n; i += len) {
            std::complex<float> w(1, 0);
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<float> u = data[i + j];
                std::complex<float> v = data[i + j + len / 2] * w;
                data[i + j] = u + v;
                data[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

//
// Hier war bisher der Code für die Biquad-Equalizer Funktionen eingebaut, gelöscht am 26.10.2026, U. Thibaut
//

// ---------------------------------------------------------
// FirCrossover Implementierung
// ---------------------------------------------------------
double blackmanHarris(int n, int filterLength) {
    double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;
    double arg = (2.0 * M_PI * n) / (filterLength - 1);
    return a0 - a1 * std::cos(arg) + a2 * std::cos(2.0 * arg) - a3 * std::cos(3.0 * arg);
}

FirCrossover::FirCrossover(int blockLength, int filterLength, int sr) 
    : M(filterLength), L(blockLength), N(65536), sampleRate(sr) 
{
    in_audio = fftwf_alloc_real(N);
    out_audio_spec = fftwf_alloc_complex(N / 2 + 1);
    
    final_audio_out_low = fftwf_alloc_real(N);
    final_audio_out_mid = fftwf_alloc_real(N);
    final_audio_out_high = fftwf_alloc_real(N);

    filter_ir_low = fftwf_alloc_real(N);
    filter_ir_mid = fftwf_alloc_real(N);
    filter_ir_high = fftwf_alloc_real(N);

    filter_spec_low = fftwf_alloc_complex(N / 2 + 1);
    filter_spec_mid = fftwf_alloc_complex(N / 2 + 1);
    filter_spec_high = fftwf_alloc_complex(N / 2 + 1);

    multiplied_spec_low = fftwf_alloc_complex(N / 2 + 1);
    multiplied_spec_mid = fftwf_alloc_complex(N / 2 + 1);
    multiplied_spec_high = fftwf_alloc_complex(N / 2 + 1);

    overlap_buffer.resize(M - 1, 0.0f);

    plan_forward = fftwf_plan_dft_r2c_1d(N, in_audio, out_audio_spec, FFTW_ESTIMATE);
    plan_inverse_low = fftwf_plan_dft_c2r_1d(N, multiplied_spec_low, final_audio_out_low, FFTW_ESTIMATE);
    plan_inverse_mid = fftwf_plan_dft_c2r_1d(N, multiplied_spec_mid, final_audio_out_mid, FFTW_ESTIMATE);
    plan_inverse_high = fftwf_plan_dft_c2r_1d(N, multiplied_spec_high, final_audio_out_high, FFTW_ESTIMATE);
    plan_filter_low = fftwf_plan_dft_r2c_1d(N, filter_ir_low, filter_spec_low, FFTW_ESTIMATE);
    plan_filter_mid = fftwf_plan_dft_r2c_1d(N, filter_ir_mid, filter_spec_mid, FFTW_ESTIMATE);
    plan_filter_high = fftwf_plan_dft_r2c_1d(N, filter_ir_high, filter_spec_high, FFTW_ESTIMATE);
}

FirCrossover::~FirCrossover() {
    fftwf_destroy_plan(plan_forward);
    fftwf_destroy_plan(plan_inverse_low);
    fftwf_destroy_plan(plan_inverse_mid);
    fftwf_destroy_plan(plan_inverse_high);
    fftwf_destroy_plan(plan_filter_low);
    fftwf_destroy_plan(plan_filter_mid);
    fftwf_destroy_plan(plan_filter_high);
    
    fftwf_free(in_audio); fftwf_free(out_audio_spec);
    fftwf_free(final_audio_out_low); fftwf_free(final_audio_out_mid); fftwf_free(final_audio_out_high);
    fftwf_free(filter_ir_low); fftwf_free(filter_ir_mid); fftwf_free(filter_ir_high);
    fftwf_free(filter_spec_low); fftwf_free(filter_spec_mid); fftwf_free(filter_spec_high);
    fftwf_free(multiplied_spec_low); fftwf_free(multiplied_spec_mid); fftwf_free(multiplied_spec_high);
}

void FirCrossover::generateWindowedSinc(float* ir_buffer, double cutoff1, double cutoff2, bool isBandpass) {
    double fc = cutoff1 / sampleRate;
    double sum = 0.0;
    int center = (M - 1) / 2;

    for (int i = 0; i < M; i++) {
        if (i == center) {
            ir_buffer[i] = 2.0 * M_PI * fc;
        } else {
            double x = 2.0 * M_PI * fc * (i - center);
            ir_buffer[i] = std::sin(x) / (i - center);
        }
        ir_buffer[i] *= blackmanHarris(i, M);
        sum += ir_buffer[i];
    }

    for (int i = 0; i < M; i++) ir_buffer[i] /= sum;

    if (cutoff2 > 0.0 && !isBandpass) { 
        for (int i = 0; i < M; i++) ir_buffer[i] = -ir_buffer[i];
        ir_buffer[center] += 1.0;
    }
    
    for (int i = M; i < N; i++) ir_buffer[i] = 0.0f;
}

void FirCrossover::updateFilters(double freqLow, double freqHigh) {
    generateWindowedSinc(filter_ir_low, freqLow, 0.0, false);
    generateWindowedSinc(filter_ir_high, freqHigh, 1.0, false); 

    int center = (M - 1) / 2;
    for (int i = 0; i < N; i++) {
        float dirac = (i == center) ? 1.0f : 0.0f;
        filter_ir_mid[i] = dirac - filter_ir_low[i] - filter_ir_high[i];
    }

    fftwf_execute(plan_filter_low);
    fftwf_execute(plan_filter_mid);
    fftwf_execute(plan_filter_high);
}

void FirCrossover::processBlock(const float* incoming_audio, float* out_low, float* out_mid, float* out_high, const float* eqCurve) {
    std::memcpy(in_audio, overlap_buffer.data(), (M - 1) * sizeof(float));
    std::memcpy(in_audio + (M - 1), incoming_audio, L * sizeof(float));
    for (int i = (M - 1 + L); i < N; i++) in_audio[i] = 0.0f; 
    
    std::memcpy(overlap_buffer.data(), incoming_audio + L - (M - 1), (M - 1) * sizeof(float));

    fftwf_execute(plan_forward);

    int numBins = N / 2 + 1;
    for (int i = 0; i < numBins; i++) {
        float a_re = out_audio_spec[i][0];
        float a_im = out_audio_spec[i][1];

        // NEU: Phasenlinearen FIR-EQ anwenden (falls aktiviert und übergeben)
        if (eqCurve != nullptr) {
            a_re *= eqCurve[i];
            a_im *= eqCurve[i];
        }

        float fl_re = filter_spec_low[i][0], fl_im = filter_spec_low[i][1];
        multiplied_spec_low[i][0] = (a_re * fl_re) - (a_im * fl_im);
        multiplied_spec_low[i][1] = (a_re * fl_im) + (a_im * fl_re);

        float fm_re = filter_spec_mid[i][0], fm_im = filter_spec_mid[i][1];
        multiplied_spec_mid[i][0] = (a_re * fm_re) - (a_im * fm_im);
        multiplied_spec_mid[i][1] = (a_re * fm_im) + (a_im * fm_re);

        float fh_re = filter_spec_high[i][0], fh_im = filter_spec_high[i][1];
        multiplied_spec_high[i][0] = (a_re * fh_re) - (a_im * fh_im);
        multiplied_spec_high[i][1] = (a_re * fh_im) + (a_im * fh_re);
    }

    fftwf_execute(plan_inverse_low);
    fftwf_execute(plan_inverse_mid);
    fftwf_execute(plan_inverse_high);

    float norm = 1.0f / N;
    for (int i = 0; i < L; i++) {
        out_low[i]  = final_audio_out_low[M - 1 + i] * norm;
        out_mid[i]  = final_audio_out_mid[M - 1 + i] * norm;
        out_high[i] = final_audio_out_high[M - 1 + i] * norm;
    }
}

// ---------------------------------------------------------
// RACEDspEngine (FIR Hauptklasse)
// ---------------------------------------------------------
RACEDspEngine::RACEDspEngine(double initialDelayUs, double initialAttenuationDb, double initialCenterP, bool initialFreqLimit, int sr) 
    : sampleRate(sr), 
      firCrossoverL(32768, 32767, sr), 
      firCrossoverR(32768, 32767, sr),
      firEqualizer(sr, 65536) // NEU: Initialisierung mit korrekter FFT-Größe (N)
{

    spectrumBuffer.resize(32768, 0.0f);
    spectrumIndex = 0;

    int maxDelaySamples = sr; 
    delayBufferL.resize(maxDelaySamples, 0.0f);
    delayBufferR.resize(maxDelaySamples, 0.0f);
    
    lowL.resize(32768, 0.0f); midL.resize(32768, 0.0f); highL.resize(32768, 0.0f);
    lowR.resize(32768, 0.0f); midR.resize(32768, 0.0f); highR.resize(32768, 0.0f);

    writeIndex = 0;
    volume = 1.0;
    makeupGainDb = 0.0;
    makeupGainLinear = 1.0;
    raceEnabled = true;
    filtersEnabled = true;
    eqEnabled = false;

    currentFreqLow = 400.0;
    currentFreqHigh = 4530.0;
    
    firCrossoverL.updateFilters(currentFreqLow, currentFreqHigh);
    firCrossoverR.updateFilters(currentFreqLow, currentFreqHigh);
    
    // NEU: EQ Bands initialisieren
    for (int i = 0; i < 3; ++i) {
        double f = (i == 0) ? 100.0 : (i == 1) ? 1000.0 : 5000.0;
        eqBands[i] = {f, 0.707, 0.0};
    }
    firEqualizer.reset();

    setParameters(initialDelayUs, initialAttenuationDb, initialCenterP, initialFreqLimit);
}

RACEDspEngine::~RACEDspEngine() {}

void RACEDspEngine::setParameters(double newDelayUs, double newAttenuationDb, double newCenterP, bool newFreqLimit) {
    std::lock_guard<std::mutex> lock(dspMutex);
    delayUs = newDelayUs;
    attenuationDb = newAttenuationDb;
    
    // KORREKTUR 1: Erzwinge mathematisch immer eine Dämpfung, egal ob GUI positiv oder negativ sendet
    double safeAtten = std::abs(attenuationDb);
    attenuationLinear = std::pow(10.0, -safeAtten / 20.0);
    
    centerP = newCenterP;
    freqLimitRACE = newFreqLimit;
    
    dn = static_cast<int>(delayUs * static_cast<double>(sampleRate) / 1000000.0);
    // KORREKTUR 2: RACE benötigt ein minimales Delay, sonst heben sich die Signale im Mittenband komplett auf
    if (dn < 5) dn = 10; 
}

void RACEDspEngine::setVolume(double newVolume) { std::lock_guard<std::mutex> lock(dspMutex); volume = newVolume; }
// neu eingefügte Volumenkorrektur-Funktion U. Thibaut 24.10.2026
    void RACEDspEngine::setMakeupGain(double gainDb) { 
    std::lock_guard<std::mutex> lock(dspMutex); 
    makeupGainDb = gainDb; 
    makeupGainLinear = std::pow(10.0, makeupGainDb / 20.0); 
}
double RACEDspEngine::getMakeupGain() { return makeupGainDb; }
// Ende der neu eingefügten Funktion 


void RACEDspEngine::setRaceEnabled(bool enabled) { std::lock_guard<std::mutex> lock(dspMutex); raceEnabled = enabled; }
void RACEDspEngine::setFiltersEnabled(bool enabled) { std::lock_guard<std::mutex> lock(dspMutex); filtersEnabled = enabled; }
void RACEDspEngine::setEqEnabled(bool enabled) { std::lock_guard<std::mutex> lock(dspMutex); eqEnabled = enabled; }

void RACEDspEngine::setEqBand(int bandIndex, double freq, double q, double gainDb) {
    std::lock_guard<std::mutex> lock(dspMutex);
    if(bandIndex >= 0 && bandIndex < 3) { 
        eqBands[bandIndex] = {freq, q, gainDb};
        
        // Komplette EQ Kurve neu berechnen
        firEqualizer.reset();
        for (int b = 0; b < 3; ++b) {
            if (eqBands[b].g != 0.0) {
                firEqualizer.addPeakingBand(eqBands[b].f, eqBands[b].g, eqBands[b].q);
            }
        }
    }
}

void RACEDspEngine::processSamples(float* interleavedSamples, size_t frames) {
    std::lock_guard<std::mutex> lock(dspMutex);

    if (!filtersEnabled) {
        for (size_t i = 0; i < frames; i++) {
            float outL = interleavedSamples[i * 2] * volume;
            float outR = interleavedSamples[i * 2 + 1] * volume;
            interleavedSamples[i * 2] = outL;
            interleavedSamples[i * 2 + 1] = outR;
            spectrumBuffer[spectrumIndex] = (outL + outR) * 0.5f;
            spectrumIndex = (spectrumIndex + 1) % 32768;
        }
        return;
    }
    
    int L = 32768;
    int bufferSize = delayBufferL.size();

    for (size_t i = 0; i < frames; i++) {
        inputFifoL.push_back(interleavedSamples[i * 2]);
        inputFifoR.push_back(interleavedSamples[i * 2 + 1]);
    }

    while (inputFifoL.size() >= (size_t)L) {
        // NEU: EQ-Kurve nur übergeben, wenn eqEnabled true ist
        const float* eqData = eqEnabled ? firEqualizer.getCurveData() : nullptr;

        firCrossoverL.processBlock(inputFifoL.data(), lowL.data(), midL.data(), highL.data(), eqData);
        firCrossoverR.processBlock(inputFifoR.data(), lowR.data(), midR.data(), highR.data(), eqData);

        inputFifoL.erase(inputFifoL.begin(), inputFifoL.begin() + L);
        inputFifoR.erase(inputFifoR.begin(), inputFifoR.begin() + L);

        for (int k = 0; k < L; k++) {
            float bypassL = 0.0f;
            float bypassR = 0.0f;
            float midInL = 0.0f;
            float midInR = 0.0f;

            // KORREKTUR 3: Restore freqLimitRACE (Schalter im GUI)
            if (freqLimitRACE) {
                // RACE nur im Mittenband
                bypassL = lowL[k] + highL[k];
                bypassR = lowR[k] + highR[k];
                midInL = midL[k];
                midInR = midR[k];
            } else {
                // RACE auf das VOLLSTÄNDIGE Signal anwenden
                bypassL = 0.0f;
                bypassR = 0.0f;
                midInL = lowL[k] + midL[k] + highL[k];
                midInR = lowR[k] + midR[k] + highR[k];
            }

            int readIndex = (writeIndex - dn + bufferSize) % bufferSize;
            float delayedMidL = delayBufferL[readIndex];
            float delayedMidR = delayBufferR[readIndex];

            float crossTalkL = 0.0f;
            float crossTalkR = 0.0f;

            if (raceEnabled) {
                crossTalkL = static_cast<float>(attenuationLinear * (delayedMidR - centerP * delayedMidL) / (1.0 + centerP));
                crossTalkR = static_cast<float>(attenuationLinear * (delayedMidL - centerP * delayedMidR) / (1.0 + centerP));
            }

            float currentRaceMidL = midInL - crossTalkL;
            float currentRaceMidR = midInR - crossTalkR;

            delayBufferL[writeIndex] = currentRaceMidL;
            delayBufferR[writeIndex] = currentRaceMidR;
            writeIndex = (writeIndex + 1) % bufferSize;

            float outL = bypassL + currentRaceMidL;
            float outR = bypassR + currentRaceMidR;

            // ALT: Hier war die if (eqEnabled) Schleife für die Biquads. 
            // Diese ist nun restlos gelöscht!

            // Makeup Gain anwenden, neu eigefügt am 24.10.2026, U. Thibaut

            outL *= static_cast<float>(makeupGainLinear);
            outR *= static_cast<float>(makeupGainLinear);

            outL *= volume;
            outR *= volume;

            // NEU: Clipping-Detektor für das Logfile, z.B. für Langzeit-Tests der Engine! U. Thibaut, 26.10.2026
            if (outL > 1.0f || outL < -1.0f || outR > 1.0f || outR < -1.0f) {
            // Statische Variable, um das Log nicht bei jedem Sample zu fluten (nur 1x pro Puffer)
            static int clipCount = 0; 
            if (clipCount++ % 100 == 0) {
               std::cerr << "[WARN] DSP Clipping erkannt! L: " << outL << ", R: " << outR << std::endl;
              }
            }

            outL = std::clamp(outL, -1.0f, 1.0f);
            outR = std::clamp(outR, -1.0f, 1.0f);

            spectrumBuffer[spectrumIndex] = (outL + outR) * 0.5f;
            spectrumIndex = (spectrumIndex + 1) % 32768;

            outputFifoL.push_back(outL);
            outputFifoR.push_back(outR);
        }
    }

    size_t available = std::min(frames, outputFifoL.size());
    for (size_t i = 0; i < available; i++) {
        interleavedSamples[i * 2] = outputFifoL[i];
        interleavedSamples[i * 2 + 1] = outputFifoR[i];
    }
    for (size_t i = available; i < frames; i++) {
        interleavedSamples[i * 2] = 0.0f;
        interleavedSamples[i * 2 + 1] = 0.0f;
    }
    
    if (available > 0) {
        outputFifoL.erase(outputFifoL.begin(), outputFifoL.begin() + available);
        outputFifoR.erase(outputFifoR.begin(), outputFifoR.begin() + available);
    }
}

std::vector<double> RACEDspEngine::getSpectrumBands() {
    std::vector<float> audioData(32768, 0.0f);
    {
        std::lock_guard<std::mutex> lock(dspMutex);
        for(int i = 0; i < 32768; i++) {
            audioData[i] = spectrumBuffer[(spectrumIndex + i) % 32768];
        }
    }
    
    std::vector<std::complex<float>> complexData(32768);
    for(int i = 0; i < 32768; i++) {
        float multiplier = 0.5f * (1.0f - std::cos(2.0f * static_cast<float>(M_PI) * i / 32767.0f));
        complexData[i] = std::complex<float>(audioData[i] * multiplier, 0.0f);
    }
    
    computeFFT(complexData);
    
    std::vector<double> bands(32, 0.0);
    float minFreq = 20.0f, maxFreq = 20000.0f;
    float logMin = std::log10(minFreq), logMax = std::log10(maxFreq);
    
    for (int b = 0; b < 32; b++) {
        float freqStart = std::pow(10.0f, logMin + (float)b / 32.0f * (logMax - logMin));
        float freqEnd   = std::pow(10.0f, logMin + (float)(b + 1) / 32.0f * (logMax - logMin));
        
        int binStart = (int)(freqStart * 32768.0f / static_cast<float>(sampleRate));
        int binEnd   = (int)(freqEnd * 32768.0f / static_cast<float>(sampleRate));
        
        if (binStart < 1) binStart = 1; 
        if (binEnd > 4095) binEnd = 4095; 
        if (binStart > binEnd) binEnd = binStart;
        
        float maxMag = 0.0f;
        for (int i = binStart; i <= binEnd; i++) {
            float mag = std::abs(complexData[i]) / 16384.0f;
            if (mag > maxMag) maxMag = mag;
        }
        
        float db = 20.0f * std::log10(maxMag + 1e-6f);
        
        // Optische Skalierung für Musik: -85 dB ist der Boden (0.0), -20 dB ist der obere rote Anschlag (1.0), U. Thibaut, 24.10.2026
        float normalized = (db + 85.0f) / 65.0f; 
        
        if (normalized < 0.0f) normalized = 0.0f; 
        if (normalized > 1.0f) normalized = 1.0f;        
        bands[b] = static_cast<double>(normalized);
    }
    return bands;
}

double RACEDspEngine::getDelayUs() { return delayUs; }
double RACEDspEngine::getAttenuationDb() { return attenuationDb; }
double RACEDspEngine::getCenterP() { return centerP; }
bool RACEDspEngine::getFreqLimit() { return freqLimitRACE; }
double RACEDspEngine::getVolume() { return volume; }
bool RACEDspEngine::getRaceEnabled() { return raceEnabled; }
bool RACEDspEngine::getFiltersEnabled() { return filtersEnabled; }
bool RACEDspEngine::getEqEnabled() { return eqEnabled; }

double RACEDspEngine::getEqFreq(int bandIndex) { if (bandIndex >= 0 && bandIndex < 3) return eqBands[bandIndex].f; return 1000.0; }
double RACEDspEngine::getEqQ(int bandIndex) { if (bandIndex >= 0 && bandIndex < 3) return eqBands[bandIndex].q; return 0.707; }
double RACEDspEngine::getEqGain(int bandIndex) { if (bandIndex >= 0 && bandIndex < 3) return eqBands[bandIndex].g; return 0.0; }

void RACEDspEngine::setCrossoverFreqLow(double freq) {
    std::lock_guard<std::mutex> lock(dspMutex);
    currentFreqLow = freq;
    firCrossoverL.updateFilters(currentFreqLow, currentFreqHigh);
    firCrossoverR.updateFilters(currentFreqLow, currentFreqHigh);
}
double RACEDspEngine::getCrossoverFreqLow() { return currentFreqLow; }

void RACEDspEngine::setCrossoverFreqHigh(double freq) {
    std::lock_guard<std::mutex> lock(dspMutex);
    currentFreqHigh = freq;
    firCrossoverL.updateFilters(currentFreqLow, currentFreqHigh);
    firCrossoverR.updateFilters(currentFreqLow, currentFreqHigh);
}
double RACEDspEngine::getCrossoverFreqHigh() { return currentFreqHigh; }
