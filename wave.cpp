#include <iostream>
#include <windows.h>
#include <vector>
#include <fstream>
#include "shyfft.h"
#include <cmath>
#pragma warning(disable:4146)

const int sampleRate = 44100;
const int waveBufferSize = 512;
const float frequency = 689.968;
const int intputBufferSize = 16;
const double M_PI = 3.14159265358979323846;

const float loTreshold = 10.0;
const float hiTreshold = 0.1;
const float loHiRatio = 0.1;

float gReadPointer = 0;
float waveBuffer[waveBufferSize];
std::vector<float> rawInput;
std::vector<float> circularInput;
bool result = 0;

// FFT
const int fftSize = 1024;
const int hopSize = 256;
const int bufferSize = 4096;
std::vector<float> lastInputPhases(fftSize);
std::vector<float> analysisMagnitudes(fftSize / 2 + 1);
std::vector<float> analysisFrequencies(fftSize / 2 + 1);
std::vector<float> inputCircularBuffer(bufferSize);
std::vector<float> outputCircularBuffer(bufferSize);
int inputPointer = 0;
int outputReadPointer = 0;
int outputWritePointer = hopSize;
int hopCounter = 0;


typedef stmlib::ShyFFT<float, fftSize> ShyFFT;
ShyFFT shyFFT;

void saveNumbersToFile(const std::vector<float>& numbers, const std::string& filename) 
{
    std::ofstream file(filename);

    if (file.is_open()) {
        for (float number : numbers) {
            file << number << ", ";
        }
        file.close();
        std::cout << "Liczby zostaly zapisane do pliku " << filename << std::endl;
    }
    else {
        std::cout << "Nie mozna otworzyc pliku do zapisu." << std::endl;
    }
}

void print()
{
    for (unsigned int n = 0; n < waveBufferSize; n++) {
        std::cout << waveBuffer[n] << ", ";
    };
    std::cout << "\n";

}

bool setup()
{
    // SIN
    for (unsigned int n = 0; n < waveBufferSize; n++) {
        waveBuffer[n] = 0.25 * sin(2.0 * 3.1415 * (float)n / (float)waveBufferSize); // +
            //0.125 * sin(2.0 * 2.0 * 3.1415 * (float)n / (float)waveBufferSize) +
            //0.125 * sin(12.0 * 2.0 * 3.1415 * (float)n / (float)waveBufferSize) + 
            //0.25 * sin(15.0 * 2.0 * 3.1415 * (float)n / (float)waveBufferSize);
    }

    // SAW
   /* for (unsigned int n = 0; n < waveBufferSize / 2; n++) {
        waveBuffer[n] = -1.0 + 4.0 * (float)n / (float)waveBufferSize;
    }
    for (unsigned int n = waveBufferSize / 2; n < waveBufferSize; n++) {
        waveBuffer[n] = -1.0 + 4.0 * (float)(n - waveBufferSize / 2) / (float)waveBufferSize;
    }*/

    // SQUARE
   /* for (unsigned int n = 0; n < waveBufferSize / 2; n++) {
        waveBuffer[n] = -1.0;
    }
    for (unsigned int n = waveBufferSize / 2; n < waveBufferSize; n++) {
        waveBuffer[n] = 1.0;
    }*/

    return true;
}

float wrapPhase(float phaseIn)
{
    if (phaseIn >= 0)
        return fmodf(phaseIn + M_PI, 2.0 * M_PI) - M_PI;
    else
        return fmodf(phaseIn - M_PI, -2.0 * M_PI) + M_PI;
}

int counter = 0;

bool processFFT() 
{
    std::vector<float> unwrappedBuffer(fftSize);
    std::vector<float> unwrappedBufferCopy(fftSize);
    int maxBinIndex = 0;
    float maxBinValue = 0;


    for (int n = 0; n < fftSize; n++) {
        int unwrappedBufferIndex = (inputPointer - fftSize + n + bufferSize) % bufferSize;
        unwrappedBuffer[n] = inputCircularBuffer[unwrappedBufferIndex];
        unwrappedBufferCopy[n] = inputCircularBuffer[unwrappedBufferIndex];
    };

    std::vector<float> bins(fftSize);
    std::vector<float> output(fftSize);
    std::vector<float> magnitudes(fftSize / 2);

    shyFFT.Direct(unwrappedBufferCopy.data(), bins.data()); // input = output repeated 2 times o_0

    //fft.fft(input.data(), re.data(), im.data());
    for (int i = 0; i < fftSize / 2; ++i)
    {
        magnitudes[i] = 2 * sqrt(bins[i] * bins[i] + bins[fftSize / 2 + i] * bins[fftSize / 2 + i]);
        float phase = atan2(bins[i], bins[fftSize / 2 + i]);

        float phaseDiff = phase - lastInputPhases[i];

        float bitCentreFrequency = 2.0 * M_PI * (float)i / (float)fftSize;
        phaseDiff = wrapPhase(phaseDiff - bitCentreFrequency * hopSize);

        float binDeviation = phaseDiff * (float)fftSize / hopSize / (2.0 * M_PI);

        analysisFrequencies[i] = (float)i + binDeviation;
        analysisMagnitudes[i] = magnitudes[i];
        lastInputPhases[i] = phase;

        if (magnitudes[i] > maxBinValue) {
            maxBinValue = magnitudes[i];
            maxBinIndex = i;
        }
    }

    shyFFT.Inverse(bins.data(), output.data()); // output = input, bins = shit

    for (int n = 0; n < fftSize; n++) {
        int circularBufferIndex = (outputWritePointer + n) % bufferSize;
        outputCircularBuffer[circularBufferIndex] += unwrappedBuffer[n];
    }

    if (counter++ > 20) {
        std::cout << maxBinIndex << ", " << analysisFrequencies[maxBinIndex] * sampleRate / fftSize << std::endl;
        saveNumbersToFile(rawInput, "raw_input.txt");
        saveNumbersToFile(circularInput, "cir_input.txt");
        saveNumbersToFile(magnitudes, "bins.txt");
        saveNumbersToFile(output, "output.txt");
        return 1;
    }

    return 0;
}

void render()
{
    for (unsigned int n = 0; n < intputBufferSize; n++) {
        int indexBelow = (int)gReadPointer;
        int indexAbove = indexBelow + 1 >= waveBufferSize ? 0 : indexBelow + 1;
        float fractionAbove = gReadPointer - indexBelow;
        float fractionBelow = 1 - fractionAbove;

        gReadPointer += waveBufferSize * frequency / sampleRate;

        if (gReadPointer >= waveBufferSize)
            gReadPointer -= waveBufferSize;

        float in = waveBuffer[indexBelow] * fractionBelow
            + waveBuffer[indexAbove] * fractionAbove;

        // input -> circular input
        inputCircularBuffer[inputPointer++] = in;
        if (inputPointer >= bufferSize) {
            inputPointer = 0;
        };

        // output <- circular output
        float out = outputCircularBuffer[outputReadPointer];
        outputCircularBuffer[outputReadPointer] = 0;

        out *= (float)hopSize / (float)fftSize;

        outputReadPointer++;
        if (outputReadPointer >= bufferSize) {
            outputReadPointer = 0;
        };

        if (++hopCounter >= hopSize) {
            result = processFFT();
            hopCounter = 0;
            outputWritePointer = (outputWritePointer + hopSize) % bufferSize;
        }

        //

        rawInput.push_back(in);
        circularInput.push_back(out);
        
    }
}

int main() {
    shyFFT.Init();

    setup();
    rawInput.reserve(300 * bufferSize);
    circularInput.reserve(300 * bufferSize);
    //print();

    unsigned int i = 0;
    while (true)
    {
        render();
        if (result == true) return 0;
    }

    return 0;
}