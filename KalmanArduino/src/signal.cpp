#include "Arduino.h"
#include "signal.h"

double* generateMask(int maskLength, double periodeLength){
    double* mask = new double[maskLength];
    for (int i = 0; i < maskLength; i++){
        mask[i] = sin(double(i * 2 * PI) / periodeLength);
    }

    return mask;
}

//Compute the correlation between two signals with a shift
double correlation(double* mask, int* signal, int length, int shift, int signalIndex){
    double sum = 0;
    double maskMean = 0;
    double signalMean = 0;

    for (int i = 0; i < length; i++){
        maskMean += mask[i];
    }
    maskMean /= length;
    for (int i = 0; i < length; i++){
        signalMean += signal[(i + signalIndex + 1) % SAMPLING_BUFFER_SIZE];
    }
    signalMean /= length;

    for (int i = 0; i < length; i++){
        double maskValue = (double)mask[i] - maskMean;
        int index = (((i + signalIndex + 1) % SAMPLING_BUFFER_SIZE) + shift) % SAMPLING_BUFFER_SIZE;
        double signalValue = (double)signal[index] - signalMean;
        sum += maskValue * signalValue;
    }
    return sum/length;
}

//Compute the correlation between two signals with a shift
double correlationShift(double* mask, int* signal, int length, int max_shift, int signalIndex){
    double max = 0;
    for (int i = 0; i < max_shift; i++){
        double c = correlation(mask, signal, length, i, signalIndex);
        if (c > max){
            max = c;
        }
    }
    return max;
}

void setupADC(ADC* adc){
    adc->adc0->setAveraging(4); // no averaging
    adc->adc0->setResolution(10); // 10 bits resolution
    adc->adc0->setConversionSpeed(ADC_CONVERSION_SPEED::VERY_HIGH_SPEED); // fastest conversion
    adc->adc0->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_HIGH_SPEED); // fastest sampling

    pinMode(A0, INPUT);
    pinMode(A1, INPUT);
    pinMode(A2, INPUT);
}