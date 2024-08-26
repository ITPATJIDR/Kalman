//Main file for Visual Light Positioning System
#include <math.h>
#include <Arduino.h>
#include <ADC.h>
#include "kalmanFilter.h"
#include "signal.h"

void readInputs();
void printCorrelations();
void printCorrelations(int led);
void printCorrelations(int led, int diode);
void printSignals();

// System parameters ***************************************************************************************************
double h = 0.35;    // Height of the receiver
double D = 0.415;    // Distance between reference points
double p[3][K] = { // Reference point positions
    {0, D, 0, D},
    {0, 0, D, D},
    {h, h, h, h}
};

double be = 60 * PI / 180; // Beam elevation angle
double al[3] = {0, 120 * PI / 180, -120 * PI / 180}; // Beam azimuth angles

double p0[3] = {0,0,0}; // Initial position guess
double theta0 = 0; // Initial angle guess

double ka = 2.178; // Placeholder for constant 'ka'
double m = 1; // Placeholder for constant 'm'
double P = 1; // Placeholder for constant 'P'


int signalFrequencies[K] = {500, 1000, 2000, 1500};
int samplingFrequency = 40000; //Sampling frequency

double correlationCompensationFactor[K][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

//Sampling variables ***************************************************************************************************
int signal[3][SAMPLING_BUFFER_SIZE]; //Buffer for the signal
IntervalTimer timerSampling; //Timer for sampling
ADC *adc = new ADC(); // adc object
int indexSample = 0; //Index of the current sample in the rolling buffer

//Correlation variables ************************************************************************************************
double *maskSignal[K];
int maskLength[K];
double correlations[K][3];
unsigned long last_millis_correlation = 0;

//Kalman filter ********************************************************************************************************
KalmanFilter kF = KalmanFilter(p0, theta0, P, al, be, ka, m, h, p);
unsigned long last_millis_ekf = 0;
unsigned long EKF_period = 50; //Period of the EKF in ms

void setup() {
    Serial.begin(115200);

    //Setting up ADC for real fast reading
    setupADC(adc);

    //Generate correlation masks
    for(int i = 0; i < K; i++){
        double periodeSize = (double)samplingFrequency/(double)signalFrequencies[i];
        maskLength[i] = SAMPLING_BUFFER_SIZE - MAXIMUM_SHIFT;
        maskSignal[i] = generateMask(maskLength[i], periodeSize);
    }



    //Setting up the sampling timer
    timerSampling.begin(readInputs, 1000000/samplingFrequency);
    timerSampling.priority(0);
}

void loop() {
 
    Serial.print("Current signal[0] value: ");
    Serial.println(signal[0][indexSample]);
    if (millis() - last_millis_ekf > EKF_period) {

        timerSampling.end(); //Pause sampling interrupt

        kF.setT((millis() - last_millis_ekf)/1000.0); //Update time step for the EKF
        last_millis_ekf = millis(); //Update last time

        //Correlation ***************************************************************************************************
        //(i is the signal, j is the sensor)
        for(int i = 0; i < K; i++){
            for(int j = 0; j < 3; j++){
                //Print signal and mask
                correlations[i][j] = correlationShift(maskSignal[i], signal[j], maskLength[i], MAXIMUM_SHIFT, indexSample);
                correlations[i][j] *= correlationCompensationFactor[i][j];
            }
        }
        //printCorrelations(0);
        
        //EKF **********************************************************************************************************

        double r[3*K] = {correlations[0][0], correlations[0][1], correlations[0][2],
                         correlations[1][0], correlations[1][1], correlations[1][2],
                         correlations[2][0], correlations[2][1], correlations[2][2],
                         correlations[3][0], correlations[3][1], correlations[3][2]};
        kF.updateMeasurements(r);
        kF.compute();
        kF.print();

        //End EKF ******************************************************************************************************
        
        timerSampling.begin(readInputs, 1000000/samplingFrequency); //Start sampling again
    }
}

void readInputs(){
    //Shift Index
    indexSample++;
    if(indexSample >= SAMPLING_BUFFER_SIZE - 1) indexSample = 0;

    //Read the new values
    signal[0][indexSample] = adc->adc0->analogRead(A0);
    signal[1][indexSample] = adc->adc0->analogRead(A1);
    signal[2][indexSample] = adc->adc0->analogRead(A2);
}

void printCorrelations(){
    for(int i = 0; i < K; i++){
        for(int j = 0; j < 3; j++){
            Serial.print(correlations[i][j]);
            Serial.print(" ");
        }
    }
    Serial.println();
}

void printCorrelations(int led){
    for(int j = 0; j < 3; j++){
        Serial.print(correlations[led][j]);
        Serial.print(" ");
    }
    Serial.println();
}

void printCorrelations(int led, int diode){
    Serial.println(correlations[led][diode]);
}

void printSignals(){
    for(int i = 0; i < SAMPLING_BUFFER_SIZE; i++){
        for(int j = 0; j < 3; j++){
            Serial.print(signal[j][i]);
            Serial.print(" ");
        }
        Serial.println();
    }
}