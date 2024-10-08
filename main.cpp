#include <ADC.h>
#include <Arduino.h>
#include <math.h>
#include <MatrixMath.h>

void sampleSignal();
double cosAngle(double a[3], double b[3]);
double sgn(double x);

ADC *adc = new ADC();
IntervalTimer sampleTimer;

const int WS=480; //window size (sample)
const int MS=400; //mask size (sample)
const int MR=80; //maximization range (sample)
const int SP=25; //sampling period (us)

double h=0.35;
double D=0.415;
double P=1.0;
double m=1.0;
double FoV=1.5708;
double ka=2.178; //ka=2.178
double sigma2=0.29; //0.29; 
double p[4][3] = {{D,D,h}, 
                  {D,0,h}, 
                  {0,0,h}, 
                  {0,D,h}};
double al[3] = {0, 120*PI/180, -120*PI/180};
double be=60*PI/180;
double v[3][3];

int sampleIndex=0;
int sig[3][WS];
double mask[4][MS];
double rss[4][3];
double r[12];
double rHat[12];

unsigned long Tms=50; //EKF iteration period (ms)
double T=Tms*0.001;
double xp[6] = {0,0,0,0,0,0};
double Pp[6][6] = {{0.01,0,0,0,0,0},
                   {0,1,0,0,0,0},
                   {0,0,0.01,0,0,0},
                   {0,0,0,1,0,0},
                   {0,0,0,0,0.01,0},
                   {0,0,0,0,0,1}};
double xm[6];
double Pm[6][6];
double A[6][6] = {{1,T,0,0,0,0},
                  {0,1,0,0,0,0},
                  {0,0,1,T,0,0},
                  {0,0,0,1,0,0},
                  {0,0,0,0,1,T},
                  {0,0,0,0,0,1}};
double Q[6][6];
double R[12][12];
double Hl[12][6];
double Kgain[6][12];


void setup() {
  int i, j;

  Serial.begin(115200); 
  delay(200); 
  Serial.println("Started");

  adc->adc0->setAveraging(0);
  adc->adc0->setResolution(10);
  adc->adc0->setConversionSpeed(ADC_CONVERSION_SPEED::VERY_HIGH_SPEED);
  adc->adc0->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_HIGH_SPEED);
  pinMode(A0,INPUT); 
  pinMode(A1,INPUT); 
  pinMode(A2,INPUT);

  for (i=0; i<MS; i++) {
    mask[0][i] = sin(2*PI*i*3/80); //1.5 kHz
    mask[1][i] = sin(2*PI*i/40); //1 kHz
    mask[2][i] = sin(2*PI*i/80); //0.5 kHz
    mask[3][i] = sin(2*PI*i/20); //2 kHz
  } 

  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) 
      Q[i][j] = (i==j) ? 0.01:0;

  for (i=0; i<12; i++) 
    for (j=0; j<12; j++) 
      R[i][j] = (i==j) ? sigma2:0;

  for (j=0; j<3; j++) { 
    v[j][0] = cos(be)*cos(al[j]); 
    v[j][1] = -cos(be)*sin(al[j]); 
    v[j][2] = sin(be);
  }
}


void loop() {
  unsigned long t1, t2;
  int i, j, k, l;
  double sigMean[3];
  double tmp0;
  double tmp[6][6];
  double tmp2[6][12];
  double tmp3[12][12];
  double tmp4[6][6];
  double tmp5[6][6];
  double tmp6[3];
  double AA, BB, CC, DD;
  double u[3];
  double u2[3] = {0,0,-1};
  double d, cosPhi, cosPsi;

  //delay(500);
  t1 = millis();
  sampleTimer.end();
  
  for (i=0; i<4; i++) 
    for (j=0; j<3; j++) { 
      rss[i][j] = 0;
      for (k=0; k<MR; k++) { 
        sigMean[j] = 0;
        for (l=0; l<MS; l++) 
          sigMean[j] += 1.0/MS * sig[j][(sampleIndex + k + l)%WS]; 
        tmp0 = 0;
        for (l=0; l<MS; l++) 
          tmp0 += mask[i][l] * (sig[j][(sampleIndex + k + l)%WS] - sigMean[j]);
        if (tmp0 > rss[i][j]) 
          rss[i][j] = tmp0; 
      }
      rss[i][j] *= 1.0/MS;
      if (j==0)
        rss[i][j] *= 1.05;
    }
  
  for (i=0; i<4; i++) 
    for (j=0; j<3; j++) {
      r[3*i+j] = rss[i][j];
      //Serial.print(r[3*i+j]); Serial.print(",");
    }
  //Serial.println("");

  // EKF ITERATION: BEGIN
  // xm = A*xp
  for (i=0; i<6; i++) { 
    xm[i] = 0;
    for (j=0; j<6; j++) 
      xm[i] += A[i][j]*xp[j]; 
  }

  // Pm = A*Pp*A' + Q
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) { 
      tmp[i][j] = 0;
      for (k=0; k<6; k++) 
        tmp[i][j] += A[i][k]*Pp[k][j]; 
    }
  
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) { 
      Pm[i][j] = Q[i][j];
      for (k=0; k<6; k++) 
        Pm[i][j] += tmp[i][k]*A[j][k]; 
    }
  
  //SUBROUTINE: COMPUTE Hl
  for (i=0; i<4; i++) 
    for (j=0; j<3; j++) { 
      for (k=0; k<6; k++) 
        Hl[3*i+j][k] = 0;
      AA = v[j][0]*cos(xm[4]) - v[j][1]*sin(xm[4]);
      BB = v[j][0]*sin(xm[4]) + v[j][1]*cos(xm[4]);
      CC = pow(p[i][0]-xm[0],2) + pow(p[i][1]-xm[2],2) + pow(h,2); 
      DD = AA*(p[i][0]-xm[0]) + BB*(p[i][1]-xm[2]) + h*v[j][2];
      u[0] = v[j][0]*cos(xm[4]) - v[j][1]*sin(xm[4]);
      u[1] = v[j][0]*sin(xm[4]) + v[j][1]*cos(xm[4]);
      u[2] = v[j][2];
      tmp6[0] = p[i][0] - xm[0]; 
      tmp6[1] = p[i][1] - xm[2]; 
      tmp6[2] = p[i][2];
      cosPsi = cosAngle(tmp6, u);
      if (cosPsi > cos(FoV)) {
        Hl[3*i+j][0] = -ka*(m+1)*pow(h,m)*sgn(DD)/pow(CC,(m+3)/2)*AA
          + ka*(m+1)*(m+3)*pow(h,m)*fabs(DD)/pow(CC,(m+5)/2)*(p[i][0]-xm[0]);
        Hl[3*i+j][2] = -ka*(m+1)*pow(h,m)*sgn(DD)/pow(CC,(m+3)/2)*BB
          + ka*(m+1)*(m+3)*pow(h,m)*fabs(DD)/pow(CC,(m+5)/2)*(p[i][1]-xm[2]);
        Hl[3*i+j][4] = ka*(m+1)*pow(h,m)*sgn(DD)/pow(CC,(m+3)/2) 
          * ( -(p[i][0]-xm[0]) * (v[j][0]*sin(xm[4]) + v[j][1]*cos(xm[4])) 
          + (p[i][1]-xm[2]) * (v[j][0]*cos(xm[4]) - v[j][1]*sin(xm[4])) );
    }
  }
  
  // Kgain = Pm * Hl' * inv(R + Hl*Pm*Hl') 
  for (i=0; i<6; i++) 
    for (j=0; j<12; j++) { 
      tmp2[i][j] = 0;
      for (k=0; k<6; k++) 
        tmp2[i][j] += Pm[i][k]*Hl[j][k]; 
    }
  
  for (i=0; i<12; i++) 
    for (j=0; j<12; j++) { 
      tmp3[i][j] = R[i][j];
      for (k=0; k<6; k++) 
        tmp3[i][j] += Hl[i][k]*tmp2[k][j]; 
    }

  Matrix.Invert((double*) tmp3, 12);

  for (i=0; i<6; i++) 
    for (j=0; j<12; j++) { 
      Kgain[i][j] = 0;
      for (k=0; k<12; k++) 
        Kgain[i][j] += tmp2[i][k]*tmp3[k][j]; 
    }
  
  // SUBROUTINE: COMPUTE rHat
  for (j=0; j<3; j++) {
    u[0] = v[j][0]*cos(xm[4]) - v[j][1]*sin(xm[4]);
    u[1] = v[j][0]*sin(xm[4]) + v[j][1]*cos(xm[4]);
    u[2] = v[j][2];
    for (i=0; i<4; i++) {
      tmp6[0] = xm[0]-p[i][0]; 
      tmp6[1] = xm[2]-p[i][1]; 
      tmp6[2] = -p[i][2];
      d = sqrt(tmp6[0]*tmp6[0] + tmp6[1]*tmp6[1] + tmp6[2]*tmp6[2]);
      cosPhi = cosAngle(tmp6, u2);
      tmp6[0] = -tmp6[0]; 
      tmp6[1] = -tmp6[1]; 
      tmp6[2] = -tmp6[2];
      cosPsi = cosAngle(tmp6, u);
      if (cosPsi > cos(FoV))
        rHat[3*i+j] = fabs( P*ka*(m+1) / pow(d,2) * pow(cosPhi,m) * cosPsi );
      else 
        rHat[3*i+j] = 0;
    }
  } 

  // xp = xm + Kgain*(r - rHat)
  for (i=0; i<6; i++) { 
    xp[i] = xm[i];
    for (j=0; j<12; j++) 
      xp[i] += Kgain[i][j] * (r[j]-rHat[j]); 
  }

  // Pp = (eye(6) - Kgain*Hl)*Pm*(eye(6) - Kgain*Hl)' + Kgain*R*Kgain'
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) { 
      tmp[i][j] = (i==j) ? 1:0;
      for (k=0; k<12; k++) 
        tmp[i][j] -= Kgain[i][k]*Hl[k][j]; 
    }
  
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) { 
      tmp4[i][j] = 0;
      for (k=0; k<6; k++) 
        tmp4[i][j] += tmp[i][k]*Pm[k][j]; 
    }
  
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) { 
      tmp5[i][j] = 0;
      for (k=0; k<6; k++) 
        tmp5[i][j] += tmp4[i][k]*tmp[j][k]; 
    }

  for (i=0; i<6; i++) 
    for (j=0; j<12; j++) { 
      tmp2[i][j] = 0;
      for (k=0; k<12; k++) 
        tmp2[i][j] += Kgain[i][k]*R[k][j]; 
    }
  
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) { 
      tmp[i][j] = 0;
      for (k=0; k<12; k++) 
        tmp[i][j] += tmp2[i][k]*Kgain[j][k]; 
    }
  
  for (i=0; i<6; i++) 
    for (j=0; j<6; j++) 
      Pp[i][j] = tmp5[i][j] + tmp[i][j];

  // EKF ITERATION: END

  Serial.print(xp[0]); Serial.print(" "); 
  Serial.print(xp[2]); Serial.print(" "); 
  Serial.print(xp[4]); Serial.println("");
  sampleTimer.begin(sampleSignal, SP);

  t2 = millis();
  if ((t2-t1) < Tms)
    delay(Tms - (t2-t1)); 
  else 
    Serial.println("computation error");
}


void sampleSignal() {
  sig[0][sampleIndex] = adc->adc0->analogRead(A0); 
  sig[1][sampleIndex] = adc->adc0->analogRead(A1); 
  sig[2][sampleIndex] = adc->adc0->analogRead(A2); 
  sampleIndex += 1;
  if (sampleIndex >= WS) 
    sampleIndex = 0;
}


double cosAngle(double a[3], double b[3]) {
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) 
         / sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) 
         / sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
}


double sgn(double x) {
  if (x > 0) 
    return 1; 
  else if (x < 0) 
    return -1; 
  else 
    return 0;
}