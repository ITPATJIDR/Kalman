#define K 4 // Number of reference points
#define FOV 1.5708 // 90 degrees in radians
#define A_BOUND 4 // RX pos in [-A, A, -A, A]

// Rotation matrix
void rot_mat(double theta, double res[3][3]);

// Vector dot product
double dot_product(double a[3], double b[3]);

// Vector norm
double norm(double a[3]);

// Cosine of angle between two vectors
double cos_ang(double a[3], double b[3]);

// Matrix-vector multiplication
void matrix_vector_mult(double mat[3][3], double vec[3], double res[3]);

// Sign function
int sign(double value);

// Derivative of the measurement model with respect to x
double drdx(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m);

// Derivative of the measurement model with respect to y
double drdy(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m);

// Derivative of the measurement model with respect to theta
double drdtheta(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m);

// Check if a reference point is out of the field of view
bool is_out_fov(int k, int l, double x, double y, double p[3][K], double v[3][3], double theta);

// Print matrix
void printMatrix(double* mat, int rows, int cols);

class KalmanFilter{
    private:
        double T = 0.05; // Time step
        double Pp[6][6] = {
            {0.01, 0, 0,    0, 0,    0},
            {0,    1, 0,    0, 0,    0},
            {0,    0, 0.01, 0, 0,    0},
            {0,    0, 0,    1, 0,    0},
            {0,    0, 0,    0, 0.01, 0},
            {0,    0, 0,    0, 0,    1}
        };
        double H0[3][K];
        double r0[3*K];
        double R[3*K][3*K];
        double A[6][6] = {
            {1, T, 0, 0, 0, 0},
            {0, 1, 0, 0, 0, 0},
            {0, 0, 1, T, 0, 0},
            {0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 1, T},
            {0, 0, 0, 0, 0, 1}
            };
        // Process noise covariance
        double Q[6][6] = {
            {0.01, 0, 0, 0, 0, 0},
            {0, 0.01, 0, 0, 0, 0},
            {0, 0, 0.01, 0, 0, 0},
            {0, 0, 0, 0.01, 0, 0},
            {0, 0, 0, 0, 0.01, 0},
            {0, 0, 0, 0, 0, 0.01}
            };
        double p[3][K];
        double v[3][3];
        double h;
        double ka;
        double m;
        double xm[6];
        double Pm[6][6];
        double r[3*K];
        double Hl[3*K][6];
        double Kg[6][3*K];
        double rH[3*K];
        void state_prediction();
        void cov_prediction();
        void linear_H0();
        void kalmanGain();
        void comp_H0(double p0[3]);
        void state_update();
        void cov_update();
        void comp_H0(double p0[3], double theta);
        void comp_R0(double P);
        void beamVectors(double al[3], double be);

    public:
        double xp[6] = {0, 0, 0, 0, 0, 0};
        void compute();
        void updateMeasurements(double r[3*K]);
        void print();
        void setT(double T);

    //Constructor
    KalmanFilter(double p0[3], double theta, double P, double al[3], double be, double ka, double m, double h, double p[3][K]){
        for (int i = 0; i < 3*K; i++) {
            for (int j = 0; j < 3*K; j++) {
                R[i][j] = (i == j) ? 0.29 : 0; //0.29 is the noise variance
            }
        }
        comp_H0(p0, theta);
        comp_R0(P);
        beamVectors(al, be);
        this->ka = ka;
        this->m = m;
        this->h = h;
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < K; j++){
                this->p[i][j] = p[i][j];
            }
        }
    }
};