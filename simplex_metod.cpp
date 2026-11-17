#include <iostream>
using namespace std;

const int M = 7;
const int VARS = 11;

double c[VARS] = {3, 4.5, 5, 6, 0, 0, 0, 0, 0, 0, 0};

void pechat_tablicy(double A[M][VARS], double b[M], int basis[M], double delta[VARS], double Z) {
    cout.precision(2);
    cout << fixed;
    cout << endl << "=========== SIMPLEKS-TABLICA ===========" << endl << endl;
    
    cout << "Bazis\t| ";
    for (int j = 0; j < VARS; j++) {
        if (j < 4) cout << "x" << j + 1 << "\t";
        else cout << "s" << j - 3 << "\t";
    }
    cout << "| b" << endl;
    cout << "\t" << endl;
    
    for (int i = 0; i < M; i++) {
        if (basis[i] < 4) cout << "x" << basis[i] + 1;
        else cout << "s" << basis[i] - 3;
        cout << "\t| ";
        for (int j = 0; j < VARS; j++) {
            cout << A[i][j] << "\t";
        }
        cout << "| " << b[i] << endl;
    }
    cout << "\t" << endl;
    cout << "Dj\t| ";
    for (int j = 0; j < VARS; j++) {
        cout << delta[j] << "\t";
    }
    cout << "| Z = " << Z << endl;
    cout << "=====================================================================" << endl << endl;
    
    cout.unsetf(ios::fixed);
    cout.precision(6);
}

void raschet_delta_z(double A[M][VARS], double b[M], int basis[M], double delta[VARS], double& Z) {
    double cb[M];
    Z = 0;
    for (int i = 0; i < M; i++) cb[i] = c[basis[i]];
    for (int i = 0; i < M; i++) Z += cb[i] * b[i];
    for (int j = 0; j < VARS; j++) {
        double zj = 0;
        for (int i = 0; i < M; i++) zj += cb[i] * A[i][j];
        delta[j] = zj - c[j];
    }
}

void pechat_odr(double A[M][VARS], double b[M], int basis[M]) {
    cout.precision(2);
    cout << fixed;
    cout << endl << "=========== ODR (TEKUSHCHEE SOSTOYANIE) ===========" << endl << endl;
    cout << "Bazis\t| ";
    for (int j = 0; j < VARS; j++) {
        if (j < 4) cout << "x" << j + 1 << "\t";
        else cout << "s" << j - 3 << "\t";
    }
    cout << "| b" << endl;
    cout << "\t" << endl;
    for (int i = 0; i < M; i++) {
        if (basis[i] < 4) cout << "x" << basis[i] + 1;
        else cout << "s" << basis[i] - 3;
        cout << "\t| ";
        for (int j = 0; j < VARS; j++) {
            cout << A[i][j] << "\t";
        }
        cout << "| " << b[i] << endl;
    }
    cout << "=====================================================================" << endl << endl;
    cout.unsetf(ios::fixed);
    cout.precision(6);
}

bool proverka_opornosti(double b[M]) {
    for (int i = 0; i < M; i++) {
        if (b[i] < -1e-9) return false;
    }
    return true;
}

void gauss(double A[M][VARS], double b[M], int basis[M]) {
    cout << endl << "=== GAUSS (nach. bazis) ===" << endl;
    for (int k = 0; k < M; k++) {
        int col = basis[k];
        double pivot = A[k][col];
        if (abs(pivot) < 1e-9) {
            cout << "Oshibka: nulevoy vedushchiy element" << endl;
            return;
        }
        cout << "Vedushchiy element: A[" << k << "][" << col << "] = " << pivot << endl;
        for (int j = 0; j < VARS; j++) A[k][j] /= pivot;
        b[k] /= pivot;
        
        for (int i = 0; i < M; i++) {
            if (i == k) continue;
            double f = A[i][col];
            for (int j = 0; j < VARS; j++) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
        pechat_odr(A, b, basis);
    }
}

bool privesti_k_odr(double A[M][VARS], double b[M], int basis[M]) {
    cout << endl << "=== PRIVEDENIE K ODR ===" << endl;
    while (true) {
        int row = -1;
        double minB = -1e-9;
        for (int i = 0; i < M; i++) {
            if (b[i] < minB) {
                minB = b[i];
                row = i;
            }
        }
        if (row == -1) {
            cout << "ODR DOSTIGNUTA" << endl;
            return true;
        }
        int col = -1;
        for (int j = 0; j < VARS; j++) {
            if (A[row][j] < -1e-9) {
                col = j;
                break;
            }
        }
        if (col == -1) {
            cout << "ODR NE SUSHCESTVUET" << endl;
            return false;
        }
        cout << "Vedushchiy element: row " << row << ", col " << col << endl;
        basis[row] = col;
        double p = A[row][col];
        for (int j = 0; j < VARS; j++) A[row][j] /= p;
        b[row] /= p;
        for (int i = 0; i < M; i++) {
            if (i == row) continue;
            double f = A[i][col];
            for (int j = 0; j < VARS; j++) A[i][j] -= f * A[row][j];
            b[i] -= f * b[row];
        }
        pechat_odr(A, b, basis);
    }
}

int main() {
    double A[M][VARS] = {
        {0.8, 1.5, 2.2, 3,   1, 0, 0, 0, 0, 0, 0},
        {1.2, 2.0, 2.5, 1.8, 0, 1, 0, 0, 0, 0, 0},
        {2.5, 3.0, 1.5, 2.0, 0, 0, 1, 0, 0, 0, 0},
        {-1,  0,   0,   0,   0, 0, 0, 1, 0, 0, 0},
        {0,  -1,   0,   0,   0, 0, 0, 0, 1, 0, 0},
        {0,   0,  -1,   0,   0, 0, 0, 0, 0, 1, 0},
        {0,   0,   0,  -1,   0, 0, 0, 0, 0, 0, 1}
    };
    
    double b[M] = {12000, 90000, 80000, -1800, -1200, -1000, -700};
    int basis[M] = {4, 5, 6, 7, 8, 9, 10};
    
    cout << "Iskhodnye dannye:" << endl;
    cout << "Modeli stulev: standart, komfort, prestizh, lyuks" << endl;
    cout << "Pribyl: 3, 4.5, 5, 6 (rub/sht)" << endl;
    cout << "Obyazatelnye dogovory: 1800, 1200, 1000, 700 sht" << endl;
    cout << "Resursy: 12000 rub, 90000 norm-chas, 80000 m^2" << endl;
    
    gauss(A, b, basis);
    if (!proverka_opornosti(b)) {
        cout << endl << "Est otricatelnye b => privodim k ODR" << endl;
        if (!privesti_k_odr(A, b, basis))
            return 0;
    }
    
    int iteraciya = 0;
    while (true) {
        cout << endl << " ITERACIYA: " << iteraciya << endl;
        double delta[VARS];
        double Z;
        raschet_delta_z(A, b, basis, delta, Z);
        
        int col = -1;
        double min = -1e-9;
        for (int j = 0; j < VARS; j++) {
            if (delta[j] < min) {
                min = delta[j];
                col = j;
            }
        }
        
        if (col == -1) {
            cout << endl << "========== OPTIMUM NAYDEN ==========" << endl;
            pechat_tablicy(A, b, basis, delta, Z);
            cout << endl << "REZULTATY:" << endl;
            cout << "Proizvodstvo stulev:" << endl;
            for (int i = 0; i < M; i++) {
                if (basis[i] < 4) {
                    cout << " x" << basis[i] + 1 << " = " << b[i] << " sht" << endl;
                }
            }
            cout << endl << "Maksimalnaya pribyl Z = " << Z << " rub" << endl;
            break;
        }
        
        if (col < 4) cout << "Vhodit: x" << col + 1 << endl;
        else cout << "Vhodit: s" << col - 3 << endl;
        
        int row = -1;
        double best = 1e18;
        for (int i = 0; i < M; i++) {
            if (A[i][col] > 1e-9) {
                double r = b[i] / A[i][col];
                if (r < best) {
                    best = r;
                    row = i;
                }
            }
        }
        if (row == -1) {
            cout << "ZADACHA NEOGRANICHENA" << endl;
            return 0;
        }
        
        cout << "Vyvoditsya: ";
        if (basis[row] < 4) cout << "x" << basis[row] + 1 << endl;
        else cout << "s" << basis[row] - 3 << endl;
        
        basis[row] = col;
        double p = A[row][col];
        for (int j = 0; j < VARS; j++) A[row][j] /= p;
        b[row] /= p;
        for (int i = 0; i < M; i++) {
            if (i == row) continue;
            double f = A[i][col];
            for (int j = 0; j < VARS; j++) A[i][j] -= f * A[row][j];
            b[i] -= f * b[row];
        }
        
        raschet_delta_z(A, b, basis, delta, Z);
        pechat_tablicy(A, b, basis, delta, Z);
        iteraciya++;
    }
    
    return 0;
}
