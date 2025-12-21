package heparansulfate;

import java.util.Random;

/**
 * Geometric programming for solving a MaxEnt problem subject to
 * linear equality constraints (in the probability simplex).
 */
public class MaxEntOptim {
    Random rand = null;
    double[][] A = null;
    double[] b = null;
    int m = 0;
    int n = 0;
    double[] x = null;
    double[][] HZ = null;
    double[][] HG = null;
    double[] gradZ = null;
    double[] gradG = null;
    double Z = 0.;
    double gamma = 0.;
    double[] delta = null;
    public double[] p = null;
    double S = 0.;

    /**
     * @param A matrix in constraint Ap = b
     * @param b vector in constraint Ap = b
     */
    public MaxEntOptim(double[][] A, double[] b) {
        this.A = A;
        this.b = b;
        m = A.length;
        n = A[0].length;
        rand = new Random(2);
        x = new double[m];
        gradZ = new double[m];
        gradG = new double[m];
        HZ = new double[m][m];
        HG = new double[m][m];
        p = new double[n];
        
        for (int i = 0; i < m; i++) {
            x[i] = 1. + rand.nextDouble();
        }

        double lambda = 1.;
        System.err.println("Potential\tNewtonStep\tEntropy\tInfeasibility");
        while (lambda > 1e-9) {
            makeDerivatives();
            makeDelta();
            double alpha = getAlpha();
            for (int i = 0; i < m; i++) {
                x[i] += alpha * delta[i];
            }
            lambda = getLambda();
            makeP();
            System.err.println(getGamma(x) + "\t" + lambda + "\t" + S + "\t" + getFeasibility());
        }
    }

    double getFeasibility() {
        double[] t = MatrixOp.multMatVec(A, p);
        double res = 0.;
        for (int i = 0; i < m; i++) {
            res += Math.abs(t[i] - b[i]);
        }
        return res;
    }

    void makeP() {
        S = 0.;
        double z = getZ(x);
        for (int s = 0; s < n; s++) {
            p[s] = Math.exp(-getDotProd(s, x)) / z;
            if (p[s] > 0) {
                S -= p[s] * Math.log(p[s]);
            }
        }
    }

    double getLambda() {
        double res = MatrixOp.geInnerProd(delta, gradG);
        return Math.abs(res);
    }

    double getAlpha() {
        double s = 0.999;
        double beta = 0.5;
        double sg = 0.001;
        double alpha = s;
        double c = -MatrixOp.geInnerProd(gradG, gradG);
        double f = getGamma(x);
        double f2 = getGamma(x, alpha);
        int iter = 0;
        while (f - f2 < -sg * alpha * c && alpha > 1e-12) {
            alpha *= beta;
            f2 = getGamma(x, alpha);
            iter++;
            if (iter > 1000) {
                System.err.println("Line search iteration limit hit");
                iter = 0;
            }
        }
        return alpha;
    }

    void makeDelta() {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = -gradG[i];
        }
        delta = MatrixOp.choleskySolve(HG, t);
    }

    double getDotProd(int s, double[] y) {
        double res = 0.;
        for (int i = 0; i < m; i++) {
            res += A[i][s] * y[i];
        }
        return res;
    }

    void makeDerivatives() {
        makeZ();
        makeGamma();
        makeGradZ();
        makeGradG();
        makeHZ();
        makeHG();
    }

    void makeHG() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                HG[i][j] = (HZ[i][j] - gradZ[i] * gradZ[j] / Z) / Z;
            }
        }
        double[][] L = MatrixOp.cholesky(HG);
        double min = Math.abs(L[0][0]);
        for (int i = 1; i < L.length; i++) {
            double d = Math.abs(L[i][i]);
            if (d < min) min = d;
        }
        if (min < 0.0001) {
            for (int i = 0; i < m; i++) HG[i][i] += 0.0001;
        }
    }

    void makeHZ() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                HZ[i][j] = 0.;
                for (int s = 0; s < n; s++) {
                    HZ[i][j] += A[i][s] * A[j][s] * Math.exp(-getDotProd(s, x));
                }
            }
        }
    }

    double getZ(double[] y) {
        double res = 0.;
        for (int s = 0; s < n; s++) {
            res += Math.exp(-getDotProd(s, y));
        }
        return res;
    }

    double getGamma(double[] y, double alpha) {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = y[i] + alpha * delta[i];
        }
        return getGamma(t);
    }

    double getGamma(double[] y) {
        double res = Math.log(getZ(y));
        for (int i = 0; i < m; i++) {
            res += y[i] * b[i];
        }
        return res;
    }

    void makeZ() { Z = getZ(x); }
    void makeGamma() { gamma = getGamma(x); }

    void makeGradG() {
        for (int i = 0; i < m; i++) {
            gradG[i] = b[i] + gradZ[i] / Z;
        }
    }

    void makeGradZ() {
        for (int i = 0; i < m; i++) {
            gradZ[i] = 0.;
            for (int s = 0; s < n; s++) {
                gradZ[i] -= A[i][s] * Math.exp(-getDotProd(s, x));
            }
        }
    }
}