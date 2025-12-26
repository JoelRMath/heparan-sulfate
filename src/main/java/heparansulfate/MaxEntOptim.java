package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Geometric programming for solving a MaxEnt problem subject to
 * linear equality constraints (in the probability simplex)
 */
public class MaxEntOptim {
    /**
     * used for seeding multipliers
     */
    Random rand = null;
    /**
     * matrix in constraints Ax = b
     */
    double[][] A = null;
    /**
     * vector in constraints Ax = b
     */
    double[] b = null;
    /**
     * number of constraints
     */
    int m = 0;
    /**
     * number of variables in primal MaxEnt problem
     */
    int n = 0;
    /**
     * Lagrange multipliers, solution to the dual
     */
    double[] x = null;
    /**
     * Hessian of partition function Z
     */
    double[][] HZ = null;
    /**
     * Hessian of Legendre potential Gamma
     */
    double[][] HG = null;
    /**
     * gradient of partition function Z
     */
    double[] gradZ = null;
    /**
     * gradient of Legendre potential Gamma
     */
    double[] gradG = null;
    /**
     * partition function
     */
    double Z = 0.;
    /**
     * Legendre potential
     */
    double gamma = 0.;
    /**
     * descent direction
     */
    double[] delta = null;
    /**
     * solution to the primal
     */
    double[] p = null;
    /**
     * entropy
     */
    double S = 0.;

    /**
     * Geometric programming for solving a MaxEnt problem subject to linear
     * equality constraints (in the probability simplex)
     * @param A matrix in constraint Ap = b
     * @param b vector in constraint Ap = b
     */
    public MaxEntOptim(double[][] A, double[] b) {
        this.A = A;
        this.b = b;
        m = A.length;
        n = A[0].length;
        rand = new Random(1);
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

    /**
     * checks feasibility of p
     * @return feasibility of p
     */
    double getFeasibility() {
        double[] t = MatrixOp.multMatVec(A, p);
        double res = 0.;
        for (int i = 0; i < m; i++) {
            res += Math.abs(t[i] - b[i]);
        }
        return res;
    }

    /**
     * computes primal solution and entropy
     */
    void makeP() {
        S = 0.;
        double z = getZ(x);
        for (int s = 0; s < n; s++) {
            p[s] = Math.exp(-getDotProd(s, x)) / z;
            S -= p[s] * Math.log(p[s]);
        }
    }

    /**
     * Newton step
     * @return Newton step
     */
    double getLambda() {
        double res = MatrixOp.geInnerProd(delta, gradG);
        return Math.abs(res);
    }

    /**
     * line search parameter (Armijo rule)
     * @return line search parameter (Armijo rule)
     */
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
                System.err.println("iter");
                iter = 0;
            }
        }
        return alpha;
    }

    /**
     * computes the descent direction
     */
    void makeDelta() {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = -gradG[i];
        }
        delta = MatrixOp.choleskySolve(HG, t);
    }

    /**
     * convenience method for computing ((As)^T)y
     * @param s column index
     * @param y vector of dim m
     * @return ((As)^T)y
     */
    double getDotProd(int s, double[] y) {
        double res = 0.;
        for (int i = 0; i < m; i++) {
            res += A[i][s] * y[i];
        }
        return res;
    }

    /**
     * computes gradient and Hessian
     */
    void makeDerivatives() {
        makeZ();
        makeGamma();
        makeGradZ();
        makeGradG();
        makeHZ();
        makeHG();
    }

    /**
     * Hessian of gamma
     */
    void makeHG() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                HG[i][j] = (HZ[i][j] - gradZ[i] * gradZ[j] / Z) / Z;
            }
        }
        double[][] L = MatrixOp.cholesky(HG);
        double min = Math.abs(L[0][0]);
        double max = Math.abs(L[0][0]);
        double smin = L[0][0];
        for (int i = 1; i < L.length; i++) {
            double d = Math.abs(L[i][i]);
            if (d < min) {
                min = d;
            }
            if (d > max) {
                max = d;
            }
            if (L[i][i] < smin) {
                smin = L[i][i];
            }
        }
        System.out.println(min + "\t" + smin);
        if (min < 0.0001) { // approximate detection of ill-conditioned problems
            for (int i = 0; i < m; i++) {
                HG[i][i] += (min * min);
            }
        }
    }

    /**
     * Hessian of Z
     */
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

    /**
     * partition function Z
     * @param y multipliers
     * @return partition function Z
     */
    double getZ(double[] y) {
        double res = 0.;
        for (int s = 0; s < n; s++) {
            double t = getDotProd(s, y);
            res += Math.exp(-t);
        }
        return res;
    }

    /**
     * Legendre potential Gamma
     * @param y multipliers
     * @param alpha line search parameter
     * @return Legendre potential Gamma
     */
    double getGamma(double[] y, double alpha) {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = y[i] + alpha * delta[i];
        }
        return getGamma(t);
    }

    /**
     * Legendre potential Gamma
     * @param y multipliers
     * @return Legendre potential Gamma
     */
    double getGamma(double[] y) {
        double res = Math.log(getZ(y));
        for (int i = 0; i < m; i++) {
            res += y[i] * b[i];
        }
        return res;
    }

    /**
     * partition function Z
     */
    void makeZ() {
        Z = getZ(x);
    }

    /**
     * Legendre potential Gamma
     */
    void makeGamma() {
        gamma = getGamma(x);
    }

    /**
     * gradient of Gamma
     */
    void makeGradG() {
        for (int i = 0; i < m; i++) {
            gradG[i] = b[i] + gradZ[i] / Z;
        }
    }

    /**
     * gradient of Z
     */
    void makeGradZ() {
        for (int i = 0; i < m; i++) {
            gradZ[i] = 0.;
            for (int s = 0; s < n; s++) {
                double t = getDotProd(s, x);
                gradZ[i] -= A[i][s] * Math.exp(-t);
            }
        }
    }

    /**
     * default linear equality constraint = composition-1 + 2(digest-1)
     * @param inDir input directory
     * @return default linear equality constraint = composition-1 + 2(digest-1)
     */
    public static LinEqCons getDefaultLEC(String inDir) {
        Species sp = new Species(2, 16);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        List<LinEqCons> v = new ArrayList<>();
        v.add(LinEqCons.removeLastRow(sp.getCompLEC(bbs.rho)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
                new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
                new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
        return new LinEqCons(v);
    }
}