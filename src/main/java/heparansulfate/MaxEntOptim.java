package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Geometric programming for solving a Maximum Entropy (MaxEnt) problem subject to
 * linear equality constraints in the probability simplex.
 * <p>
 * This class implements a Newton-based optimization method to find the 
 * probability distribution that maximizes entropy while satisfying linear 
 * constraints of the form {@code Ax = b}.
 */
public class MaxEntOptim {
    /** Random number generator for seeding initial Lagrange multipliers. */
    Random rand = null;
    /** Matrix representing the linear equality constraints {@code Ax = b}. */
    double[][] A = null;
    /** Constant vector representing the targets in the constraints {@code Ax = b}. */
    double[] b = null;
    /** Number of constraints (rows in matrix {@code A}). */
    int m = 0;
    /** Number of variables in the primal MaxEnt problem (columns in matrix {@code A}). */
    int n = 0;
    /** Lagrange multipliers; the numerical solution to the dual problem. */
    double[] x = null;
    /** Hessian of the partition function {@code Z}. */
    double[][] HZ = null;
    /** Hessian of the Legendre potential {@code Gamma}. */
    double[][] HG = null;
    /** Gradient of the partition function {@code Z}. */
    double[] gradZ = null;
    /** Gradient of the Legendre potential {@code Gamma}. */
    double[] gradG = null;
    /** The partition function value. */
    double Z = 0.;
    /** The Legendre potential (dual objective) value. */
    double gamma = 0.;
    /** The calculated descent direction for the current Newton step. */
    double[] delta = null;
    /** The calculated solution to the primal problem (probabilities). */
    double[] p = null;
    /** The calculated Gibbs entropy of the distribution. */
    double S = 0.;

    /**
     * Solves the MaxEnt problem subject to linear equality constraints using Newton's method.
     * 
     * @param A Matrix of the constraint system {@code Ap = b}.
     * @param b Target vector of the constraint system {@code Ap = b}.
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
     * Checks the feasibility of the current primal solution {@code p}.
     * @return The absolute sum of the residual {@code |Ap - b|}.
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
     * Computes the primal probability distribution {@code p} and its entropy {@code S} 
     * based on current dual variables.
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
     * Computes the Newton decrement.
     * @return The absolute value of the inner product of the descent direction and gradient.
     */
    double getLambda() {
        double res = MatrixOp.geInnerProd(delta, gradG);
        return Math.abs(res);
    }

    /**
     * Calculates the line search parameter using the Armijo rule.
     * @return The line search step size {@code alpha}.
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
     * Computes the Newton descent direction by solving the system involving the Hessian.
     */
    void makeDelta() {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = -gradG[i];
        }
        delta = MatrixOp.choleskySolve(HG, t);
    }

    /**
     * Computes the inner product between the {@code s}-th column of matrix {@code A} 
     * and vector {@code y}.
     * @param s Column index of matrix {@code A}.
     * @param y Input vector.
     * @return The calculated dot product.
     */
    double getDotProd(int s, double[] y) {
        double res = 0.;
        for (int i = 0; i < m; i++) {
            res += A[i][s] * y[i];
        }
        return res;
    }

    /**
     * Computes the dual objective function values, gradients, and Hessians.
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
     * Computes the Hessian of the Legendre potential {@code Gamma}.
     */
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
            if (d < min) {
                min = d;
            }
        }
        if (min < 0.0001) { // Tikhonov regularization for ill-conditioned Hessians
            for (int i = 0; i < m; i++) {
                HG[i][i] += (min * min);
            }
        }
    }

    /**
     * Computes the Hessian of the partition function {@code Z}.
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
     * Computes the partition function {@code Z} for a given set of multipliers.
     * @param y Lagrange multipliers.
     * @return The partition function value.
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
     * Evaluates the Legendre potential {@code Gamma} at a stepped position.
     * @param y Base Lagrange multipliers.
     * @param alpha Step size.
     * @return The potential at {@code y + alpha * delta}.
     */
    double getGamma(double[] y, double alpha) {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = y[i] + alpha * delta[i];
        }
        return getGamma(t);
    }

    /**
     * Evaluates the Legendre potential {@code Gamma}.
     * @param y Lagrange multipliers.
     * @return The Legendre potential value.
     */
    double getGamma(double[] y) {
        double res = Math.log(getZ(y));
        for (int i = 0; i < m; i++) {
            res += y[i] * b[i];
        }
        return res;
    }

    /** Updates the current partition function {@code Z}. */
    void makeZ() {
        Z = getZ(x);
    }

    /** Updates the current dual objective {@code Gamma}. */
    void makeGamma() {
        gamma = getGamma(x);
    }

    /** Updates the current gradient of the Legendre potential {@code Gamma}. */
    void makeGradG() {
        for (int i = 0; i < m; i++) {
            gradG[i] = b[i] + gradZ[i] / Z;
        }
    }

    /** Updates the current gradient of the partition function {@code Z}. */
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
     * Generates default constraints for the MaxEnt problem.
     * @param inDir Input directory for abundance and specificity files.
     * @return A merged set of linear constraints.
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