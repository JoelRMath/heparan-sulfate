package heparansulfate;

import java.util.Random;

/**
 * Linear inequality constraints on a H&amp;C model (homogeneous Markov model).
 * Includes constraints for nonnegativity, balance equations, and the stochastic matrix property.
 */
public class HCModelLIC {
    /**
     * Matrix in {@code Ax <= b}.
     */
    double[][] A = null;
    /**
     * Vector in {@code Ax <= b}.
     */
    double[] b = null;
    /**
     * Number of building blocks.
     */
    int m = 0;

    /**
     * Linear inequality constraints on a H&amp;C model (homogeneous Markov model).
     * Nonnegativity, balance equation and stochastic matrix.
     * @param bbs Set of building blocks.
     */
    public HCModelLIC(BBSet bbs) {
        m = bbs.m;
        A = new double[m * (m + 4)][m * m];
        b = new double[m * (m + 4)];
        // Nonnegativity
        for (int rr = 1; rr <= m * m; rr++) {
            int r = rr - 1;
            A[r][r] = -1.;
            b[r] = 0.;
        }
        // Balance equations: <= rho
        for (int rr = m * m + 1; rr <= m * m + m; rr++) {
            int r = rr - 1;
            b[r] = bbs.rho[r - m * m];
            for (int kk = 1; kk <= m * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == rr - m * m) {
                    A[r][k] = bbs.rho[getII(kk) - 1];
                }
            }
        }
        // Balance equations: <= -rho
        for (int rr = m * m + m + 1; rr <= m * m + 2 * m; rr++) {
            int r = rr - 1;
            b[r] = -bbs.rho[r - m * m - m];
            for (int kk = 1; kk <= m * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == rr - m * m - m) {
                    A[r][k] = -bbs.rho[getII(kk) - 1];
                }
            }
        }
        // Stochastic matrix: <= 1
        for (int rr = m * m + 2 * m + 1; rr <= m * m + 3 * m; rr++) {
            int r = rr - 1;
            b[r] = 1.;
            for (int kk = 1; kk <= m * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == rr - m * m - 2 * m) {
                    A[r][k] = 1.;
                }
            }
        }
        // Stochastic matrix: <= -1
        for (int rr = m * m + 3 * m + 1; rr <= m * m + 4 * m; rr++) {
            int r = rr - 1;
            b[r] = -1.;
            for (int kk = 1; kk <= m * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == rr - m * m - 3 * m) {
                    A[r][k] = -1.;
                }
            }
        }
    }

    /**
     * Index in the vector representation of matrix {@code P} for entry {@code Pij}.
     * @param ii Row index (1 to {@code m}).
     * @param jj Column index (1 to {@code m}).
     * @return Index {@code kk} (1 to {@code m^2}) in the vector representation of matrix
     * {@code P} for entry {@code Pij}.
     */
    int getKK(int ii, int jj) {
        int res = jj + (ii - 1) * m;
        return res;
    }

    /**
     * Row index (1 to {@code m}) in matrix {@code P} for component {@code p_kk} of its vector representation.
     * @param kk Index of {@code Pij} in its vector representation {@code p_kk} (1 to {@code m^2}).
     * @return Row index (1 to {@code m}) in matrix {@code P} for component {@code p_kk} of its vector
     * representation.
     */
    int getII(int kk) {
        int res = 1 + (kk - 1) / m;
        return res;
    }

    /**
     * Column index (1 to {@code m}) in matrix {@code P} for component {@code p_kk} of its vector representation.
     * @param kk Index of {@code Pij} in its vector representation {@code p_kk} (1 to {@code m^2}).
     * @return Column index (1 to {@code m}) in matrix {@code P} for component {@code p_kk} of its vector
     * representation.
     */
    int getJJ(int kk) {
        int res = kk - (getII(kk) - 1) * m;
        return res;
    }

    /**
     * For testing.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        HCModelLIC mlic = new HCModelLIC(bbs);
        Random rand = new Random();
        double[][] P = new double[bbs.m][bbs.m];
        for (int i = 0; i < bbs.m; i++) {
            for (int j = 0; j < bbs.m; j++) {
                P[i][j] = rand.nextDouble();
            }
        }
        double[] p = HCModel.toVector(P);
        ProjOnPolyHSet proj = new ProjOnPolyHSet(p, mlic.A, mlic.b, rand);
        P = HCModel.toMatrix(bbs.m, proj.optimum);
        MatrixOp.printMat(P);
        System.out.println("stochastic matrix");
        for (int i = 0; i < bbs.m; i++) {
            double t = 0.;
            for (int j = 0; j < bbs.m; j++) {
                t += P[i][j];
            }
            System.out.println(t);
        }
        System.out.println("balance equations");
        for (int i = 0; i < bbs.m; i++) {
            double t = 0.;
            for (int j = 0; j < bbs.m; j++) {
                t += bbs.rho[j] * P[j][i];
            }
            System.out.println(t + "\t" + bbs.rho[i]);
        }
    }
}