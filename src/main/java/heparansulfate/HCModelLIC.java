package heparansulfate;

import java.util.Random;

/**
 * Linear inequality constraints on a H&C model (homogeneous Markov model):
 * nonnegativity, balance equation and stochastic matrix.
 */
public class HCModelLIC {
    /**
     * matrix in Ax <= b
     */
    double[][] A = null;
    /**
     * vector in Ax <= b
     */
    double[] b = null;
    /**
     * number of building blocks
     */
    int m = 0;

    /**
     * Linear inequality constraints on a H&C model.
     * @param bbs set of building blocks
     */
    public HCModelLIC(BBSet bbs) {
        m = bbs.m;
        // Total rows: m*m (non-neg) + m (bal <=) + m (bal >=) + m (stoch <=) + m (stoch >=)
        A = new double[m * (m + 4)][m * m];
        b = new double[m * (m + 4)];

        // 1. Non-negativity: -Pij <= 0
        for (int rr = 1; rr <= m * m; rr++) {
            int r = rr - 1;
            A[r][r] = -1.;
            b[r] = 0.;
        }

        // 2. Balance equations: sum_i(rho_i * P_ij) <= rho_j
        for (int rr = m * m + 1; rr <= m * m + m; rr++) {
            int r = rr - 1;
            int targetJJ = rr - m * m; // 1 to m
            b[r] = bbs.rho[targetJJ - 1];
            for (int kk = 1; kk <= m * m; kk++) {
                if (getJJ(kk) == targetJJ) {
                    A[r][kk - 1] = bbs.rho[getII(kk) - 1];
                }
            }
        }

        // 3. Balance equations: -sum_i(rho_i * P_ij) <= -rho_j
        for (int rr = m * m + m + 1; rr <= m * m + 2 * m; rr++) {
            int r = rr - 1;
            int targetJJ = rr - m * m - m;
            b[r] = -bbs.rho[targetJJ - 1];
            for (int kk = 1; kk <= m * m; kk++) {
                if (getJJ(kk) == targetJJ) {
                    A[r][kk - 1] = -bbs.rho[getII(kk) - 1];
                }
            }
        }

        // 4. Stochastic matrix: sum_j(P_ij) <= 1
        for (int rr = m * m + 2 * m + 1; rr <= m * m + 3 * m; rr++) {
            int r = rr - 1;
            int targetII = rr - m * m - 2 * m;
            b[r] = 1.;
            for (int kk = 1; kk <= m * m; kk++) {
                if (getII(kk) == targetII) {
                    A[r][kk - 1] = 1.;
                }
            }
        }

        // 5. Stochastic matrix: -sum_j(P_ij) <= -1
        for (int rr = m * m + 3 * m + 1; rr <= m * m + 4 * m; rr++) {
            int r = rr - 1;
            int targetII = rr - m * m - 3 * m;
            b[r] = -1.;
            for (int kk = 1; kk <= m * m; kk++) {
                if (getII(kk) == targetII) {
                    A[r][kk - 1] = -1.;
                }
            }
        }
    }

    /**
     * index in vector representation (1 to m^2) for entry Pij
     */
    int getKK(int ii, int jj) {
        return jj + (ii - 1) * m;
    }

    /**
     * row index (1 to m) for vector index kk
     */
    int getII(int kk) {
        return 1 + (kk - 1) / m;
    }

    /**
     * column index (1 to m) for vector index kk
     */
    int getJJ(int kk) {
        return kk - (getII(kk) - 1) * m;
    }

    public static void main(String[] args) {
        String inDir = "input/";
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
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

        System.out.println("Resulting Matrix P:");
        MatrixOp.printMat(P);

        System.out.println("\nStochasticity Check (Rows should sum to 1.0):");
        for (int i = 0; i < bbs.m; i++) {
            double t = 0.;
            for (int j = 0; j < bbs.m; j++) t += P[i][j];
            System.out.println("Row " + i + ": " + t);
        }

        System.out.println("\nBalance Equation Check (Calculated vs Target):");
        for (int i = 0; i < bbs.m; i++) {
            double t = 0.;
            for (int j = 0; j < bbs.m; j++) t += bbs.rho[j] * P[j][i];
            System.out.println("Col " + i + ": " + t + " (Target: " + bbs.rho[i] + ")");
        }
    }
}