package heparansulfate;

import java.util.Random;

/**
 * Linear inequality constraints on a N&I model: nonnegativity, sum to 1 at each position
 * and overall disaccharide abundances.
 */
public class NIModelLIC {
    /**
     * matrix in Ax <= b
     */
    double[][] A = null;
    /**
     * vector in Ax <= b
     */
    double[] b = null;
    /**
     * chain length
     */
    int n = 0;
    /**
     * number of building blocks
     */
    int m = 0;

    /**
     * Linear inequality constraints on a N&I model.
     * @param n chain length
     * @param bbs set of disaccharides
     * @param pos positions with composition fixed by zeta
     * @param zeta composition at positions in pos
     */
    public NIModelLIC(int n, BBSet bbs, double[][] zeta, int[] pos) {
        this.n = n;
        m = bbs.m;
        // Calculation of row dimensions for the constraint matrix
        int rows = m * n + 2 * m + 2 * n + 2 * (pos.length * (m - 1));
        A = new double[rows][n * m];
        b = new double[rows];

        // 1. Non-negativity constraints: -x <= 0
        for (int rr = 1; rr <= n * m; rr++) {
            int r = rr - 1;
            b[r] = 0.;
            A[r][r] = -1.;
        }

        // 2. Overall disaccharide abundances: sum_pos(rho_pos) <= target
        for (int rr = n * m + 1; rr <= n * m + m; rr++) {
            int r = rr - 1;
            b[r] = bbs.rho[r - n * m] * (double) (n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == (r - n * m + 1) && getII(kk) != 1) {
                    A[r][k] = 1.;
                }
            }
        }

        // 3. Overall disaccharide abundances: -sum_pos(rho_pos) <= -target
        for (int rr = n * m + m + 1; rr <= n * m + 2 * m; rr++) {
            int r = rr - 1;
            b[r] = -bbs.rho[r - n * m - m] * (double) (n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == (r - n * m - m + 1) && getII(kk) != 1) {
                    A[r][k] = -1.;
                }
            }
        }

        // 4. Sum to 1 at each position: sum_bb(p_bb) <= 1
        for (int rr = n * m + 2 * m + 1; rr <= n * m + 2 * m + n; rr++) {
            int r = rr - 1;
            b[r] = 1.;
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == (rr - n * m - 2 * m)) {
                    A[r][k] = 1.;
                }
            }
        }

        // 5. Sum to 1 at each position: -sum_bb(p_bb) <= -1
        for (int rr = n * m + 2 * m + n + 1; rr <= n * m + 2 * m + 2 * n; rr++) {
            int r = rr - 1;
            b[r] = -1.;
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == (rr - n * m - 2 * m - n)) {
                    A[r][k] = -1.;
                }
            }
        }

        // 6. Zeta constraints: fixing specific positions
        int rr = n * m + 2 * m + 2 * n;
        for (int p = 0; p < pos.length; p++) {
            int pp = pos[p] + 1;
            for (int z = 0; z < m - 1; z++) {
                int zz = z + 1;
                // <= zeta
                rr++;
                int r1 = rr - 1;
                b[r1] = zeta[p][z];
                for (int kk = 1; kk <= n * m; kk++) {
                    if (getII(kk) == pp && getJJ(kk) == zz) {
                        A[r1][kk - 1] = 1.;
                    }
                }
                // >= zeta
                rr++;
                int r2 = rr - 1;
                b[r2] = -zeta[p][z];
                for (int kk = 1; kk <= n * m; kk++) {
                    if (getII(kk) == pp && getJJ(kk) == zz) {
                        A[r2][kk - 1] = -1.;
                    }
                }
            }
        }
    }

    /**
     * Simple constructor without zeta constraints.
     */
    public NIModelLIC(int n, BBSet bbs) {
        this.n = n;
        m = bbs.m;
        int rows = m * n + 2 * m + 2 * n;
        A = new double[rows][n * m];
        b = new double[rows];

        for (int rr = 1; rr <= n * m; rr++) {
            int r = rr - 1;
            b[r] = 0.;
            A[r][r] = -1.;
        }
        for (int rr = n * m + 1; rr <= n * m + m; rr++) {
            int r = rr - 1;
            b[r] = bbs.rho[r - n * m] * (double) (n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                if (getJJ(kk) == (r - n * m + 1) && getII(kk) != 1) {
                    A[r][kk - 1] = 1.;
                }
            }
        }
        for (int rr = n * m + m + 1; rr <= n * m + 2 * m; rr++) {
            int r = rr - 1;
            b[r] = -bbs.rho[r - n * m - m] * (double) (n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                if (getJJ(kk) == (r - n * m - m + 1) && getII(kk) != 1) {
                    A[r][kk - 1] = -1.;
                }
            }
        }
        for (int rr = n * m + 2 * m + 1; rr <= n * m + 2 * m + n; rr++) {
            int r = rr - 1;
            b[r] = 1.;
            for (int kk = 1; kk <= n * m; kk++) {
                if (getII(kk) == (rr - n * m - 2 * m)) {
                    A[r][kk - 1] = 1.;
                }
            }
        }
        for (int rr = n * m + 2 * m + n + 1; rr <= n * m + 2 * m + 2 * n; rr++) {
            int r = rr - 1;
            b[r] = -1.;
            for (int kk = 1; kk <= n * m; kk++) {
                if (getII(kk) == (rr - n * m - 2 * m - n)) {
                    A[r][kk - 1] = -1.;
                }
            }
        }
    }

    int getKK(int ii, int jj) {
        return ii + (jj - 1) * n;
    }

    int getJJ(int kk) {
        return 1 + (kk - 1) / n;
    }

    int getII(int kk) {
        return kk - (getJJ(kk) - 1) * n;
    }

    public static void main(String[] args) {
        String inDir = "input/";
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        int n = 16;
        Random rand = new Random();
        int m = bbs.m;
        int[] pos = {n - 1};
        double[][] zeta = new double[1][m];
        zeta[0][0] = 0.95;

        NIModelLIC plic = new NIModelLIC(n, bbs, zeta, pos);
        double[][] gam = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) gam[i][j] = rand.nextDouble();
        }

        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gam), plic.A, plic.b, rand);
        gam = NIModel.toMatrix(proj.optimum, n);
        
        System.out.println("Resulting Gamma Matrix:");
        MatrixOp.printMat(gam);
        System.out.println("\nNote: Projection is corrected in NIModelSA if entries exceed bounds.");
    }
}