package heparansulfate;

import java.util.Random;

/**
 * linear inequality constraints on a N&I model: nonnegativity, sum to 1 at each position
 * and overall disaccharide abundances
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
     * linear inequality constraints on a N&I model:
     * Ax <= b including nonnegativity, sum to 1 at each position and
     * overall disaccharide abundances
     * @param n chain length in PModel
     * @param bbs set of disaccharides
     * @param pos positions with composition fixed by zeta
     * @param zeta composition at positions in pos
     */
    public NIModelLIC(int n, BBSet bbs, double[][] zeta, int[] pos) {
        this.n = n;
        m = bbs.m;
        A = new double[m * n + 2 * m + 2 * n + 2 * (pos.length * (m - 1))][n * m];
        b = new double[m * n + 2 * m + 2 * n + 2 * (pos.length * (m - 1))];

        for (int rr = 1; rr <= n * m; rr++) {
            int r = rr - 1;
            b[r] = 0.;
            A[r][r] = -1.;
        }
        for (int rr = n * m + 1; rr <= n * m + m; rr++) {
            int r = rr - 1;
            b[r] = bbs.rho[r - n * m] * (double)(n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == rr - n * m && getII(kk) != 1) {
                    A[r][k] = 1.;
                }
            }
        }
        for (int rr = n * m + m + 1; rr <= n * m + 2 * m; rr++) {
            int r = rr - 1;
            b[r] = -bbs.rho[r - n * m - m] * (double)(n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == rr - n * m - m && getII(kk) != 1) {
                    A[r][k] = -1.;
                }
            }
        }
        for (int rr = n * m + 2 * m + 1; rr <= n * m + 2 * m + n; rr++) {
            int r = rr - 1;
            b[r] = 1.;
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == rr - n * m - 2 * m) {
                    A[r][k] = 1.;
                }
            }
        }
        for (int rr = n * m + 2 * m + n + 1; rr <= n * m + 2 * m + 2 * n; rr++) {
            int r = rr - 1;
            b[r] = -1.;
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == rr - n * m - 2 * m - n) {
                    A[r][k] = -1.;
                }
            }
        }
        int rr = n * m + 2 * m + 2 * n;
        for (int p = 0; p < pos.length; p++) {
            int pp = pos[p] + 1;
            for (int z = 0; z < m - 1; z++) {
                int zz = z + 1;
                rr++;
                int r = rr - 1;
                b[r] = zeta[p][z];
                for (int kk = 1; kk <= n * m; kk++) {
                    int k = kk - 1;
                    if (getII(kk) == pp && getJJ(kk) == zz) {
                        A[r][k] = 1.;
                    }
                }
                rr++;
                r = rr - 1;
                b[r] = -zeta[p][z];
                for (int kk = 1; kk <= n * m; kk++) {
                    int k = kk - 1;
                    if (getII(kk) == pp && getJJ(kk) == zz) {
                        A[r][k] = -1.;
                    }
                }
            }
        }
    }

    /**
     * linear inequality constraints on a N&I model:
     * Ax <= b including nonnegativity, sum to 1 at each position and overall
     * disaccharide abundances
     * @param n chain length in PModel
     * @param bbs set of building blocks (disaccharides)
     */
    public NIModelLIC(int n, BBSet bbs) {
        this.n = n;
        m = bbs.m;
        A = new double[m * n + 2 * m + 2 * n][n * m];
        b = new double[m * n + 2 * m + 2 * n];
        for (int rr = 1; rr <= n * m; rr++) {
            int r = rr - 1;
            b[r] = 0.;
            A[r][r] = -1.;
        }
        for (int rr = n * m + 1; rr <= n * m + m; rr++) {
            int r = rr - 1;
            b[r] = bbs.rho[r - n * m] * (double)(n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == rr - n * m && getII(kk) != 1) {
                    A[r][k] = 1.;
                }
            }
        }
        for (int rr = n * m + m + 1; rr <= n * m + 2 * m; rr++) {
            int r = rr - 1;
            b[r] = -bbs.rho[r - n * m - m] * (double)(n - 1);
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getJJ(kk) == rr - n * m - m && getII(kk) != 1) {
                    A[r][k] = -1.;
                }
            }
        }
        for (int rr = n * m + 2 * m + 1; rr <= n * m + 2 * m + n; rr++) {
            int r = rr - 1;
            b[r] = 1.;
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == rr - n * m - 2 * m) {
                    A[r][k] = 1.;
                }
            }
        }
        for (int rr = n * m + 2 * m + n + 1; rr <= n * m + 2 * m + 2 * n; rr++) {
            int r = rr - 1;
            b[r] = -1.;
            for (int kk = 1; kk <= n * m; kk++) {
                int k = kk - 1;
                if (getII(kk) == rr - n * m - 2 * m - n) {
                    A[r][k] = -1.;
                }
            }
        }
    }

    /**
     * index kk (1 to n*m) in vector gamma corresponding to entry (ii,jj) in matrix Gamma
     * @param ii row index in matrix Gamma (1 to n)
     * @param jj column index in matrix Gamma (1 to m)
     * @return index kk (1 to n*m) in vector gamma corresponding to entry
     * (ii,jj) in matrix Gamma
     */
    int getKK(int ii, int jj) {
        int res = ii + (jj - 1) * n;
        return res;
    }

    /**
     * column index jj (1 to m) in matrix Gamma corresponding to index kk
     * (1 to n*m) in vector gamma
     * @param kk index (1 to n*m) in vector gamma
     * @return column index jj (1 to m) in matrix Gamma corresponding to
     * index kk (1 to n*m) in vector gamma
     */
    int getJJ(int kk) {
        int res = 1 + (kk - 1) / n;
        return res;
    }

    /**
     * row index ii (1 to n) in matrix Gamma corresponding to index kk
     * (1 to n*m) in vector gamma
     * @param kk index (1 to n*m) in vector gamma
     * @return row index ii (1 to n) in matrix Gamma corresponding to
     * index kk (1 to n*m) in vector gamma
     */
    int getII(int kk) {
        int res = kk - (getJJ(kk) - 1) * n;
        return res;
    }

    /**
     * for testing
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        int n = 16;
        Random rand = new Random();
        int m = bbs.m;
        int[] pos = new int[1];
        pos[0] = n - 1;
        double[][] zeta = new double[pos.length][m];
        zeta[0][0] = 0.95;
        NIModelLIC plic = new NIModelLIC(n, bbs, zeta, pos);
        for (int sim = 0; sim < 1; sim++) {
            double[][] gam = new double[n][m];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    gam[i][j] = rand.nextDouble();
                }
            }
            ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gam),
                    plic.A, plic.b, rand);
            gam = NIModel.toMatrix(proj.optimum, n);
            MatrixOp.printMat(gam);
        }
        MatrixOp.printVec(bbs.rho);
        System.out.println("Note: projection is not perfect, "
                + "yielding sometimes entries < 0 or > 1."
                + " This is accounted for (corrected) in class NIModelSA");
    }
}