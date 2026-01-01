package heparansulfate;

import java.util.Random;

/**
 * Linear inequality constraints on a Nonhomogeneity and Independence (N&amp;I) model.
 * Includes constraints for nonnegativity, sum to 1 at each position, and overall 
 * disaccharide abundances.
 */
public class NIModelLIC {
    /**
     * Matrix in the constraint system {@code Ax <= b}.
     */
    double[][] A = null;
    /**
     * Vector in the constraint system {@code Ax <= b}.
     */
    double[] b = null;
    /**
     * Chain length.
     */
    int n = 0;
    /**
     * Number of building block types.
     */
    int m = 0;

    /**
     * Linear inequality constraints on a N&amp;I model: {@code Ax <= b} including 
     * nonnegativity, sum to 1 at each position, and overall disaccharide abundances.
     * @param n Chain length in the model.
     * @param bbs Set of disaccharides.
     * @param zeta Composition values at the fixed positions.
     * @param pos Array of specific positions with composition fixed by {@code zeta}.
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
     * Linear inequality constraints on a N&amp;I model: {@code Ax <= b} including 
     * nonnegativity, sum to 1 at each position, and overall disaccharide abundances.
     * @param n Chain length in the model.
     * @param bbs Set of building blocks (disaccharides).
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
     * Index {@code kk} (1 to {@code n * m}) in vector {@code gamma} corresponding to 
     * entry {@code (ii, jj)} in matrix {@code Gamma}.
     * @param ii Row index in matrix {@code Gamma} (1 to {@code n}).
     * @param jj Column index in matrix {@code Gamma} (1 to {@code m}).
     * @return Index {@code kk} in the vector representation.
     */
    int getKK(int ii, int jj) {
        int res = ii + (jj - 1) * n;
        return res;
    }

    /**
     * Column index {@code jj} (1 to {@code m}) in matrix {@code Gamma} corresponding 
     * to index {@code kk} (1 to {@code n * m}) in vector {@code gamma}.
     * @param kk Index in vector {@code gamma}.
     * @return Column index {@code jj} in matrix {@code Gamma}.
     */
    int getJJ(int kk) {
        int res = 1 + (kk - 1) / n;
        return res;
    }

    /**
     * Row index {@code ii} (1 to {@code n}) in matrix {@code Gamma} corresponding 
     * to index {@code kk} (1 to {@code n * m}) in vector {@code gamma}.
     * @param kk Index in vector {@code gamma}.
     * @return Row index {@code ii} in matrix {@code Gamma}.
     */
    int getII(int kk) {
        int res = kk - (getJJ(kk) - 1) * n;
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