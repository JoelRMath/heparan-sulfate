package heparansulfate;

/**
 * This class represents a linear constraint Ax = b in standard form:
 * matrix A (m x n) has rank {@code m < n} and {@code b >= 0}.
 */
public class SFormCons {
    /**
     * matrix for linear equality Ax = b
     */
    double[][] A = null;
    /**
     * vector in constraint Ax = b, matrix A being transformed
     * so that {@code b >= 0}
     */
    double[] b = null;
    /**
     * number of columns
     */
    int n = 0;
    /**
     * number of rows (rank of A)
     */
    int m = 0;

    /**
     * creates a linear equality constraint in standard form (Ax = b with
     * {@code b >= 0}) starting from a m x n matrix A (of rank m) and a vector b
     * @param mA m x n matrix of rank {@code m < n}
     * @param mb m dimensional vector
     */
    public SFormCons(double[][] mA, double[] mb) {
        m = mA.length;
        n = mA[0].length;
        A = new double[m][n];
        b = new double[m];
        for (int i = 0; i < m; i++) {
            if (mb[i] >= 0.) {
                b[i] = mb[i];
                for (int j = 0; j < n; j++) {
                    A[i][j] = mA[i][j];
                }
            } else {
                b[i] = -mb[i];
                for (int j = 0; j < n; j++) {
                    A[i][j] = -mA[i][j];
                }
            }
        }
    }

    /**
     * prints out the constraint
     */
    void print() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(A[i][j] + "\t");
            }
            System.out.println(b[i]);
        }
    }

    /**
     * for testing
     * @param args command line arguments
     */
    public static void main(String[] args) {
        double[][] mA = new double[2][3];
        double[] mb = new double[2];
        mA[0][0] = 1.;
        mA[0][1] = 2.;
        mA[0][2] = -3.;
        mb[0] = -4.;
        mA[1][0] = 5.;
        mA[1][1] = 6.;
        mA[1][2] = -7.;
        mb[1] = -8.;
        SFormCons sfc = new SFormCons(mA, mb);
        sfc.print();
    }
}