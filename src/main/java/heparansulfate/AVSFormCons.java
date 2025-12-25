import java.util.Random;

/**
 * creates an affine constraint in standard form (Ax = b, with rank(A) < dim(b)
 * and b >= 0)and with artificial variables
 * (for phase I of the simplex)
 */
public class AVSFormCons {
    /**
     * affine constraint in standard form, before adding artificial variables
     */
    SFormCons sfc = null;
    /**
     * number of variables, after adding artificial ones
     */
    int n = 0;
    /**
     * number of equality constraints
     */
    int m = 0;
    /**
     * matrix for the constraints, after adding artificial variables
     */
    double[][] A = null;
    /**
     * vector of constraints, in standard form
     */
    double[] b = null;

    /**
     * creates an affine constraint in standard canonical form (with added
     * artificial
     * variables) for phase I of the simplex: [I A] x = b, with b >= 0
     * 
     * @param mA constraint matrix
     * @param mb constraint values
     */
    public AVSFormCons(double[][] mA, double[] mb) {
        sfc = new SFormCons(mA, mb);
        m = sfc.m;
        n = sfc.m + sfc.n;
        b = new double[m];
        A = new double[m][n];
        for (int i = 0; i < m; i++) {
            b[i] = sfc.b[i];
            A[i][i] = 1.;
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < sfc.n; j++) {
                A[i][j + m] = sfc.A[i][j];
            }
        }
    }
     
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
     * 
     * @param args
     */
    public static void main(String[] args) {
        int M = 10;
        int N = 15;
        Random rand = new Random();
        double[][] mA = new double[M][N];
        double[] mb = new double[M];
        double[] x = new double[N];
        for (int j = 0; j < N; j++)
            x[j] = rand.nextDouble();
        for (int i = 0; i < M; i++) {
            mb[i] = 0.;
            for (int j = 0; j < N; j++) {
                mA[i][j] = (0.5 - rand.nextDouble());
                mb[i] += (mA[i][j] * x[j]);
            }
        }
        new AVSFormCons(mA, mb);
    }
}