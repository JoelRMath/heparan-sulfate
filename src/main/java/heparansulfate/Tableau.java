package heparansulfate;

/**
 * Simplex tableau used for representing linear programming problems.
 * <p>
 * The tableau stores constraints and cost coefficients, facilitating pivot operations
 * used in both Phase I and Phase II of the simplex method.
 */
public class Tableau {
    /**
     * Number of columns in the tableau: {@code n + 1}.
     */
    int ncols = 0;
    /**
     * Number of variables.
     */
    int n = 0;
    /**
     * Number of rows in the tableau: {@code m + 1}.
     */
    int nrows = 0;
    /**
     * Number of constraints.
     */
    int m = 0;
    /**
     * Simplex tableau matrix data.
     */
    double[][] y = null;
    /**
     * Cost coefficients for the objective function.
     */
    double[] c = null;

    /**
     * Constructor for Phase I: creates a simplex tableau based on a constraint in
     * canonical form (linear constraint in standard form with added artificial variables).
     * 
     * @param avsfc Linear constraint in standard form with added artificial variables.
     */
    public Tableau(AVSFormCons avsfc) {
        m = avsfc.m;
        nrows = m + 1;
        n = avsfc.n;
        ncols = n + 1;
        this.c = new double[n];
        for (int j = 0; j < m; j++) {
            c[j] = 1.;
        }
        y = new double[nrows][ncols];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                y[i][j] = avsfc.A[i][j];
            }
            y[i][n] = avsfc.b[i];
        }
        for (int j = 0; j < n; j++)
            y[m][j] = c[j];
        reduce();
    }

    /**
     * Constructor for Phase II based on the results of Phase I.
     * @param sp1 Phase I solver results containing the feasible tableau.
     * @param coef Cost coefficients for the Phase II objective function.
     */
    public Tableau(SimplexPhaseI sp1, double[] coef) {
        m = sp1.m;
        nrows = m + 1;
        n = sp1.n;
        ncols = n + 1;
        c = new double[n];
        for (int j = 0; j < n; j++) {
            Integer col = sp1.var2col.get(j);
            c[col] = coef[j];
        }
        y = new double[nrows][ncols];
        for (int i = 0; i < sp1.m; i++) {
            for (int j = 0; j < ncols; j++) {
                y[i][j] = sp1.phaseIITableau[i][j];
            }
        }
        for (int j = 0; j < n; j++)
            y[m][j] = c[j];
        reduce();
    }

    /**
     * Constructor for Phase II based on an array-tableau already in canonical form, 
     * with the exception of the last row (cost coefficients) which is reduced here.
     * @param t Array-tableau in canonical form.
     */
    public Tableau(double[][] t) {
        m = t.length - 1;
        nrows = m + 1;
        n = t[0].length - 1;
        ncols = n + 1;
        c = new double[n];
        for (int j = 0; j < n; j++) {
            c[j] = t[m][j];
        }
        y = new double[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                y[i][j] = t[i][j];
            }
        }
        reduce();
    }

    /**
     * Performs pivot operations to reduce the last row of the tableau, ensuring 
     * that coefficients of the basic variables are equal to 0.
     */
    void reduce() {
        for (int i = 0; i < m; i++) {
            double coeff = y[m][i];
            for (int j = 0; j < ncols; j++) {
                y[m][j] = y[m][j] - y[i][j] * coeff;
            }
        }
    }

    /**
     * Prints the tableau to the console.
     */
    void print() {
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < n; j++) {
                double d = y[i][j];
                if (Math.abs(d) < 1e-12)
                    d = 0.;
                System.out.print(d + "\t");
            }
            double d = y[i][n];
            if (Math.abs(d) < 1e-12)
                d = 0.;
            System.out.println(d);
        }
        System.out.println("-------------");
    }

    /**
     * Main entry point for testing.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        double[][] A = new double[2][3];
        double[] b = new double[2];
        A[0][0] = 2.;
        A[0][1] = 1.;
        A[0][2] = 2.;
        b[0] = 4.;
        A[1][0] = 3.;
        A[1][1] = 3.;
        A[1][2] = 1.;
        b[1] = 3.;
        AVSFormCons avsfc = new AVSFormCons(A, b);
        Tableau tab = new Tableau(avsfc);
        tab.print();
    }
}