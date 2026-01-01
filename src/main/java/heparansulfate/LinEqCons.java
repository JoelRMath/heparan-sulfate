package heparansulfate;

import java.util.List;

/**
 * Class that represents linear equality or inequality constraints: {@code Ax = b} or {@code Ax <= b},
 * where {@code A} is an {@code m} by {@code n} matrix and {@code b} is an {@code m}-dimensional vector. 
 * The main utility of this class is the merging of multiple linear constraints.
 */
public class LinEqCons {
    /**
     * Number of rows (constraints).
     */
    int m = 0;
    /**
     * Number of columns (variables).
     */
    int n = 0;
    /**
     * Matrix in the linear system {@code Ax = b}.
     */
    public double[][] A = null;
    /**
     * Constant vector in the linear system {@code Ax = b}.
     */
    public double[] b = null;
    /**
     * Labels for each row, identifying the type of constraint.
     */
    public String[] type = null;

    /**
     * Constructs a linear equality constraint: {@code Ax = b}.
     * @param A Matrix in {@code Ax = b}.
     * @param b Vector in {@code Ax = b}.
     * @param label Row labels for identification.
     */
    public LinEqCons(double[][] A, double[] b, String[] label) {
        m = A.length;
        n = A[0].length;
        this.A = new double[m][n];
        this.b = new double[m];
        type = new String[m];
        for (int i = 0; i < m; i++) {
            this.b[i] = b[i];
            type[i] = label[i];
            for (int j = 0; j < n; j++) {
                this.A[i][j] = A[i][j];
            }
        }
    }

    /**
     * Merges an array of linear constraints into a single {@code LinEqCons} object.
     * @param lec Array of linear constraints to merge.
     */
    public LinEqCons(LinEqCons[] lec) {
        m = 0;
        n = lec[0].n;
        for (int k = 0; k < lec.length; k++) {
            m += lec[k].m;
        }
        A = new double[m][n];
        b = new double[m];
        type = new String[m];
        int i = -1;
        for (int k = 0; k < lec.length; k++) {
            for (int l = 0; l < lec[k].m; l++) {
                i++;
                b[i] = lec[k].b[l];
                type[i] = lec[k].type[l];
                for (int j = 0; j < lec[k].n; j++) {
                    A[i][j] = lec[k].A[l][j];
                }
            }
        }
    }

    /**
     * Merges a {@code List} of linear constraints into a single {@code LinEqCons} object.
     * @param lec List containing linear constraints to merge.
     */
    public LinEqCons(List<LinEqCons> lec) {
        m = 0;
        LinEqCons lc = lec.get(0);
        n = lc.n;
        for (int k = 0; k < lec.size(); k++) {
            lc = lec.get(k);
            m += lc.m;
        }
        A = new double[m][n];
        b = new double[m];
        type = new String[m];
        int i = -1;
        for (int k = 0; k < lec.size(); k++) {
            lc = lec.get(k);
            for (int l = 0; l < lc.m; l++) {
                i++;
                b[i] = lc.b[l];
                type[i] = lc.type[l];
                for (int j = 0; j < lc.n; j++) {
                    A[i][j] = lc.A[l][j];
                }
            }
        }
    }

    /**
     * Returns a tab-delimited String representation of this constraint set:
     * {@code type \t b \t A[][0] \t A[][1] \t ... \t A[][n-1]}.
     * @return A String representation of this {@code LinEqCons}.
     */
    public String getStringRep() {
        String res = "";
        for (int i = 0; i < m; i++) {
            res += type[i];
            res += "\t" + b[i];
            for (int j = 0; j < n; j++) {
                res += "\t" + A[i][j];
            }
            res += "\n";
        }
        return res;
    }

    /**
     * Transforms equality constraints {@code Ax = b} into a set of equivalent inequalities: 
     * {@code Ax <= b} and {@code -Ax <= -b}.
     * @param lec Set of equality constraints.
     * @return Inequality version of the provided constraints.
     */
    public static LinEqCons equalitiesToInequalities(LinEqCons lec) {
        int n = lec.n;
        int m = 2 * lec.m;
        double[][] A = new double[m][n];
        double[] b = new double[m];
        String[] tp = new String[m];
        int row = -1;
        for (int i = 0; i < lec.m; i++) {
            row++;
            b[row] = lec.b[i];
            tp[row] = lec.type[i];
            for (int j = 0; j < n; j++) {
                A[row][j] = lec.A[i][j];
            }
        }
        for (int i = 0; i < lec.m; i++) {
            row++;
            b[row] = -lec.b[i];
            tp[row] = lec.type[i];
            for (int j = 0; j < n; j++) {
                A[row][j] = -lec.A[i][j];
            }
        }
        return new LinEqCons(A, b, tp);
    }

    /**
     * Removes the last row (constraint) from a {@code LinEqCons} object.
     * @param lec Original linear constraint set.
     * @return A new {@code LinEqCons} object excluding the last row of the input.
     */
    public static LinEqCons removeLastRow(LinEqCons lec) {
        int m = lec.m - 1;
        int n = lec.n;
        double[][] A = new double[m][n];
        double[] b = new double[m];
        String[] lab = new String[m];
        for (int i = 0; i < m; i++) {
            b[i] = lec.b[i];
            lab[i] = lec.type[i];
            for (int j = 0; j < n; j++) {
                A[i][j] = lec.A[i][j];
            }
        }
        return new LinEqCons(A, b, lab);
    }
}