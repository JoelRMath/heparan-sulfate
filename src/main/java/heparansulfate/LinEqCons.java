package heparansulfate;

import java.util.List;

/**
 * Class that represents linear equality/inequality constraints: Ax = b or Ax <= b,
 * where A is an m by n matrix and b an m vector. The main utility of this class
 * is merging of linear constraints
 */
public class LinEqCons {
    /**
     * number of rows
     */
    int m = 0;
    /**
     * number of columns
     */
    int n = 0;
    /**
     * matrix in Ax = b
     */
    public double[][] A = null;
    /**
     * vector in Ax = b
     */
    public double[] b = null;
    /**
     * label for each row (type of constraint)
     */
    public String[] type = null;

    /**
     * represents a linear equality constraint: Ax = b
     * @param A matrix in Ax = b
     * @param b vector in Ax = b
     * @param label row labels
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
     * Merges several linear constraints (array) into one
     * @param lec array of linear constraints
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
     * Merges several linear constraints (List<LinEqCons>) into one
     * @param lec List containing linear constraints
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
     * Returns a String representation of this LinEqCons:
     * type \t b \t A[][0] \t A[][1] \t ... \t A[][n-1]
     * @return a String representation of this LinEqCons
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
     * transforms Ax = b into Ax <= b and -Ax <= -b
     * @param lec set of equality constraints
     * @return inequality version of Ax = b (Ax <= b and -Ax <= -b)
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
     * removes one row (constraint) in lec
     * @param lec linear constraint
     * @return lec minus its last row (constraint)
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