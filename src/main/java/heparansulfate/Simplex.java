package heparansulfate;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.TreeSet;

/**
 * simplex method for linear programming (phase II)
 */
public class Simplex {
    /**
     * simplex tableau
     */
    Tableau tab = null;
    /**
     * values of certificate are: -1 for no feasible point, 0 for unbounded
     * below in phase II and 1 for optimum found
     */
    int certificate = -3;
    /**
     * certificate in words
     */
    String certificateString = null;
    /**
     * buffer copy of the tableau
     */
    double[][] y = null;
    /**
     * basic variables are always kept in the first m columns,
     * this hashtable maps variable indices to column indices
     */
    Hashtable<Integer, Integer> var2col = new Hashtable<Integer, Integer>();
    /**
     * basic variables are always kept in the first m columns,
     * this hashtable maps column indices to variable indices
     */
    Hashtable<Integer, Integer> col2var = new Hashtable<Integer, Integer>();
    /**
     * used to keep track of which variables indices are basic,
     * in order to implement Bland’s rule
     */
    TreeSet<Integer> basic = new TreeSet<Integer>();
    /**
     * used to keep track of which variables indices are nonbasic,
     * in order to implement Bland’s rule
     */
    TreeSet<Integer> nonbasic = new TreeSet<Integer>();
    /**
     * feasible point obtained at the end of phaseI
     */
    double[] feasiblePoint = null;
    /**
     * the point at the end of phase II
     */
    public double[] optimum = null;
    /**
     * inconsistency threshold for phase II
     */
    double epsilon = 1e-5;
    /**
     * number of constraints in phase II
     */
    int m = 0;
    /**
     * number of variables in phase II
     */
    int n = 0;
    /**
     * solution, i.e. final cost
     */
    public double finalCost = 0.;

    /**
     * simplex method for linear programming
     * @param sp1 phase I of the simplex (tableau)
     * @param c vector of cost coefficients
     */
    public Simplex(SimplexPhaseI sp1, double[] c) {
        init(sp1, c);
        boolean phaseII = true;
        while (phaseII) {
            int q = getQ();
            int p = -1;
            if (q == -1) {
                phaseII = false;
                certificate = 1; // Optimum found
            } else {
                p = getP(q);
                if (p == -1) {
                    certificate = 0; // Unbounded
                    phaseII = false;
                } else {
                    pivot(p, q);
                    saveTab();
                }
            }
        }
        optimum = new double[n];
        for (Integer var : basic) {
            Integer col = var2col.get(var);
            optimum[var] = tab.y[col][tab.n];
        }
        finalCost = 0.;
        for (int i = 0; i < c.length; i++) {
            finalCost += c[i] * optimum[i];
        }
    }

    public double getCost() {
        return -tab.y[tab.nrows - 1][tab.ncols - 1];
    }

    void pivot(int p, int q) {
        for (int i = 0; i < tab.nrows; i++) {
            if (i == p) {
                for (int j = 0; j < tab.ncols; j++) {
                    tab.y[p][j] = y[p][j] / y[p][q];
                }
            } else {
                for (int j = 0; j < tab.ncols; j++) {
                    tab.y[i][j] = y[i][j] - y[p][j] * y[i][q] / y[p][q];
                }
            }
        }
        swapColumns(p, q);
    }

    void swapColumns(int p, int q) {
        double[] buf = new double[tab.nrows];
        for (int i = 0; i < tab.nrows; i++) {
            buf[i] = tab.y[i][p];
        }
        for (int i = 0; i < tab.nrows; i++) {
            tab.y[i][p] = tab.y[i][q];
        }
        for (int i = 0; i < tab.nrows; i++) {
            tab.y[i][q] = buf[i];
        }

        Integer varp = col2var.get(p);
        Integer varq = col2var.get(q);

        col2var.put(p, varq);
        var2col.put(varq, p);
        col2var.put(q, varp);
        var2col.put(varp, q);

        basic.remove(varp);
        nonbasic.add(varp);
        nonbasic.remove(varq);
        basic.add(varq);
    }

    int getP(int q) {
        int res = -1;
        double min = 0.;
        boolean first = true;
        for (Integer var : basic) {
            int row = var2col.get(var);
            if (tab.y[row][q] > 0.) {
                double val = tab.y[row][tab.n] / tab.y[row][q];
                if (first) {
                    min = val;
                    first = false;
                    res = row;
                } else {
                    if (val < min) {
                        min = val;
                        res = row;
                    }
                }
            }
        }
        return res;
    }

    int getQ() {
        int res = -1;
        double minr = 0.;
        for (Integer var : nonbasic) {
            int col = var2col.get(var);
            if (tab.y[tab.m][col] < minr) {
                minr = tab.y[tab.m][col];
                res = col;
            }
        }
        return res;
    }

    void init(SimplexPhaseI sp1, double[] c) {
        m = sp1.m;
        n = sp1.n;
        tab = new Tableau(sp1, c);
        y = new double[tab.nrows][tab.ncols];
        saveTab();
        col2var = new Hashtable<Integer, Integer>();
        var2col = new Hashtable<Integer, Integer>();
        Enumeration<Integer> en = sp1.col2var.keys();
        while (en.hasMoreElements()) {
            Integer col = en.nextElement();
            Integer var = sp1.col2var.get(col);
            col2var.put(col, var);
            var2col.put(var, col);
        }
        basic = new TreeSet<Integer>();
        nonbasic = new TreeSet<Integer>();
        basic.addAll(sp1.basic);
        nonbasic.addAll(sp1.nonbasic);

        feasiblePoint = new double[n];
        System.arraycopy(sp1.feasiblePoint, 0, feasiblePoint, 0, n);
    }

    void saveTab() {
        for (int i = 0; i < tab.nrows; i++) {
            System.arraycopy(tab.y[i], 0, y[i], 0, tab.ncols);
        }
    }
}