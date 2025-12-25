package heparansulfate;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
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
    Map<Integer, Integer> var2col = new HashMap<>();
    /**
     * basic variables are always kept in the first m columns,
     * this hashtable maps column indices to variable indices
     */
    Map<Integer, Integer> col2var = new HashMap<>();
    /**
     * used to keep track of which variables indices are basic,
     * in order to implement Bland’s rule
     */
    TreeSet<Integer> basic = new TreeSet<>();
    /**
     * used to keep track of which variables indices are nonbasic,
     * in order to implement Bland’s rule
     */
    TreeSet<Integer> nonbasic = new TreeSet<>();
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
    double finalCost = 0.;

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
            } else {
                p = getP(q);
                if (p == -1) {
                    certificate = -3;
                    phaseII = false;
                } else {
                    pivot(p, q);
                    saveTab();
                }
            }
        }
        optimum = new double[n];
        Iterator<Integer> it = basic.iterator();
        while (it.hasNext()) {
            Integer var = it.next();
            Integer col = var2col.get(var);
            optimum[var] = tab.y[col][tab.n];
        }
        finalCost = 0.;
        for (int i = 0; i < c.length; i++) {
            finalCost += c[i] * optimum[i];
        }
    }

    /**
     * cost = objective function value
     * @return objective function value
     */
    public double getCost() {
        return -tab.y[tab.y.length - 1][tab.y[0].length - 1];
    }

    /**
     * pivoting operation in the tableau
     * @param p index of the column to leave the basis
     * @param q index of the column to enter the basis
     */
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

    /**
     * swaps columns p and q in the tableau and updates mappings
     * between column and variable indices
     * @param p index of the column to leave the basis
     * @param q index of the column to enter the basis
     */
    void swapColumns(int p, int q) {
        double[] buf = new double[tab.y.length];
        for (int i = 0; i < buf.length; i++) {
            buf[i] = tab.y[i][p];
        }
        for (int i = 0; i < buf.length; i++) {
            tab.y[i][p] = tab.y[i][q];
        }
        for (int i = 0; i < buf.length; i++) {
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

    /**
     * Given the index q of the column to enter next the basic set, finds the
     * index q of the column to leave the basic set, or -1 if the linear
     * problem is unbounded below. This method also implements Bland’s rule
     * @param q index of the column to enter the basic set
     * @return index of the column to leave the basis
     */
    int getP(int q) {
        int res = -1;
        double min = 0.;
        boolean first = true;
        Iterator<Integer> it = basic.iterator();
        while (it.hasNext()) {
            Integer var = it.next();
            Integer I = var2col.get(var);
            int row = I.intValue();
            if (tab.y[row][q] > 0.) {
                if (first) {
                    min = tab.y[row][tab.n] / tab.y[row][q];
                    first = false;
                    res = row;
                } else {
                    if (tab.y[row][tab.n] / tab.y[row][q] < min) {
                        min = tab.y[row][tab.n] / tab.y[row][q];
                        res = row;
                    }
                }
            }
        }
        return res;
    }

    /**
     * returns the index q of the column to enter the basic set in
     * order to lower the cost, or -1 if the current set of basic
     * variables is optimal
     * @return the index q of the column to enter the basis
     */
    int getQ() {
        int res = -1;
        double minr = 0.;
        Iterator<Integer> it = nonbasic.iterator();
        while (it.hasNext()) {
            Integer var = it.next();
            Integer I = var2col.get(var);
            int col = I.intValue();
            if (tab.y[tab.m][col] < minr) {
                minr = tab.y[tab.m][col];
                res = col;
            }
        }
        return res;
    }

    /**
     * initialization
     * @param sp1 output of phase I
     * @param c vector of cost coefficients
     */
    void init(SimplexPhaseI sp1, double[] c) {
        m = sp1.m;
        n = sp1.n;
        tab = new Tableau(sp1, c);
        y = new double[tab.y.length][tab.y[0].length];
        saveTab();
        col2var = new HashMap<>();
        var2col = new HashMap<>();
        for (Integer col : sp1.col2var.keySet()) {
            Integer var = sp1.col2var.get(col);
            col2var.put(col, var);
            var2col.put(var, col);
        }
        basic = new TreeSet<>();
        nonbasic = new TreeSet<>();
        Iterator<Integer> it = sp1.basic.iterator();
        while (it.hasNext()) {
            basic.add(it.next());
        }
        it = sp1.nonbasic.iterator();
        while (it.hasNext()) {
            nonbasic.add(it.next());
        }
        feasiblePoint = new double[n];
        for (int j = 0; j < n; j++) {
            feasiblePoint[j] = sp1.feasiblePoint[j];
        }
    }

    /**
     * buffers the current tableau
     */
    void saveTab() {
        for (int i = 0; i < tab.nrows; i++) {
            for (int j = 0; j < tab.ncols; j++) {
                y[i][j] = tab.y[i][j];
            }
        }
    }
}