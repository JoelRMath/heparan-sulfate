package heparansulfate;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * implementation of the Simplex method, but limited to phase I implemented
 * via addition of m artificial variables
 */
public class SimplexPhaseI {
    /**
     * contains the linear constraint Ax = b in standard for, that is,
     * if some elements of b were negative, they are turned positive and the
     * corresponding row of A is multiplied by -1
     */
    SFormCons sfc = null;
    /**
     * simplex tableau
     */
    Tableau tab = null;
    /**
     * values of certificate are:
     * -1 for no feasible point, 0 for unbounded below in phase II and 1
     * for optimum found
     */
    int certificate = -3;
    /**
     * certificate in words
     */
    public String certificateString = null;
    /**
     * the largest absolute value among artificial variables at the end
     * of phase I (feasible problem only if absAVmax < epsilon)
     */
    double absAVmax = 0.;
    /**
     * buffer of the tableau
     */
    double[][] y = null;
    /**
     * basic variables are always kept in the first m columns, this hashtable
     * maps variable indices to column indices
     */
    Map<Integer, Integer> var2col = new HashMap<>();
    /**
     * basic variables are always kept in the first m columns, this hashtable
     * maps column indices to variable indices
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
     * keeps track of indices of the artificial variables
     */
    Set<Integer> artificial = new HashSet<>();
    /**
     * the point obtained at the end of phaseI, feasible or not
     */
    double[] phaseIPoint = null;
    /**
     * same as phaseIpoint but omitting the artificial variables
     */
    public double[] feasiblePoint = null;
    /**
     * inconsistency threshold for phase I
     */
    double epsilon = 1e-14;
    /**
     * number of constraints
     */
    int m = 0;
    /**
     * number of variables in phase II (i.e. n is n+m in phase I)
     */
    public int n = 0;
    /**
     * at the end of phase I an initial tableau for phase II is created;
     * note that this tableau does not have a last row for relative costs
     */
    double[][] phaseIITableau = null;
    /**
     * final infeasibility
     */
    double finalCost = 0.;

    /**
     * implementation of the Simplex method, limited to phase I
     * implemented via addition of m artificial variables
     * @param A matrix in constraint Ax = b
     * @param b vector in constraint Ax = b
     */
    public SimplexPhaseI(double[][] A, double[] b) {
        sfc = new SFormCons(A, b);
        m = A.length;
        n = A[0].length; // for phase II only, phase I uses tab.n
        // Phase I
        AVSFormCons avsfc = new AVSFormCons(A, b);
        tab = new Tableau(avsfc);
        y = new double[tab.y.length][tab.y[0].length];
        saveTab();
        for (int i = 0; i < tab.n; i++) {
            col2var.put(i, i);
            var2col.put(i, i);
            if (i < tab.m) {
                basic.add(i);
                artificial.add(i);
            } else
                nonbasic.add(i);
        }
        boolean phaseI = true;
        while (phaseI) {
            int q = getQ();
            int p = -1;
            if (q == -1) {
                phaseI = false;
            } else {
                p = getP(q);
                if (p == -1) {
                    certificate = -3;
                    phaseI = false;
                } else {
                    pivot(p, q);
                    saveTab();
                }
            }
            double cost = getCost();
            finalCost = Math.abs(cost);
            if (cost < 1e-13) {
                phaseI = false;
            }
        }
        fixDegeneracy();
        makePhaseIPoint();
        makeFeasiblePoint();
        updateForPhaseII();
    }

    /**
     * creates an initial tableau for phase II and updates mapping
     * between variables and columns
     */
    void updateForPhaseII() {
        Map<Integer, Integer> col2varB = new HashMap<>();
        TreeSet<Integer> basicB = new TreeSet<>();
        TreeSet<Integer> nonbasicB = new TreeSet<>();
        phaseIITableau = new double[m][n + 1];
        // copy of the tableau (without last row) and removal of artificial variables
        int col = -1;
        for (int j = 0; j < tab.n; j++) {
            Integer var = col2var.get(j);
            if (!artificial.contains(var)) {
                col++;
                int varB = var - tab.m;
                col2varB.put(col, varB);
                for (int i = 0; i < tab.m; i++) {
                    phaseIITableau[i][col] = tab.y[i][j];
                }
                if (j < tab.m) {
                    basicB.add(varB);
                } else {
                    nonbasicB.add(varB);
                }
            }
        }
        for (int i = 0; i < tab.m; i++) {
            phaseIITableau[i][n] = tab.y[i][tab.n];
        }
        // update of index mappings
        col2var = new HashMap<>();
        var2col = new HashMap<>();
        basic = new TreeSet<>();
        nonbasic = new TreeSet<>();
        for (Integer cl : col2varB.keySet()) {
            Integer var = col2varB.get(cl);
            col2var.put(cl, var);
            var2col.put(var, cl);
        }
        Iterator<Integer> it = basicB.iterator();
        while (it.hasNext()) {
            basic.add(it.next());
        }
        it = nonbasicB.iterator();
        while (it.hasNext()) {
            nonbasic.add(it.next());
        }
    }

    /**
     * saves the found feasible point; note that fixDegeneracy() and
     * makePhaseIPoint() must have been called before
     */
    void makeFeasiblePoint() {
        certificateString = "feasible";
        feasiblePoint = new double[n];
        for (int j = 0; j < n; j++) {
            feasiblePoint[j] = phaseIPoint[m + j];
        }
    }

    /**
     * saves the last point of phase I (including artificial variables)
     */
    void makePhaseIPoint() {
        phaseIPoint = new double[tab.n];
        Iterator<Integer> it = basic.iterator();
        while (it.hasNext()) {
            Integer var = it.next();
            Integer col = var2col.get(var);
            phaseIPoint[var] = tab.y[col][tab.n];
        }
    }

    /**
     * in case of degeneracy some of the final basic variables might still be artificial
     * variables; if so, these basic variables must be swapped with nonbasic variables
     * that are not artificial variables; this is performed by this method, updating also
     * mappings between tableau columns and variable indices
     */
    void fixDegeneracy() {
        Set<Integer> toswap = new HashSet<>();
        Iterator<Integer> it = basic.iterator();
        while (it.hasNext()) {
            Integer var = it.next();
            if (artificial.contains(var)) {
                toswap.add(var);
            }
        }
        if (toswap.size() > 0) {
            it = toswap.iterator();
            while (it.hasNext()) {
                Integer I = it.next();
                Integer P = var2col.get(I);
                Integer Q = null;
                Iterator<Integer> it2 = nonbasic.iterator();
                while (it2.hasNext()) {
                    Integer J = it2.next();
                    if (!artificial.contains(J)) {
                        Q = var2col.get(J);
                        break;
                    }
                }
                if (Q != null) {
                    pivot(P, Q);
                }
            }
        }
    }

    /**
     * the current cost (objective function value)
     * @return the current cost
     */
    double getCost() {
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
            Integer row = var2col.get(var);
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
     * variables is optimal.
     * @return the index q of the column to enter the basis
     */
    int getQ() {
        int res = -1;
        double minr = 0.;
        Iterator<Integer> it = nonbasic.iterator();
        while (it.hasNext()) {
            Integer var = it.next();
            Integer col = var2col.get(var);
            if (tab.y[tab.m][col] < minr) {
                minr = tab.y[tab.m][col];
                res = col;
            }
        }
        return res;
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

    /**
     * for testing
     * @param args
     */
    public static void main(String[] args) {
        int n = 100;
        int m = 8;
        Random rand = new Random();
        double[][] A = new double[m][n];
        double[] b = new double[m];
        double[] x = new double[n];
        double sum = 0.;
        for (int j = 0; j < n; j++) {
            x[j] = rand.nextDouble();
            sum += x[j];
        }
        for (int j = 0; j < n; j++) {
            x[j] /= sum;
        }
        for (int i = 0; i < m; i++) {
            b[i] = 0.;
            for (int j = 0; j < n; j++) {
                A[i][j] = rand.nextDouble();
                b[i] += x[j] * A[i][j];
            }
        }
        SimplexPhaseI sp1 = new SimplexPhaseI(A, b);
        System.out.println(sp1.finalCost);
        double[] c = new double[n];
        for (int j = 0; j < n; j++) {
            c[j] = 1.;
        }
        Simplex sp = new Simplex(sp1, c);
        System.out.println(sp.finalCost);
    }
}