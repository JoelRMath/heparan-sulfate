package heparansulfate;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * Implementation of the Simplex method, but limited to Phase I.
 * Phase I is implemented via the addition of {@code m} artificial variables to find 
 * an initial feasible solution.
 */
public class SimplexPhaseI {
    /**
     * Contains the linear constraint {@code Ax = b} in standard form.
     * Elements of {@code b} are ensured to be positive.
     */
    SFormCons sfc = null;

    /**
     * Simplex tableau.
     */
    Tableau tab = null;

    /**
     * Values of certificate are: {@code -1} for no feasible point, {@code 0} for 
     * unbounded below in Phase II, and {@code 1} for optimum found.
     */
    int certificate = -3;

    /**
     * Certificate description in words.
     */
    public String certificateString = null;

    /**
     * The largest absolute value among artificial variables at the end of Phase I.
     * The problem is feasible only if {@code absAVmax < epsilon}.
     */
    double absAVmax = 0.;

    /**
     * Buffer copy of the tableau.
     */
    double[][] y = null;

    /**
     * Maps variable indices to column indices. Basic variables are kept in the first {@code m} columns.
     */
    Map<Integer, Integer> var2col = new HashMap<>();

    /**
     * Maps column indices to variable indices. Basic variables are kept in the first {@code m} columns.
     */
    Map<Integer, Integer> col2var = new HashMap<>();

    /**
     * Keeps track of basic variable indices to implement Bland's rule.
     */
    TreeSet<Integer> basic = new TreeSet<>();

    /**
     * Keeps track of nonbasic variable indices to implement Bland's rule.
     */
    TreeSet<Integer> nonbasic = new TreeSet<>();

    /**
     * Keeps track of indices of the artificial variables.
     */
    Set<Integer> artificial = new HashSet<>();

    /**
     * The point obtained at the end of Phase I, including artificial variables.
     */
    double[] phaseIPoint = null;

    /**
     * The feasible point obtained at the end of Phase I, omitting artificial variables.
     */
    public double[] feasiblePoint = null;

    /**
     * Inconsistency threshold for Phase I.
     */
    double epsilon = 1e-14;

    /**
     * Number of constraints.
     */
    int m = 0;

    /**
     * Number of variables for Phase II.
     */
    public int n = 0;

    /**
     * Initial tableau created for Phase II at the end of Phase I.
     */
    double[][] phaseIITableau = null;

    /**
     * Final infeasibility value.
     */
    double finalCost = 0.;

    /**
     * Implementation of the Simplex method, limited to Phase I.
     * Find a feasible starting point via the addition of {@code m} artificial variables.
     * 
     * @param A Matrix in constraint {@code Ax = b}.
     * @param b Vector in constraint {@code Ax = b}.
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
     * Creates an initial tableau for Phase II and updates mappings between 
     * variables and columns.
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
     * Saves the found feasible point. 
     * {@code fixDegeneracy()} and {@code makePhaseIPoint()} must be called first.
     */
    void makeFeasiblePoint() {
        certificateString = "feasible";
        feasiblePoint = new double[n];
        for (int j = 0; j < n; j++) {
            feasiblePoint[j] = phaseIPoint[m + j];
        }
    }

    /**
     * Saves the final point of Phase I, including values for artificial variables.
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
     * Handles degeneracy by swapping any artificial variables remaining in the 
     * basis with non-artificial nonbasic variables.
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
     * Returns the current cost (objective function value).
     * @return The current cost.
     */
    double getCost() {
        return -tab.y[tab.y.length - 1][tab.y[0].length - 1];
    }

    /**
     * Performs a pivoting operation in the tableau.
     * @param p Index of the row (column to leave basis).
     * @param q Index of the entering column.
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
     * Swaps columns {@code p} and {@code q} in the tableau and updates mappings.
     * @param p Index of the column to leave the basis.
     * @param q Index of the column to enter the basis.
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
     * Finds the index of the row/column to leave the basis using Bland's rule.
     * @param q Index of the column entering the basic set.
     * @return Index of the column to leave the basis, or {@code -1} if unbounded.
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
     * Returns the index of the column to enter the basic set.
     * @return The index {@code q} of the entering column, or {@code -1} if optimal.
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
     * Buffers the current tableau.
     */
    void saveTab() {
        for (int i = 0; i < tab.nrows; i++) {
            for (int j = 0; j < tab.ncols; j++) {
                y[i][j] = tab.y[i][j];
            }
        }
    }

    /**
     * For testing.
     * @param args Command line arguments.
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