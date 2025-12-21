package heparansulfate;

import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Random;
import java.util.TreeSet;

/**
 * implementation of the Simplex method, but limited to phase I implemented
 * via addition of m artificial variables
 */
public class SimplexPhaseI {
    SFormCons sfc = null;
    Tableau tab = null;
    int certificate = -3;
    public String certificateString = null;
    double absAVmax = 0.;
    double[][] y = null;
    
    Hashtable<Integer, Integer> var2col = new Hashtable<Integer, Integer>();
    Hashtable<Integer, Integer> col2var = new Hashtable<Integer, Integer>();
    TreeSet<Integer> basic = new TreeSet<Integer>();
    TreeSet<Integer> nonbasic = new TreeSet<Integer>();
    HashSet<Integer> artificial = new HashSet<Integer>();
    
    double[] phaseIPoint = null;
    public double[] feasiblePoint = null;
    double epsilon = 1e-14;
    int m = 0;
    public int n = 0;
    public double[][] phaseIITableau = null;
    public double finalCost = 0.;

    public SimplexPhaseI(double[][] A, double[] b) {
        sfc = new SFormCons(A, b);
        m = A.length;
        n = A[0].length; 
        
        AVSFormCons avsfc = new AVSFormCons(A, b);
        tab = new Tableau(avsfc);
        y = new double[tab.nrows][tab.ncols];
        saveTab();
        
        for (int i = 0; i < tab.n; i++) {
            col2var.put(Integer.valueOf(i), Integer.valueOf(i));
            var2col.put(Integer.valueOf(i), Integer.valueOf(i));
            if (i < tab.m) {
                basic.add(Integer.valueOf(i));
                artificial.add(Integer.valueOf(i));
            } else {
                nonbasic.add(Integer.valueOf(i));
            }
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
            if (finalCost < 1e-13) {
                phaseI = false;
            }
        }
        fixDegeneracy();
        makePhaseIPoint();
        makeFeasiblePoint();
        updateForPhaseII();
    }

    void updateForPhaseII() {
        Hashtable<Integer, Integer> col2varB = new Hashtable<Integer, Integer>();
        TreeSet<Integer> basicB = new TreeSet<Integer>();
        TreeSet<Integer> nonbasicB = new TreeSet<Integer>();
        phaseIITableau = new double[m][n + 1];

        int col = -1;
        for (int j = 0; j < tab.n; j++) {
            Integer var = col2var.get(Integer.valueOf(j));
            if (!artificial.contains(var)) {
                col++;
                Integer varB = Integer.valueOf(var.intValue() - tab.m);
                col2varB.put(Integer.valueOf(col), varB);
                for (int i = 0; i < tab.m; i++) {
                    phaseIITableau[i][col] = tab.y[i][j];
                }
                if (basic.contains(var)) {
                    basicB.add(varB);
                } else {
                    nonbasicB.add(varB);
                }
            }
        }
        for (int i = 0; i < tab.m; i++) {
            phaseIITableau[i][n] = tab.y[i][tab.n];
        }

        col2var = new Hashtable<Integer, Integer>();
        var2col = new Hashtable<Integer, Integer>();
        basic = new TreeSet<Integer>();
        nonbasic = new TreeSet<Integer>();
        
        Enumeration<Integer> en = col2varB.keys();
        while (en.hasMoreElements()) {
            Integer cl = en.nextElement();
            Integer var = col2varB.get(cl);
            col2var.put(cl, var);
            var2col.put(var, cl);
        }
        basic.addAll(basicB);
        nonbasic.addAll(nonbasicB);
    }

    void makeFeasiblePoint() {
        certificateString = "feasible";
        feasiblePoint = new double[n];
        for (int j = 0; j < n; j++) {
            feasiblePoint[j] = phaseIPoint[m + j];
        }
    }

    void makePhaseIPoint() {
        phaseIPoint = new double[tab.n];
        for (Integer var : basic) {
            Integer col = var2col.get(var);
            phaseIPoint[var.intValue()] = tab.y[col.intValue()][tab.n];
        }
    }

    void fixDegeneracy() {
        HashSet<Integer> toswap = new HashSet<Integer>();
        for (Integer var : basic) {
            if (artificial.contains(var)) toswap.add(var);
        }
        if (toswap.size() > 0) {
            for (Integer I : toswap) {
                Integer P = var2col.get(I);
                Integer Q = null;
                for (Integer J : nonbasic) {
                    if (!artificial.contains(J)) {
                        Q = var2col.get(J);
                        break;
                    }
                }
                if (Q != null) pivot(P.intValue(), Q.intValue());
            }
        }
    }

    double getCost() {
        return -tab.y[tab.nrows - 1][tab.ncols - 1];
    }

    void pivot(int p, int q) {
        for (int i = 0; i < tab.nrows; i++) {
            if (i == p) {
                for (int j = 0; j < tab.ncols; j++) tab.y[p][j] = y[p][j] / y[p][q];
            } else {
                for (int j = 0; j < tab.ncols; j++) tab.y[i][j] = y[i][j] - y[p][j] * y[i][q] / y[p][q];
            }
        }
        swapColumns(p, q);
    }

    void swapColumns(int p, int q) {
        double[] buf = new double[tab.nrows];
        for (int i = 0; i < tab.nrows; i++) buf[i] = tab.y[i][p];
        for (int i = 0; i < tab.nrows; i++) tab.y[i][p] = tab.y[i][q];
        for (int i = 0; i < tab.nrows; i++) tab.y[i][q] = buf[i];

        Integer varp = col2var.get(Integer.valueOf(p));
        Integer varq = col2var.get(Integer.valueOf(q));
        
        col2var.put(Integer.valueOf(p), varq);
        var2col.put(varq, Integer.valueOf(p));
        col2var.put(Integer.valueOf(q), varp);
        var2col.put(varp, Integer.valueOf(q));
        
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
            int row = var2col.get(var).intValue();
            if (tab.y[row][q] > 0.) {
                double val = tab.y[row][tab.n] / tab.y[row][q];
                if (first) {
                    min = val;
                    first = false;
                    res = row;
                } else if (val < min) {
                    min = val;
                    res = row;
                }
            }
        }
        return res;
    }

    int getQ() {
        int res = -1;
        double minr = 0.;
        for (Integer var : nonbasic) {
            int col = var2col.get(var).intValue();
            if (tab.y[tab.m][col] < minr) {
                minr = tab.y[tab.m][col];
                res = col;
            }
        }
        return res;
    }

    void saveTab() {
        for (int i = 0; i < tab.nrows; i++) {
            System.arraycopy(tab.y[i], 0, y[i], 0, tab.ncols);
        }
    }

    public static void main(String[] args) {
        int n = 100, m = 8;
        Random rand = new Random();
        double[][] A = new double[m][n];
        double[] b = new double[m], x = new double[n];
        double sum = 0.;
        for (int j = 0; j < n; j++) {
            x[j] = rand.nextDouble();
            sum += x[j];
        }
        for (int j = 0; j < n; j++) x[j] /= sum;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = rand.nextDouble();
                b[i] += x[j] * A[i][j];
            }
        }
        SimplexPhaseI sp1 = new SimplexPhaseI(A, b);
        System.out.println("Phase I Infeasibility: " + sp1.finalCost);
        new Simplex(sp1, new double[n]);
        System.out.println("Simplex Optimization Complete.");
    }
}