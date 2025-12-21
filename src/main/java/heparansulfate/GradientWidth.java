package heparansulfate;

import java.util.Vector;

/**
 * Utilizes linear programming to estimate bounds of S/U proportions at each position
 * in chains; two different constructors: one for all BKHS chains having same length
 * and one for a mixture of chain lengths.
 */
public class GradientWidth {
    Species sp = null;
    double[][] A = null;
    double[] b = null;
    LinEqCons lec = null;
    SimplexPhaseI sp1 = null;
    double infeasibility = 0.;
    String[] lab = null;
    MixSpecies msp = null;

    /**
     * Utilizes linear programming for a mixture of BKHS chain lengths.
     */
    public GradientWidth(String[] lb, int lmi, int lmx, double sig, double mu, String inDir) {
        this.lab = lb;
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        msp = new MixSpecies(lmi, lmx, sig, mu);
        String[] csFile = new String[2];
        String[] fragFile = new String[2];
        csFile[0] = inDir + "US.hepI.txt";
        csFile[1] = inDir + "US.hepIII.txt";
        fragFile[0] = inDir + "hepI.f.txt";
        fragFile[1] = inDir + "hepIII.f.txt";
        lec = msp.getCompleteLEC(bbs, csFile, fragFile);
        A = lec.A;
        b = lec.b;
        sp1 = new SimplexPhaseI(A, b);
        infeasibility = sp1.finalCost;
        System.err.println("Infeasibility: " + infeasibility);
    }

    /**
     * Constructor for the case of all BKHS chains having same length n.
     */
    public GradientWidth(int m, int n, String[] lab, String inDir) {
        this.lab = lab;
        sp = new Species(m, n);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        String[] csFile = new String[2];
        String[] fragFile = new String[2];
        csFile[0] = inDir + "US.hepI.txt";
        csFile[1] = inDir + "US.hepIII.txt";
        fragFile[0] = inDir + "hepI.f.txt";
        fragFile[1] = inDir + "hepIII.f.txt";
        lec = sp.getCompleteLEC(bbs, csFile, fragFile);
        A = lec.A;
        b = lec.b;
        sp1 = new SimplexPhaseI(A, b);
        infeasibility = sp1.finalCost;
    }

    double[] getCostCoeffNRE(int pos, int bb) {
        double[] res = new double[msp.N];
        for (int s = 0; s < msp.N; s++) {
            if (pos < msp.seq[s].length) {
                if (msp.seq[s][pos] == bb) {
                    res[s] = 1.;
                }
            }
        }
        return res;
    }

    double[] getCostCoeffRE(int pos, int bb) {
        double[] res = new double[msp.N];
        for (int s = 0; s < msp.N; s++) {
            if (pos < msp.seq[s].length) {
                if (msp.seq[s][msp.seq[s].length - 1 - pos] == bb) {
                    res[s] = 1.;
                }
            }
        }
        return res;
    }

    double[] getCostCoeff(int pos, int bb) {
        double[] res = new double[sp.N];
        for (int s = 0; s < sp.N; s++) {
            if (sp.seq[s][pos] == bb) {
                res[s] = 1.;
            }
        }
        return res;
    }

    double getLowerNRE(int pos, int bb) {
        double[] c = getCostCoeffNRE(pos, bb);
        Simplex spx = new Simplex(sp1, c);
        double res = spx.finalCost;
        return Math.min(1.0, Math.max(0.0, res));
    }

    double getLowerRE(int pos, int bb) {
        double[] c = getCostCoeffRE(pos, bb);
        Simplex spx = new Simplex(sp1, c);
        double res = spx.finalCost;
        return Math.min(1.0, Math.max(0.0, res));
    }

    double getUpperNRE(int pos, int bb) {
        double[] c = getCostCoeffNRE(pos, bb);
        for (int i = 0; i < c.length; i++) c[i] *= -1.;
        Simplex spx = new Simplex(sp1, c);
        double res = -spx.finalCost;
        return Math.min(1.0, Math.max(0.0, res));
    }

    double getUpperRE(int pos, int bb) {
        double[] c = getCostCoeffRE(pos, bb);
        for (int i = 0; i < c.length; i++) c[i] *= -1.;
        Simplex spx = new Simplex(sp1, c);
        double res = -spx.finalCost;
        return Math.min(1.0, Math.max(0.0, res));
    }

    double getLower(int pos, int bb) {
        double[] c = getCostCoeff(pos, bb);
        Simplex spx = new Simplex(sp1, c);
        double res = spx.finalCost;
        return Math.min(1.0, Math.max(0.0, res));
    }

    double getUpper(int pos, int bb) {
        double[] c = getCostCoeff(pos, bb);
        for (int i = 0; i < c.length; i++) c[i] *= -1.;
        Simplex spx = new Simplex(sp1, c);
        double res = -spx.finalCost;
        return Math.min(1.0, Math.max(0.0, res));
    }

    void saveBounds(String file) {
        Vector<String> v = new Vector<>();
        StringBuilder s = new StringBuilder("pos");
        for (int i = 0; i < sp.m; i++) {
            s.append("\tlower").append(lab[i]).append("\tupper").append(lab[i]);
        }
        v.add(s.toString());
        for (int pos = 0; pos < sp.n; pos++) {
            s = new StringBuilder(String.valueOf(pos + 1));
            for (int j = 0; j < sp.m; j++) {
                s.append("\t").append(getLower(pos, j));
                s.append("\t").append(getUpper(pos, j));
            }
            v.add(s.toString());
            System.err.println(s);
            Utils.saveFile(v, file);
        }
    }

    void saveBoundsMixSpNRE(String file) {
        Vector<String> v = new Vector<>();
        StringBuilder s = new StringBuilder("posNRE");
        for (int i = 0; i < msp.m; i++) {
            s.append("\tlower").append(lab[i]).append("\tupper").append(lab[i]);
        }
        v.add(s.toString());
        for (int pos = 0; pos < msp.lmax; pos++) {
            s = new StringBuilder(String.valueOf(pos + 1));
            for (int j = 0; j < msp.m; j++) {
                s.append("\t").append(getLowerNRE(pos, j) / msp.lengthAtLeast[pos]);
                s.append("\t").append(getUpperNRE(pos, j) / msp.lengthAtLeast[pos]);
            }
            v.add(s.toString());
            System.err.println(s);
        }
        Utils.saveFile(v, file);
    }

    void saveBoundsMixSpRE(String file) {
        Vector<String> v = new Vector<>();
        StringBuilder s = new StringBuilder("posRE");
        for (int i = 0; i < msp.m; i++) {
            s.append("\tlower").append(lab[i]).append("\tupper").append(lab[i]);
        }
        v.add(s.toString());
        for (int pos = 0; pos < msp.lmax; pos++) {
            s = new StringBuilder(String.valueOf(pos + 1));
            for (int j = 0; j < msp.m; j++) {
                s.append("\t").append(getLowerRE(pos, j) / msp.lengthAtLeast[pos]);
                s.append("\t").append(getUpperRE(pos, j) / msp.lengthAtLeast[pos]);
            }
            v.add(s.toString());
            System.err.println(s);
            Utils.saveFile(v, file);
        }
    }

    public static void bounds(String inDir, String outDir) {
        String[] lab = {"U", "S"};
        GradientWidth gw = new GradientWidth(2, 16, lab, inDir);
        gw.saveBounds(outDir + "GradientWidth.res");
    }

    public static void boundsAndN(String inDir, String outDir) {
        String[] lab = {"U", "S"};
        for (int n = 13; n <= 18; n++) {
            GradientWidth gw = new GradientWidth(2, n, lab, inDir);
            gw.saveBounds(outDir + "GradientWidth.N" + n + ".res");
        }
    }

    public static void defaultMixSpec(String inDir, String outDir) {
        String[] lab = {"U", "S"};
        GradientWidth gw = new GradientWidth(lab, 10, 20, 3.5, 16., inDir);
        gw.saveBoundsMixSpNRE(outDir + "BoundsMixSpNRE.res");
        gw.saveBoundsMixSpRE(outDir + "BoundsMixSpRE.res");
    }

    public static void main(String[] args) {
        // Recommended Mac paths
        String inDir = "input/";
        String outDir = "output/";
        bounds(inDir, outDir);
        boundsAndN(inDir, outDir);
        defaultMixSpec(inDir, outDir);
    }
}