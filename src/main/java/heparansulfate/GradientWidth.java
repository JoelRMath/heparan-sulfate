package heparansulfate;

import java.util.ArrayList;
import java.util.List;

/**
 * utilizes linear programming to estimate bounds of S/U proportions at each position
 * in chains; two different constructors: one for all BKHS chains having same length
 * and one for a mixture of chain lengths
 */
public class GradientWidth {
    /**
     * enumeration of all possible sequences, one chain length only
     */
    Species sp = null;
    /**
     * matrix in constraint Ap = b
     */
    double[][] A = null;
    /**
     * vector in constraint Ap = b
     */
    double[] b = null;
    /**
     * convenience wrapper for the constraints
     */
    LinEqCons lec = null;
    /**
     * phase I of the simplex
     */
    SimplexPhaseI sp1 = null;
    /**
     * output of phase I
     */
    double infeasibility = 0.;
    /**
     * disaccharide labels
     */
    String[] lab = null;
    /**
     * enumeration of all sequences, mixture of chain lengths
     */
    MixSpecies msp = null;

    /**
     * utilizes linear programming to estimate bounds of S/U proportions at each position;
     * constructor for a mixture of BKHS chain length (average mu and spread sigma),
     * Note: the second constructor requires large memory (e.g. -Xmx8000M)
     * @param lb disaccharide labels
     * @param lmi smallest chain length
     * @param lmx largest chain length
     * @param sig spread of chain length distribution
     * @param mu average chain length
     * @param inDir input directory (ends with "\\")
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
        System.err.println(infeasibility);
    }

    /**
     * utilizes linear programming to estimate bounds of S/U proportions at each position;
     * constructor for the case of all BKHS chains having same length n
     * @param m number of disaccharides
     * @param n BKHS chain length
     * @param lab disaccharide labels
     * @param inDir input directory (ends with "\\")
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

    /**
     * linear cost coefficients for presence of bb at position pos,
     * when chains are aligned by the NRE; mixture of chain lengths
     * @param pos position in chains (from NRE)
     * @param bb disaccharide
     * @return linear cost coefficients for presence of bb at position pos,
     * when chains are aligned by the NRE
     */
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

    /**
     * linear cost coefficients for presence of bb at position pos,
     * when chains are aligned by the RE; mixture of chain lengths
     * @param pos position from RE
     * @param bb disaccharide
     * @return linear cost coefficients for presence of bb at position pos,
     * when chains are aligned by the RE
     */
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

    /**
     * linear cost coefficients for bb at position pos
     * when all BKHS have same length
     * @param pos position from NRE
     * @param bb disaccharide
     * @return linear cost coefficients for bb at position pos
     * when all BKHS have same length
     */
    double[] getCostCoeff(int pos, int bb) {
        double[] res = new double[sp.N];
        for (int s = 0; s < sp.N; s++) {
            if (sp.seq[s][pos] == bb) {
                res[s] = 1.;
            }
        }
        return res;
    }

    /**
     * lower bound for proportion of bb at position pos
     * when chains are aligned by NRE; mixture of chain lengths
     * @param pos position from the NRE
     * @param bb disaccharide
     * @return lower bound for proportion of bb at position pos
     * when chains are aligned by NRE; mixture of chain lengths
     */
    double getLowerNRE(int pos, int bb) {
        double[] c = getCostCoeffNRE(pos, bb);
        Simplex sp = new Simplex(sp1, c);
        double res = sp.finalCost;
        if (res > 1.) {
            res = 1.;
        }
        if (res < 0.) {
            res = 0.;
        }
        return res;
    }

    /**
     * lower bound for proportion of bb at position pos
     * when chains are aligned by RE; mixture of chain lengths
     * @param pos position from RE
     * @param bb disaccharide
     * @return lower bound for proportion of bb at position pos
     * when chains are aligned by RE; mixture of chain lengths
     */
    double getLowerRE(int pos, int bb) {
        double[] c = getCostCoeffRE(pos, bb);
        Simplex sp = new Simplex(sp1, c);
        double res = sp.finalCost;
        if (res > 1.) {
            res = 1.;
        }
        if (res < 0.) {
            res = 0.;
        }
        return res;
    }

    /**
     * upper bound for proportion of bb at position pos
     * when chains are aligned by NRE; mixture of chain lengths
     * @param pos position from NRE
     * @param bb disaccharide
     * @return upper bound for proportion of bb at position pos
     * when chains are aligned by NRE; mixture of chain lengths
     */
    double getUpperNRE(int pos, int bb) {
        double[] c = getCostCoeffNRE(pos, bb);
        for (int i = 0; i < c.length; i++) {
            c[i] *= -1.;
        }
        Simplex sp = new Simplex(sp1, c);
        double res = -sp.finalCost;
        if (res > 1.) {
            res = 1.;
        }
        if (res < 0.) {
            res = 0.;
        }
        return res;
    }

    /**
     * upper bound for proportion of bb at position pos
     * when chains are aligned by RE; mixture of chain lengths
     * @param pos position from RE
     * @param bb disaccharide
     * @return upper bound for proportion of bb at position pos
     * when chains are aligned by RE; mixture of chain lengths
     */
    double getUpperRE(int pos, int bb) {
        double[] c = getCostCoeffRE(pos, bb);
        for (int i = 0; i < c.length; i++) {
            c[i] *= -1.;
        }
        Simplex sp = new Simplex(sp1, c);
        double res = -sp.finalCost;
        if (res > 1.) {
            res = 1.;
        }
        if (res < 0.) {
            res = 0.;
        }
        return res;
    }

    /**
     * lower bound for proportion of bb at position pos
     * when all BKHS chains have same length
     * @param pos position from the NRE
     * @param bb disaccharide
     * @return lower bound for proportion of bb at position pos
     * when all BKHS chains have same length
     */
    double getLower(int pos, int bb) {
        double[] c = getCostCoeff(pos, bb);
        Simplex sp = new Simplex(sp1, c);
        double res = sp.finalCost;
        if (res > 1.) {
            res = 1.;
        }
        if (res < 0.) {
            res = 0.;
        }
        return res;
    }

    /**
     * upper bound for proportion of bb at position pos
     * when all BKHS chains have same length
     * @param pos position from NRE
     * @param bb disaccharide
     * @return upper bound for proportion of bb at position pos
     * when all BKHS chains have same length
     */
    double getUpper(int pos, int bb) {
        double[] c = getCostCoeff(pos, bb);
        for (int i = 0; i < c.length; i++) {
            c[i] *= -1.;
        }
        Simplex sp = new Simplex(sp1, c);
        double res = -sp.finalCost;
        if (res > 1.) {
            res = 1.;
        }
        if (res < 0.) {
            res = 0.;
        }
        return res;
    }

    /**
     * saves bounds for S/U proportions at each position
     * when all BKHS chains have same length
     * @param file output file
     */
    void saveBounds(String file) {
        List<String> v = new ArrayList<>();
        String s = "pos";
        for (int i = 0; i < sp.m; i++) {
            s += "\tlower" + lab[i] + "\tupper" + lab[i];
        }
        v.add(s);
        for (int pos = 0; pos < sp.n; pos++) {
            s = Integer.toString(pos + 1);
            for (int j = 0; j < sp.m; j++) {
                s += "\t" + getLower(pos, j);
                s += "\t" + getUpper(pos, j);
            }
            v.add(s);
            System.err.println(s);
        }
        Utils.saveFile(v, file);
    }

    /**
     * saves bounds of S/U proportions at each position from the NRE
     * for a mixture of BKHS chain lengths
     * @param file output file
     */
    void saveBoundsMixSpNRE(String file) {
        List<String> v = new ArrayList<>();
        String s = "posNRE";
        for (int i = 0; i < msp.m; i++) {
            s += "\tlower" + lab[i] + "\tupper" + lab[i];
        }
        v.add(s);
        for (int pos = 0; pos < msp.lmax; pos++) {
            s = Integer.toString(pos + 1);
            for (int j = 0; j < msp.m; j++) {
                s += "\t" + (getLowerNRE(pos, j) / msp.lengthAtLeast[pos]);
                s += "\t" + (getUpperNRE(pos, j) / msp.lengthAtLeast[pos]);
            }
            v.add(s);
            System.err.println(s);
        }
        Utils.saveFile(v, file);
    }

    /**
     * saves bounds of S/U proportions at each position from the RE
     * for a mixture of BKHS chain lengths
     * @param file output file
     */
    void saveBoundsMixSpRE(String file) {
        List<String> v = new ArrayList<>();
        String s = "posRE";
        for (int i = 0; i < msp.m; i++) {
            s += "\tlower" + lab[i] + "\tupper" + lab[i];
        }
        v.add(s);
        for (int pos = 0; pos < msp.lmax; pos++) {
            s = Integer.toString(pos + 1);
            for (int j = 0; j < msp.m; j++) {
                s += "\t" + (getLowerRE(pos, j) / msp.lengthAtLeast[pos]);
                s += "\t" + (getUpperRE(pos, j) / msp.lengthAtLeast[pos]);
            }
            v.add(s);
            System.err.println(s);
        }
        Utils.saveFile(v, file);
    }

    /**
     * makes and saves S/U bounds at each position
     * when all BKHS chains have same length
     * @param inDir input directory (ends with "\\")
     * @param outDir output directory (ends with "\\")
     */
    public static void bounds(String inDir, String outDir) {
        String[] lab = new String[2];
        lab[0] = "U";
        lab[1] = "S";
        GradientWidth gw = new GradientWidth(2, 16, lab, inDir);
        gw.saveBounds(outDir + "GradientWidth.res");
    }

    /**
     * for supplementary material; bounds for different values of chain length n
     * @param inDir input directory (ends with "\\")
     * @param outDir output directory (ends with "\\")
     */
    public static void boundsAndN(String inDir, String outDir) {
        String[] lab = new String[2];
        lab[0] = "U";
        lab[1] = "S";
        for (int n = 13; n <= 18; n++) {
            GradientWidth gw = new GradientWidth(2, n, lab, inDir);
            gw.saveBounds(outDir + "GradientWidth.N" + n + ".res");
        }
    }

    /**
     * bounds of S/U proportions at each position for mixture of BKHS chain lengths
     */
    public static void defaultMixSpec(String inDir, String outDir) {
        String[] lab = new String[2];
        lab[0] = "U";
        lab[1] = "S";
        GradientWidth gw = new GradientWidth(lab, 10, 20, 3.5, 16., inDir);
        gw.saveBoundsMixSpNRE(outDir + "BoundsMixSpNRE.res");
        gw.saveBoundsMixSpRE(outDir + "BoundsMixSpRE.res");
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        String outDir = "output\\";
        bounds(inDir, outDir);
        boundsAndN(inDir, outDir);
        defaultMixSpec(inDir, outDir);
    }
}