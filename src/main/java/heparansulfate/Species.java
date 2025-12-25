package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * enumeration of all possible sequences of length n with m disaccharides
 */
public class Species {
    /**
     * chain length
     */
    int n = 0;
    /**
     * number of building blocks
     */
    int m = 0;
    /**
     * number of species = m^n
     */
    int N = 0;
    /**
     * seq[i] is the sequence of species i
     */
    int[][] seq = null;

    /**
     * enumeration of all possible sequences of length n with m disaccharides
     * @param m number of disaccharides
     * @param n BKHS chain length
     */
    public Species(int m, int n) {
        this.m = m;
        this.n = n;
        N = (int) Math.pow((double) m, (double) n);
        seq = new int[N][n];
        int[] buff = new int[n];
        int pos = 0;
        for (int s = 1; s < N; s++) {
            pos = 0;
            boolean kg = true;
            while (kg) {
                buff[pos]++;
                if (buff[pos] < m) {
                    kg = false;
                } else {
                    buff[pos] = 0;
                    pos++;
                }
            }
            for (int i = 0; i < n; i++) {
                seq[s][i] = buff[i];
            }
        }
    }

    /**
     * linear constraint for homogeneity of disaccharide composition along chains
     * @param bbs disaccharide set
     * @return linear constraint for homogeneity of disaccharide composition
     * along chains
     */
    LinEqCons getHomogeneityLEC(BBSet bbs) {
        double[][] A = new double[(n - 1) * (m - 1)][N];
        double[] b = new double[(n - 1) * (m - 1)];
        String[] tp = new String[(n - 1) * (m - 1)];
        int r = -1;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < m - 1; j++) {
                r++;
                b[r] = bbs.rho[j];
                tp[r] = "hom";
                for (int s = 0; s < N; s++) {
                    if (seq[s][i] == j) {
                        A[r][s] = 1.;
                    }
                }
            }
        }
        return new LinEqCons(A, b, tp);
    }

    /**
     * linear constraint for one heparinase digest (fragments lengths f[])
     * @param f distribution of fragment lengths
     * @param cs cleavage specificities of one heparinase
     * @param bbs disaccharide set
     * @return linear constraint for one heparinase digest
     */
    LinEqCons getFragLEC(double[] f, CSpec cs, BBSet bbs) {
        double[][] A = new double[f.length][N];
        double[] b = new double[f.length];
        String[] lab = new String[f.length];
        for (int i = 0; i < b.length; i++) {
            b[i] = f[i];
            lab[i] = "frag";
        }
        int lm = f.length;
        for (int l = 0; l < lm - 1; l++) {
            for (int s = 0; s < N; s++) {
                A[l][s] = getFragContrib(s, l + 1, cs);
            }
        }
        for (int s = 0; s < N; s++) {
            A[lm - 1][s] = 0.;
            for (int l = lm - 1; l < n - 1; l++) {
                A[lm - 1][s] += getFragContrib(s, l + 1, cs);
            }
        }
        double denom = 0.;
        for (int i = 0; i < bbs.m; i++) {
            denom += bbs.rho[i] * cs.c[i];
        }
        denom *= (double)(n - 1);
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                A[i][j] /= denom;
            }
        }
        return new LinEqCons(A, b, lab);
    }

    /**
     * contribution of species s to fragments of length ll
     * @param s species (sequence) index
     * @param ll fragment length
     * @param cs cleavage specificities
     * @return contribution of species s to fragments of length ll
     */
    double getFragContrib(int s, int ll, CSpec cs) {
        double res = 0.;
        if (ll > n) {
            return res;
        }
        for (int i = 1; i < n - ll; i++) {
            double t = cs.c[seq[s][i]];
            for (int j = 1; j <= ll - 1; j++) {
                t *= (1. - cs.c[seq[s][i + j]]);
            }
            t *= cs.c[seq[s][i + ll]];
            res += t;
        }
        double t = cs.c[seq[s][n - ll]];
        for (int j = n - ll + 1; j < n; j++) {
            t *= (1. - cs.c[seq[s][j]]);
        }
        res += t;
        return res;
    }

    /**
     * linear constraint for overall disaccharide composition
     * @param rho overall disaccharide composition
     * @return linear constraint for overall disaccharide composition
     */
    LinEqCons getCompLEC(double[] rho) {
        int m = rho.length;
        double[][] A = new double[m][N];
        double[] b = new double[m];
        String[] lab = new String[m];
        for (int i = 0; i < m; i++) {
            b[i] = rho[i];
            lab[i] = "comp";
            for (int s = 0; s < N; s++) {
                for (int j = 1; j < n; j++) {
                    if (seq[s][j] == i) {
                        A[i][s] += 1.;
                    }
                }
                A[i][s] /= (double)(n - 1);
            }
        }
        return new LinEqCons(A, b, lab);
    }

    /**
     * linear constraint for normalization (sum to 1)
     * @return linear constraint for normalization (sum to 1)
     */
    LinEqCons getNormLEC() {
        double[][] A = new double[1][N];
        double[] b = new double[1];
        String[] lab = new String[1];
        b[0] = 1.;
        lab[0] = "norm";
        for (int i = 0; i < N; i++) {
            A[0][i] = 1.;
        }
        return new LinEqCons(A, b, lab);
    }

    /**
     * vector of individual species abundances; this method is based on a
     * N&I model and utilized to numerically check linear constraints
     * @param gamma matrix Gamma in model N&I
     * @return vector of individual species abundances under N&I with Gamma
     */
    double[] getP(double[][] gamma) {
        double[] res = new double[N];
        for (int s = 0; s < N; s++) {
            res[s] = 1.;
            for (int i = 0; i < n; i++) {
                res[s] *= gamma[i][seq[s][i]];
            }
        }
        return res;
    }

    /**
     * complete set of linear constraints (sum-to-1, overall disaccharide
     * composition and fragment length distributions)
     * @param bbs disaccharide set
     * @param csFile files of cleavage specificities
     * @param fragFile files of fragment length distributions
     * @return complete set of linear constraints (sum-to-1, overall
     * disaccharide composition and fragment length distributions)
     */
    LinEqCons getCompleteLEC(BBSet bbs, String[] csFile, String[] fragFile) {
        List<LinEqCons> v = new ArrayList<>();
        v.add(getNormLEC());
        v.add(LinEqCons.removeLastRow(getCompLEC(bbs.rho)));
        for (int i = 0; i < csFile.length; i++) {
            CSpec cs = new CSpec(csFile[i], bbs);
            double[] f = loadFragAbund(fragFile[i]);
            v.add(LinEqCons.removeLastRow(getFragLEC(f, cs, bbs)));
        }
        return new LinEqCons(v);
    }

    /**
     * loads matrix Gamma of model N&I from a file;
     * used for testing of linear constraints
     * @param file ASCII file with one header row: position \t u \t s
     * @return matrix Gamma of a N&I model
     */
    public static double[][] loadGamma(String file) {
        List<String> v = Utils.loadFileNoheader(file);
        StringTokenizer st = new StringTokenizer(v.get(0), "\t");
        int m = st.countTokens() - 1;
        int n = v.size();
        double[][] res = new double[n][m];
        for (int i = 0; i < n; i++) {
            st = new StringTokenizer(v.get(i), "\t");
            st.nextToken();
            int j = -1;
            while (st.hasMoreTokens()) {
                j++;
                res[i][j] = Double.parseDouble(st.nextToken());
            }
        }
        return res;
    }

    /**
     * loads fragment length abundances for one heparinase digest
     * @param file ASCII file with header row: length \t abundance;
     * note that lengths are assumed to be in ascending order, the last row
     * ll represents abundances of lengths at least ll and abundances must sum to 1
     * @return fragment length abundances for one heparinase digest
     */
    public static double[] loadFragAbund(String file) {
        List<String> v = Utils.loadFileNoheader(file);
        StringTokenizer st = new StringTokenizer(v.get(0), "\t");
        int n = v.size();
        double[] res = new double[n];
        for (int i = 0; i < n; i++) {
            st = new StringTokenizer(v.get(i), "\t");
            st.nextToken();
            res[i] = Double.parseDouble(st.nextToken());
        }
        return res;
    }

    /**
     * testing of the sum-to-1 constraint, with a N&I model
     * @param inDir input directory (ends with "\\")
     */
    public static void testNormCons(String inDir) {
        System.out.println("testNormCons");
        Species sp = new Species(2, 16);
        double[] p = sp.getP(loadGamma(inDir + "GammaExample.txt"));
        LinEqCons lec = sp.getNormLEC();
        double[] b = MatrixOp.multMatVec(lec.A, p);
        for (int i = 0; i < b.length; i++) {
            System.out.println(b[i] + "\t" + lec.b[i]);
        }
        System.out.println("*****");
    }

    /**
     * testing of the constraint for overall disaccharide composition,
     * with a N&I model
     * @param inDir input directory (ends with "\\")
     */
    public static void testCompCons(String inDir) {
        System.out.println("testCompCons");
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        Species sp = new Species(2, 16);
        double[] p = sp.getP(loadGamma(inDir + "GammaExample.txt"));
        LinEqCons lec = sp.getCompLEC(bbs.rho);
        double[] b = MatrixOp.multMatVec(lec.A, p);
        for (int i = 0; i < b.length; i++) {
            System.out.println(b[i] + "\t" + lec.b[i]);
        }
        System.out.println("*****");
    }

    /**
     * testing of the constraints for fragment length distributions,
     * with a N&I model; note that since N&I cannot reproduce fragment
     * length distributions, values of constraints are calculated
     * from the N&I model
     * @param inDir input directory (ends with "\\")
     */
    public static void testFragCons(String inDir) {
        System.out.println("testFragCons");
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        Species sp = new Species(2, 16);
        double[][] gamma = loadGamma(inDir + "GammaExample.txt");
        double[] p = sp.getP(gamma);
        System.out.println("hepI (based on N&I model)");
        CSpec cs = new CSpec(inDir + "US.hepI.txt", bbs);
        NIModel nim = new NIModel(bbs, cs, gamma, 11);
        double[] f = nim.h;
        LinEqCons lec = sp.getFragLEC(f, cs, bbs);
        double[] b = MatrixOp.multMatVec(lec.A, p);
        for (int i = 0; i < b.length; i++) {
            System.out.println(b[i] + "\t" + lec.b[i]);
        }
        System.out.println("*****");
        System.out.println("hepIII (based on N&I model)");
        cs = new CSpec(inDir + "US.hepIII.txt", bbs);
        nim = new NIModel(bbs, cs, gamma, 6);
        f = nim.h;
        lec = sp.getFragLEC(f, cs, bbs);
        b = MatrixOp.multMatVec(lec.A, p);
        for (int i = 0; i < b.length; i++) {
            System.out.println(b[i] + "\t" + lec.b[i]);
        }
        System.out.println("*****");
    }

    /**
     * testing of homogeneity constraint with model H&I
     * @param inDir input directory (ends with "\\")
     */
    public static void testHomogeneityCons(String inDir) {
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        Species sp = new Species(2, 16);
        LinEqCons lec = sp.getHomogeneityLEC(bbs);
        SimplexPhaseI sp1 = new SimplexPhaseI(lec.A, lec.b);
        System.out.println("infeasibility = " + sp1.finalCost);
        System.out.println("*****");
        int N = sp.N;
        double[] p = new double[N];
        for (int s = 0; s < N; s++) {
            p[s] = 1.;
            for (int i = 0; i < sp.n; i++) {
                p[s] *= bbs.rho[sp.seq[s][i]];
            }
        }
        double[] b = MatrixOp.multMatVec(lec.A, p);
        for (int i = 0; i < b.length; i++) {
            System.out.println(b[i] + "\t" + lec.b[i]);
        }
    }

    /**
     * testing of computed linear constraints with N&I or H&I model
     * @param inDir input directory (ends with "\\")
     */
    public static void testing(String inDir) {
        testNormCons(inDir);
        testCompCons(inDir);
        testFragCons(inDir);
        testHomogeneityCons(inDir);
    }

    /**
     * for testing
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        testing(inDir);
    }
}
