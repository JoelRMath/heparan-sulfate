package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * Enumeration of all possible sequences of length {@code n} with {@code m} disaccharides.
 * This class facilitates the generation of linear equality constraints based on 
 * molecular species sequences.
 */
public class Species {
    /**
     * Chain length.
     */
    int n = 0;
    /**
     * Number of building blocks.
     */
    int m = 0;
    /**
     * Total number of unique species, calculated as {@code m^n}.
     */
    int N = 0;
    /**
     * Array of sequences where {@code seq[i]} represents the sequence of species {@code i}.
     */
    int[][] seq = null;

    /**
     * Enumerates all possible sequences of length {@code n} with {@code m} disaccharides.
     * @param m Number of disaccharides.
     * @param n BKHS chain length.
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
     * Creates linear constraints for the homogeneity of disaccharide composition along chains.
     * @param bbs Disaccharide set.
     * @return Linear constraint for disaccharide homogeneity along chains.
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
     * Creates linear constraints for one heparinase digest based on fragment length distribution.
     * @param f Distribution of fragment lengths.
     * @param cs Cleavage specificities of the heparinase enzyme.
     * @param bbs Disaccharide set.
     * @return Linear constraint for the heparinase digest.
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
     * Calculates the contribution of species {@code s} to fragments of length {@code ll}.
     * 
     * @param s Species (sequence) index.
     * @param ll Fragment length.
     * @param cs Cleavage specificities.
     * @return Contribution of species {@code s} to fragments of length {@code ll}.
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
     * Creates linear constraints for the overall disaccharide composition.
     * @param rho Overall disaccharide composition vector.
     * @return Linear constraint for overall disaccharide composition.
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
     * Creates linear constraints for normalization (ensuring abundance sums to 1).
     * @return Linear constraint for normalization.
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
     * Calculates the vector of individual species abundances based on a Nonhomogeneity 
     * and Independence (N&amp;I) model. Utilized for numerical checking of linear constraints.
     * @param gamma Matrix {@code Gamma} in model N&amp;I.
     * @return Vector of species abundances under the N&amp;I model with {@code Gamma}.
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
     * Combines all linear constraints including sum-to-1, disaccharide composition, 
     * and fragment length distributions.
     * @param bbs Disaccharide set.
     * @param csFile Array of files specifying cleavage specificities.
     * @param fragFile Array of files specifying fragment length distributions.
     * @return Complete set of merged linear equality constraints.
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
     * Loads the matrix {@code Gamma} of a N&amp;I model from a tab-delimited ASCII file.
     * @param file ASCII file with header row: {@code position \t u \t s}.
     * @return Matrix {@code Gamma} representing the N&amp;I model.
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
     * Loads fragment length abundances for a heparinase digest from a tab-delimited ASCII file.
     * @param file ASCII file with header row: {@code length \t abundance}.
     * @return Array of fragment length abundances.
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
     * Performs testing of the sum-to-1 constraint using a N&amp;I model.
     * @param inDir Input directory.
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
     * Performs testing of the overall disaccharide composition constraint using a N&amp;I model.
     * @param inDir Input directory.
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
     * Performs testing of the fragment length distribution constraints using a N&amp;I model.
     * @param inDir Input directory.
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
     * Performs testing of the homogeneity constraint using a model with overall 
     * disaccharide composition.
     * @param inDir Input directory.
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
     * Runs all linear constraint tests using N&amp;I or overall disaccharide composition models.
     * @param inDir Input directory.
     */
    public static void testing(String inDir) {
        testNormCons(inDir);
        testCompCons(inDir);
        testFragCons(inDir);
        testHomogeneityCons(inDir);
    }

    /**
     * Main entry point for the test runner.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        testing(inDir);
    }
}