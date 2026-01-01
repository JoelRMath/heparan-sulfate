package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Enumeration of sequences (species) with a specified chain length distribution.
 * The model of the chain length distribution is obtained by projection.
 */
public class MixSpecies {
    /** Number of disaccharides. */
    int m = 2;
    /** Smallest chain length. */
    int lmin = 0;
    /** Longest chain length. */
    int lmax = 0;
    /** Number of unique chain lengths. */
    int nc = 0;
    /** Array of chain lengths. */
    int[] len = null;
    /** Chain length distribution (weights). */
    double[] w = null;
    /** Total number of species across all lengths. */
    int N = 0;
    /** Enumerated sequences: {@code seq[i]} represents sequence index {@code i}. */
    int[][] seq = null;
    /** * Based on {@code w}, {@code lengthAtLeast[l]} is the proportion of chains 
     * having length at least {@code l+1}. 
     */
    double[] lengthAtLeast = null;

    /**
     * Mixture model of BKHS chain lengths and sequence enumeration.
     * @param lmin Smallest chain length.
     * @param lmax Largest chain length.
     * @param sigma Spread (standard deviation) of chain lengths.
     * @param mu Average chain length.
     */
    public MixSpecies(int lmin, int lmax, double sigma, double mu) {
        this.lmin = lmin;
        this.lmax = lmax;
        nc = lmax - lmin + 1;
        w = getW(lmin, lmax, sigma, mu);
        lengthAtLeast = new double[lmax];
        for (int l = 0; l < lmax; l++) {
            int j = l + 1 - lmin;
            if (j < 0) {
                j = 0;
            }
            for (int i = j; i < w.length; i++) {
                lengthAtLeast[l] += w[i];
            }
        }
        len = new int[nc];
        for (int i = 0; i < nc; i++) {
            len[i] = lmin + i;
        }
        N = 0;
        for (int l = 0; l < nc; l++) {
            N += (int) Math.pow((double) m, (double) len[l]);
        }
        seq = new int[N][];
        int s = -1;
        for (int l = 0; l < nc; l++) {
            int max = (int) Math.pow((double) m, (double) len[l]);
            for (int i = 0; i < max; i++) {
                s++;
                seq[s] = new int[len[l]];
            }
        }
        s = -1;
        for (int l = 0; l < nc; l++) {
            Species sp = new Species(m, len[l]);
            for (int i = 0; i < sp.N; i++) {
                s++;
                for (int j = 0; j < sp.n; j++) {
                    seq[s][j] = sp.seq[i][j];
                }
            }
        }
        MatrixOp.printVec(w);
    }

    /**
     * Combines constraints on the vector of species abundances: normalization, 
     * composition, and fragment length distributions after digestion.
     * @param bbs Disaccharide set with composition in {@code bbs.rho}.
     * @param csFile Files containing cleavage specificities.
     * @param fragFile Files containing distributions of fragment lengths.
     * @return Combined linear equality constraints.
     */
    LinEqCons getCompleteLEC(BBSet bbs, String[] csFile, String[] fragFile) {
        List<LinEqCons> v = new ArrayList<>();
        v.add(getNormLEC());
        v.add(LinEqCons.removeLastRow(getRhoLEC(bbs)));
        for (int i = 0; i < csFile.length; i++) {
            CSpec cs = new CSpec(csFile[i], bbs);
            double[] f = Species.loadFragAbund(fragFile[i]);
            v.add(LinEqCons.removeLastRow(getFragLEC(f, cs, bbs)));
        }
        v.add(LinEqCons.removeLastRow(getWLEC()));
        return new LinEqCons(v);
    }

    /**
     * Constraint for homogeneity of composition with chain length combined 
     * with digest constraints.
     * @param bbs Disaccharide set.
     * @param csFile Files containing cleavage specificities.
     * @param fragFile Files containing distributions of fragment lengths.
     * @return Combined linear equality constraints for homogeneity feasibility.
     */
    LinEqCons getRhoLHomFeasLEC(BBSet bbs, String[] csFile, String[] fragFile) {
        List<LinEqCons> v = new ArrayList<>();
        v.add(getNormLEC());
        for (int i = 0; i < csFile.length; i++) {
            CSpec cs = new CSpec(csFile[i], bbs);
            double[] f = Species.loadFragAbund(fragFile[i]);
            v.add(LinEqCons.removeLastRow(getFragLEC(f, cs, bbs)));
        }
        v.add(LinEqCons.removeLastRow(getWLEC()));
        v.add(getRhoLHomLEC(bbs));
        return new LinEqCons(v);
    }

    /**
     * Constraint for homogeneity of composition with chain length.
     * @param bbs Disaccharide set.
     * @return Linear equality constraint for homogeneity.
     */
    LinEqCons getRhoLHomLEC(BBSet bbs) {
        double[][] A = new double[nc][N];
        double[] b = new double[nc];
        String[] tp = new String[nc];
        for (int i = 0; i < nc; i++) {
            b[i] = bbs.rho[1] * w[i];
            tp[i] = "rhol";
            for (int s = 0; s < N; s++) {
                if (seq[s].length == len[i]) {
                    for (int k = 1; k < seq[s].length; k++) {
                        if (seq[s][k] == 1) {
                            A[i][s] += 1.;
                        }
                    }
                    A[i][s] /= (double)(len[i] - 1);
                }
            }
        }
        return new LinEqCons(A, b, tp);
    }

    /**
     * Linear equality constraint on species abundances for chain length distribution.
     * @return Linear equality constraint for the distribution {@code w}.
     */
    LinEqCons getWLEC() {
        double[][] A = new double[nc][N];
        double[] b = new double[nc];
        String[] tp = new String[nc];
        for (int i = 0; i < nc; i++) {
            b[i] = w[i];
            tp[i] = "W";
            for (int s = 0; s < N; s++) {
                if (seq[s].length == len[i]) {
                    A[i][s] = 1.;
                }
            }
        }
        return new LinEqCons(A, b, tp);
    }

    /**
     * Sum-to-1 (normalization) constraint on the vector of species abundances.
     * @return Normalization linear equality constraint.
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
     * Linear equality constraint for the distribution {@code f} of fragment lengths 
     * after digestion by an enzyme.
     * @param f Distribution of fragment lengths.
     * @param cs Cleavage specificities.
     * @param bbs Disaccharide set (composition in {@code bbs.rho}).
     * @return Linear equality constraint for fragment distribution.
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
            for (int l = lm - 1; l < seq[s].length - 1; l++) {
                A[lm - 1][s] += getFragContrib(s, l + 1, cs);
            }
        }
        double denom = 0.;
        for (int i = 0; i < bbs.m; i++) {
            denom += bbs.rho[i] * cs.c[i];
        }
        double avL = 0.;
        for (int l = 0; l < nc; l++) {
            avL += w[l] * (double) len[l];
        }
        denom *= (avL - 1.);
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                A[i][j] /= denom;
            }
        }
        return new LinEqCons(A, b, lab);
    }

    /**
     * Contribution of species {@code s} to the expected number of fragments of 
     * length {@code ll} after digestion.
     * @param s Species index.
     * @param ll Fragment length.
     * @param cs Cleavage specificities.
     * @return Expected fragment contribution.
     */
    double getFragContrib(int s, int ll, CSpec cs) {
        int n = seq[s].length;
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
     * Linear equality constraint on species abundances for overall disaccharide 
     * composition between positions 2 and {@code n}.
     * @param bbs Set of disaccharides with composition in {@code bbs.rho}.
     * @return Composition linear equality constraint.
     */
    LinEqCons getRhoLEC(BBSet bbs) {
        double[][] A = new double[m][N];
        double[] b = new double[m];
        String[] tp = new String[m];
        double avL = 0.;
        for (int l = 0; l < nc; l++) {
            avL += w[l] * (double) len[l];
        }
        for (int j = 0; j < m; j++) {
            b[j] = bbs.rho[j];
            tp[j] = "rho";
            for (int s = 0; s < N; s++) {
                for (int k = 1; k < seq[s].length; k++) {
                    if (seq[s][k] == j) {
                        A[j][s] += 1.;
                    }
                }
                A[j][s] /= (avL - 1.);
            }
        }
        return new LinEqCons(A, b, tp);
    }

    /**
     * Tests that homogeneity of composition across chain lengths is feasible 
     * after incorporating digest constraints.
     * <p>Note: Requires large memory (e.g., {@code -Xmx8000M}).</p>
     * @param inDir Input directory.
     */
    public static void testRhoLHomFeas(String inDir) {
        int lmin = 10;
        int lmax = 20;
        double sigma = 3.5;
        double mu = 16.;
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        MixSpecies msp = new MixSpecies(lmin, lmax, sigma, mu);
        String[] csFile = new String[2];
        String[] fragFile = new String[2];
        csFile[0] = inDir + "US.hepI.txt";
        csFile[1] = inDir + "US.hepIII.txt";
        fragFile[0] = inDir + "hepI.f.txt";
        fragFile[1] = inDir + "hepIII.f.txt";
        LinEqCons lec = msp.getRhoLHomFeasLEC(bbs, csFile, fragFile);
        SimplexPhaseI sp1 = new SimplexPhaseI(lec.A, lec.b);
        System.err.println(sp1.finalCost);
    }

    /**
     * Computes the chain length distribution. A Gaussian ({@code mu}, {@code sigma}) 
     * restricted between {@code lmin} and {@code lmax} is used as an approximation, 
     * followed by projection on the probability simplex to preserve {@code mu}.
     * @param lmin Smallest chain length.
     * @param lmax Longest chain length.
     * @param sigma Standard deviation of Gaussian approximation.
     * @param mu Desired average chain length.
     * @return Chain length distribution between {@code lmin} and {@code lmax} with mean {@code mu}.
     */
    public static double[] getW(int lmin, int lmax, double sigma, double mu) {
        int nc = lmax - lmin + 1;
        double[] res = new double[nc];
        double[] t = new double[nc];
        for (int l = lmin; l <= lmax; l++) {
            int i = l - lmin;
            t[i] = 0.;
            double dx = 0.001;
            double x = (double) l;
            while (x < (double)(l + 1)) {
                x += dx;
                double fx = Math.exp(-0.5 * (x - mu) * (x - mu) / (sigma * sigma))
                / (sigma * Math.sqrt(2. * Math.PI));
                t[i] += fx * dx;
            }
        }
        double sum = 0.;
        for (int i = 0; i < nc; i++) {
            sum += t[i];
        }
        for (int i = 0; i < nc; i++) {
            res[i] = t[i] / sum;
        }
        double[][] A = new double[nc + 4][nc];
        double[] b = new double[nc + 4];
        int r = -1;
        for (int i = 0; i < nc; i++) {
            b[i] = 0.;
            r++;
            A[r][i] = -1.;
        }
        r++;
        b[r] = 1.;
        for (int i = 0; i < nc; i++) {
            A[r][i] = 1.;
        }
        r++;
        b[r] = -1.;
        for (int i = 0; i < nc; i++) {
            A[r][i] = -1.;
        }
        r++;
        b[r] = mu;
        for (int i = 0; i < nc; i++) {
            A[r][i] = (double)(lmin + i);
        }
        r++;
        b[r] = -mu;
        for (int i = 0; i < nc; i++) {
            A[r][i] = -(double)(lmin + i);
        }
        ProjOnPolyHSet proj = new ProjOnPolyHSet(res, A, b, new Random());
        sum = 0.;
        for (int i = 0; i < res.length; i++) {
            res[i] = proj.optimum[i];
            if (res[i] < 0.) {
                res[i] = 0.;
            }
            sum += res[i];
        }
        double av = 0.;
        for (int i = 0; i < res.length; i++) {
            res[i] /= sum;
            av += res[i] * (double)(lmin + i);
        }
        System.err.println(av);
        return res;
    }

    /**
     * Saves a distribution of chain lengths to a file.
     * @param lmin Smallest length.
     * @param lmax Largest length.
     * @param sigma Spread parameter.
     * @param mu Mean length.
     * @param filePath Path to the output file.
     */
    public void saveFullDistribution(int lmin, int lmax, double sigma, double mu, String filePath) {
        List<String> lines = new ArrayList<>();
        lines.add("l\tw");
        double[] weights = getW(lmin, lmax, sigma, mu); 
        for (int i = 0; i < weights.length; i++) {
            lines.add((lmin + i) + "\t" + weights[i]);
        }
        Utils.saveFile(lines, filePath);
    }

    /**
     * Static wrapper to compute and save chain length distribution.
     * @param lmin Smallest length.
     * @param lmax Largest length.
     * @param sigma Spread parameter.
     * @param mu Mean length.
     * @param filePath Path to the output file.
     */
    public static void saveDistribution(int lmin, int lmax, double sigma, double mu, String filePath) {
        MixSpecies msp = new MixSpecies(lmin, lmax, sigma, mu);
        msp.saveFullDistribution(lmin, lmax, sigma, mu, filePath);
    }

    /**
     * Main entry point for testing feasibility.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        testRhoLHomFeas(inDir);
    }
}