package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * species (enumeration of sequences) with specified chain length distribution;
 * the model of chain length distribution is obtained by projection
 */
public class MixSpecies {
    /**
     * number of disaccharides
     */
    int m = 2;
    /**
     * smallest chain length
     */
    int lmin = 0;
    /**
     * longest chain length
     */
    int lmax = 0;
    /**
     * number of chain lengths
     */
    int nc = 0;
    /**
     * chain lengths
     */
    int[] len = null;
    /**
     * chain length distribution
     */
    double[] w = null;
    /**
     * number of species
     */
    int N = 0;
    /**
     * sequences: seq[i] = sequence # i
     */
    int[][] seq = null;
    /**
     * based on this.w, lengthAtLeast[l] is the proportion of chains having
     * length at least l+1
     */
    double[] lengthAtLeast = null;

    /**
     * mixture model of BKHS chain lengths; sequence enumeration
     * @param lmin smallest chain length
     * @param lmax largest chain length
     * @param sigma spread of chain lengths
     * @param mu average chain length
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
     * combines constraints on vector of species abundances: normalization, composition
     * and fragment length distributions (fragFile[]) after digestion (csFile[])
     * @param bbs disaccharide set, with composition in bbs.rho
     * @param csFile files containing cleavage specificities
     * @param fragFile files containing distributions of fragment lengths
     * @return combined constraints on vector of species abundances: normalization,
     * composition and fragment length distributions (fragFile[]) after digestion (csFile[])
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
     * constraint for homogeneity of composition with chain length combined to
     * digest constraints
     * @param bbs disaccharide set
     * @param csFile files containing cleavage specificities
     * @param fragFile files containing distributions of fragment lengths
     * @return constraint for homogeneity of composition with chain length
     * combined to digest constraints
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
     * constraint for homogeneity of composition with chain length
     * @param bbs disaccharide set
     * @return constraint for homogeneity of composition with chain length
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
     * linear equality constraint on species abundances for chain length distribution
     * @return linear equality constraint on species abundances for
     * chain length distribution
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
     * sum-to-1 constraint on vector of species abundances
     * @return sum-to-1 constraint on vector of species abundances
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
     * linear equality constraint on the vector of species abundances for distribution
     * f of fragment length after digestion by enzyme (cleavage specificities in cs.c)
     * @param f distribution of fragment length
     * @param cs cleavage specificities
     * @param bbs disaccharide set (composition in bbs.rho)
     * @return linear equality constraint on the vector of species abundances for
     * distribution f of fragment length after digestion by enzyme
     * (cleavage specificities in cs.c)
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
     * contribution of species s to expected number of fragments of length
     * ll after digestion by enzyme with cleavage specificities given by cs.c
     * @param s species index
     * @param ll fragment length
     * @param cs cleavage specificities
     * @return contribution of species s to expected number of fragments of length
     * ll after digestion by enzyme with cleavage specificities given by cs.c
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
     * linear equality constraint on species abundances for overall disaccharide
     * composition between positions 2 and n
     * @param bbs set of disaccharides with composition in bbs.rho
     * @return linear equality constraint on species abundances for overall
     * disaccharide composition between positions 2 and n
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
     * testing that homogeneity of composition across chain lengths is
     * feasible after incorporating digest constraints;
     * Note: requires large memory (-Xmx8000M)
     * @param inDir input directory (ends with "\\")
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
     * computes chain length distribution. First a Gaussian (mu, sigma) restricted
     * between lmin and lmax is utilized as approximation. Then, projection on the
     * probability simplex with constraint to preserve mu is performed.
     * @param lmin smallest chain length
     * @param lmax longest chain length
     * @param sigma standard deviation of Gaussian approximation
     * @param mu desired average chain length
     * @return chain length distribution: between lmin and lmax with mean mu
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
     *
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        testRhoLHomFeas(inDir);
    }
}