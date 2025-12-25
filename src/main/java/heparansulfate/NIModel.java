package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * represents a nonhomogeneity and independence (N&I) model:
 * different building block composition at each position in a chain
 * and independence between positions; field this.g contains the expected
 * heparinase fragment length distribution under model N&I
 */
public class NIModel {
    /**
     * length of a chain
     */
    int n = 0;
    /**
     * number of building block types
     */
    int m = 0;
    /**
     * gamma[i][j] is the proportion of building block j at position i from
     * the non-reducing end
     */
    double[][] gamma = null;
    /**
     * cumulative version of gamma: gamF[i][j] = sum_{k=0}^j gamma[i][k]
     */
    double[][] gamF = null;
    /**
     * overall cleavage probability
     */
    double c = 0;
    /**
     * see equation of g(l) in paper
     */
    double alpha = 0.;
    /**
     * cleavage probability at position i
     */
    double[] ci = null;
    /**
     * set of disaccharides and their overall proportions
     */
    BBSet bbs = null;
    /**
     * cleavage specificities/yields for one heparinase
     */
    CSpec cs = null;
    /**
     * expected distribution of fragment length in on heparinase digest
     */
    double[] g = null;
    /**
     * same as this.g but cumulative for length ll >= lm
     */
    double[] h = null;
    /**
     * maximum measured fragment length
     */
    int lm = 0;

    /**
     * represents a nonhomogeneity and independence (N&I) model:
     * different building block composition at each position in a chain
     * and independence between positions
     * @param bbs disaccharides and their overall proportions
     * @param cs cleavage specificities/yields for one heparinase
     * @param gam initial (random or optimized) matrix Gamma of disaccharide
     * proportions along chains;
     * this initial matrix is then projected to satisfy constraints
     * @param lm maximum measured fragment length
     */
    public NIModel(BBSet bbs, CSpec cs, double[][] gam, int lm) {
        this.lm = lm;
        this.bbs = bbs;
        this.cs = cs;
        m = bbs.m;
        n = gam.length;
        gamma = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamma[i][j] = gam[i][j];
                // projection is not always perfect due to numerical accuracy
                if (gamma[i][j] > 1.) {
                    gamma[i][j] = 1.;
                }
                if (gamma[i][j] < 0.) {
                    gamma[i][j] = 0.;
                }
            }
        }
        gamF = new double[n][m];
        for (int i = 0; i < n; i++) {
            gamF[i][0] = gamma[i][0];
            for (int j = 1; j < m; j++) {
                gamF[i][j] = gamF[i][j - 1] + gamma[i][j];
            }
        }
        makeC();
        makeCi();
        makeG();
        makeH();
    }

    /**
     * this.h: fragment length distribution expected under N&I
     * and cumulative for length ll >= lm
     */
    void makeH() {
        h = new double[lm];
        int max = lm;
        if (max > n - 1) {
            max = n - 1;
        }
        for (int i = 0; i < max; i++) {
            h[i] = g[i];
        }
        for (int i = lm; i < n; i++) {
            h[lm - 1] += g[i];
        }
    }

    /**
     * generates one random BKHS sequence based on this.gamma
     * @param rand
     * @return one random BKHS sequence based on this.gamma
     */
    int[] getSequence(Random rand) {
        int[] res = new int[n];
        for (int i = 0; i < n; i++) {
            res[i] = Utils.getRandIndexInF(gamF[i], rand);
        }
        return res;
    }

    /**
     * generate a random set of cleavage positions in seq based on
     * heparinase cleavage specificities/yields
     * @param seq BKHS sequence
     * @param rand
     * @return a random set of cleavage positions in seq based on
     * heparinase cleavage specificities/yields
     */
    int[] getCleavages(int[] seq, Random rand) {
        List<Integer> cuts = new ArrayList<>();
        for (int i = 1; i < n; i++) {
            double d = rand.nextDouble();
            if (d < cs.c[seq[i]]) {
                cuts.add(i);
            }
        }
        cuts.add(n);
        int[] res = new int[cuts.size()];
        for (int i = 0; i < cuts.size(); i++) {
            res[i] = cuts.get(i);
        }
        return res;
    }

    /**
     * generates one BKHS sequence and cleaves;
     * see getFragments() of CleavedSequence to obtain fragments
     * @param rand
     * @return one cleaved BKHS chain
     */
    CleavedSequence getCleavageSequence(Random rand) {
        int[] seq = null;
        int[] cuts = new int[1];
        while (cuts.length < 2) {
            seq = getSequence(rand);
            cuts = getCleavages(seq, rand);
        }
        return new CleavedSequence(seq, cuts);
    }

    /**
     * distribution of fragment length after heparinase digestion
     * and expected under N&I
     */
    void makeG() {
        g = new double[n];
        double[][] omega = new double[n][n];
        for (int i = 0; i < n; i++) {
            omega[i][i] = 1. - ci[i];
            for (int j = i + 1; j < n; j++) {
                omega[i][j] = omega[i][j - 1] * (1. - ci[j]);
            }
        }
        alpha = 0.;
        for (int i = 1; i < n; i++) {
            alpha += ci[i];
        }
        // ll = 1
        g[0] = ci[n - 1];
        for (int i = 1; i < n - 1; i++) {
            g[0] += ci[i] * ci[i + 1];
        }
        g[0] /= alpha;
        // 2 <= ll <= n-2
        for (int ll = 2; ll <= n - 2; ll++) {
            int l = ll - 1;
            g[l] = ci[n - ll] * omega[n - ll + 1][n - 1];
            for (int i = 1; i < n - ll; i++) {
                g[l] += ci[i] * ci[i + ll] * omega[i + 1][i + ll - 1];
            }
            g[l] /= alpha;
        }
        // ll = n-1
        g[n - 2] = ci[1] * omega[2][n - 1];
        g[n - 2] /= alpha;
    }

    /**
     * this.ci: cleavage probabilities at each position i
     * (average quantity in the mixture, defined by this.gamma and this.cs)
     */
    void makeCi() {
        ci = new double[n];
        for (int i = 1; i < n; i++) {
            ci[i] = 0.;
            for (int j = 0; j < m; j++) {
                ci[i] += gamma[i][j] * cs.c[j];
            }
        }
    }

    /**
     * overall cleavage probability c
     */
    void makeC() {
        c = 0.;
        for (int i = 0; i < m; i++) {
            c += bbs.rho[i] * cs.c[i];
        }
    }

    /**
     * numerical check for this.g
     * @param rand
     */
    void checkGL(Random rand) {
        int nsim = 10000000;
        int totfrag = 0;
        double[] f = new double[n];
        for (int sim = 0; sim < nsim; sim++) {
            CleavedSequence cseq = getCleavageSequence(rand);
            int[][] frag = cseq.getFragments();
            totfrag += frag.length;
            for (int i = 0; i < frag.length; i++) {
                f[frag[i].length - 1] += 1.;
            }
        }
        for (int i = 0; i < n; i++) {
            double d = 0.;
            if (i < lm) {
                d = h[i];
            }
            System.out.println((i + 1) + "\t" + g[i] + "\t" +
                    (f[i] / (double) totfrag) + "\t" + f[i] + "\t" + d);
        }
    }

    /**
     * turns a gamma matrix (n by m) in vector form
     * @param G gamma matrix (n by m)
     * @return a vector representing matrix G
     */
    public static double[] toVector(double[][] G) {
        double[] res = new double[G.length * G[0].length];
        int k = -1;
        for (int j = 0; j < G[0].length; j++) {
            for (int i = 0; i < G.length; i++) {
                k++;
                res[k] = G[i][j];
            }
        }
        return res;
    }

    /**
     * turns a vector into a gamma matrix
     * @param g vector of dimension (m*nn)
     * @param nn BKHS chain length
     * @return a matrix representing a gamma vector
     */
    public static double[][] toMatrix(double[] g, int nn) {
        int mm = g.length / nn;
        double[][] res = new double[nn][mm];
        int k = -1;
        for (int i = 0; i < mm; i++) {
            for (int j = 0; j < nn; j++) {
                k++;
                res[j][i] = g[k];
            }
        }
        return res;
    }
}