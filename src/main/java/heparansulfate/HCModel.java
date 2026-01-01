package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Homogeneous Markov model of species (chain sequence) abundances (model H&amp;C).
 * Field {@code this.g} contains the expected fragment length distribution after digestion
 * by one heparinase.
 */
public class HCModel {
    /**
     * Number of building blocks.
     */
    int m = 0;
    /**
     * Chain length.
     */
    int n = 0;
    /**
     * Matrix of transition probabilities.
     */
    double[][] P = null;
    /**
     * Cumulative P: {@code pF[i][j] = sum_{k=0}^j P[i][k]}, utilized for random drawing.
     */
    double[][] pF = null;
    /**
     * Powers of P.
     */
    double[][][] pi = null;
    /**
     * Overall cleavage probability.
     */
    double c = 0.;
    /**
     * Set of building blocks.
     */
    BBSet bbs = null;
    /**
     * Cleavage specificities/yields.
     */
    CSpec cs = null;
    /**
     * Fragment length distribution.
     */
    double[] g = null;
    /**
     * Fragment length distribution with cumulative abundance for {@code length >= lm}.
     */
    double[] h = null;
    /**
     * Maximum fragment length in experimental data.
     */
    int lm = 0;
    /**
     * Cumulative version of overall building-block composition in {@code this.bbs.rho}.
     */
    double[] rhoF = null;

    /**
     * Homogeneous Markov model of species (chain sequence) abundances (model H&amp;C).
     * Field {@code this.g} contains the expected fragment length distribution after digestion
     * by one heparinase.
     * @param n BKHS chain length.
     * @param bbs Disaccharides and their overall proportions.
     * @param cs Cleavage specificities/yields for one heparinase.
     * @param tp Initial matrix of transition probabilities (random or already optimized)
     * which is then projected to satisfy constraints: output in {@code this.P}.
     * @param lm Maximum fragment length.
     */
    public HCModel(int n, BBSet bbs, CSpec cs, double[][] tp, int lm) {
        this.n = n;
        this.bbs = bbs;
        m = bbs.m;
        this.cs = cs;
        this.lm = lm;
        P = new double[m][m];
        pF = new double[m][m];
        for (int i = 0; i < m; i++) {
            double s = 0.;
            for (int j = 0; j < m; j++) {
                s += tp[i][j];
            }
            for (int j = 0; j < m; j++) { // projection is not always perfect due to numerical accuracy
                P[i][j] = tp[i][j] / s;
            }
            pF[i][0] = P[i][0];
            for (int j = 1; j < m; j++) {
                pF[i][j] = pF[i][j - 1] + P[i][j];
            }
        }
        rhoF = new double[m];
        rhoF[0] = bbs.rho[0];
        for (int i = 1; i < m; i++) {
            rhoF[i] = rhoF[i - 1] + bbs.rho[i];
        }
        makeC();
        makeG();
        makeH();
    }

    /**
     * Computes expected distribution of fragment length ({@code this.h}) cumulative for {@code length >= lm}.
     */
    void makeH() {
        h = new double[n];
        for (int i = 0; i < lm; i++) {
            h[i] = g[i];
        }
        for (int i = lm; i < n - 1; i++) {
            h[lm - 1] += g[i];
        }
    }

    /**
     * Computes expected distribution of fragment length ({@code this.g}).
     */
    void makeG() {
        pi = new double[n][m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                pi[0][i][j] = P[i][j] * (1. - cs.c[j]);
            }
        }
        for (int i = 1; i < n; i++) {
            pi[i] = MatrixOp.multMat(pi[i - 1], pi[0]);
        }
        g = new double[n];
        // ll = 1
        g[0] = 0.;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                g[0] += bbs.rho[i] * cs.c[i] * cs.c[j] * P[i][j];
            }
        }
        g[0] /= c;
        g[0] *= (double)(n - 2);
        g[0] += 1.;
        g[0] /= (double)(n - 1);
        // 2 <= ll <= n-1
        for (int l = 1; l < n - 1; l++) {
            int ll = l + 1;
            g[l] = 0.;
            for (int u = 0; u < m; u++) {
                double tu = 0.;
                for (int v = 0; v < m; v++) {
                    double tv = 0.;
                    for (int w = 0; w < m; w++) {
                        tv += P[v][w] * cs.c[w];
                    }
                    tv *= (double)(n - ll - 1);
                    tv += 1.;
                    tv *= pi[l - 1][u][v];
                    tu += tv;
                }
                g[l] += bbs.rho[u] * cs.c[u] * tu;
            }
            g[l] /= (double)(n - 1);
            g[l] /= c;
        }
    }

    /**
     * Computes parameter {@code c} of H&amp;C.
     */
    void makeC() {
        c = 0.;
        for (int i = 0; i < m; i++) {
            c += bbs.rho[i] * cs.c[i];
        }
    }

    /**
     * Generates a random BKHS sequence based on {@code this.P}.
     * @param rand Random number generator.
     * @return A random BKHS sequence based on {@code this.P}.
     */
    int[] getSequence(Random rand) {
        int[] res = new int[n];
        res[0] = Utils.getRandIndexInF(rhoF, rand);
        for (int i = 1; i < n; i++) {
            res[i] = Utils.getRandIndexInF(pF[res[i - 1]], rand);
        }
        return res;
    }

    /**
     * Generates a random set of cleavage positions in BKHS sequence {@code seq} based
     * on heparinase cleavage specificities.
     * @param seq BKHS sequence.
     * @param rand Random number generator.
     * @return A random set of cleavage positions in BKHS sequence {@code seq} based
     * on heparinase cleavage specificities.
     */
    int[] getCuts(int[] seq, Random rand) {
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
     * Generates a random {@code CleavedSequence} (generates a BKHS chain and randomly cleaves
     * it based on cleavage specificities); see method {@code getFragments()} of class
     * {@code CleavedSequence} to access resulting fragments.
     * @param rand Random number generator.
     * @return A random {@code CleavedSequence} (generates a BKHS chain and randomly cleaves
     * it based on cleavage specificities).
     */
    CleavedSequence getCutSequence(Random rand) {
        int[] seq = null;
        int[] cuts = new int[1];
        while (cuts.length < 2) {
            seq = getSequence(rand);
            cuts = getCuts(seq, rand);
        }
        return new CleavedSequence(seq, cuts);
    }

    /**
     * Numerical check for {@code this.g}.
     * @param rand Random number generator.
     */
    void checkGL(Random rand) {
        int nsim = 10000000;
        int totfrag = 0;
        double[] f = new double[n];
        for (int sim = 0; sim < nsim; sim++) {
            CleavedSequence cseq = getCutSequence(rand);
            int[][] frag = cseq.getFragments();
            totfrag += frag.length;
            for (int i = 0; i < frag.length; i++) {
                f[frag[i].length - 1] += 1.;
            }
        }
        for (int i = 0; i < n; i++) {
            System.out.println((i + 1) + "\t" + g[i] + "\t" + (f[i] / (double) totfrag));
        }
    }

    /**
     * Returns a vector version of matrix {@code P}.
     * @param P Matrix of transition probabilities.
     * @return A vector version of matrix {@code P}.
     */
    public static double[] toVector(double[][] P) {
        double[] res = new double[P.length * P.length];
        int k = -1;
        for (int i = 0; i < P.length; i++) {
            for (int j = 0; j < P.length; j++) {
                k++;
                res[k] = P[i][j];
            }
        }
        return res;
    }

    /**
     * Returns a matrix version of vector {@code p}.
     * @param m Number of rows/columns.
     * @param p Vector.
     * @return A {@code m} by {@code m} matrix version of vector {@code p}.
     */
    public static double[][] toMatrix(int m, double[] p) {
        double[][] res = new double[m][m];
        int k = -1;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                k++;
                res[i][j] = p[k];
            }
        }
        return res;
    }

    /**
     * For testing.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        file = inDir + "US.hepI.txt";
        CSpec cs = new CSpec(file, bbs);
        int m = cs.m;
        int n = 16;
        double[][] P = new double[m][m];
        Random rand = new Random();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = rand.nextDouble();
            }
        }
        HCModelLIC mlic = new HCModelLIC(bbs);
        ProjOnPolyHSet proj = new ProjOnPolyHSet(toVector(P), mlic.A, mlic.b, rand);
        P = toMatrix(m, proj.optimum);
        int lm = 11;
        HCModel mm = new HCModel(n, bbs, cs, P, lm);
        mm.checkGL(rand);
    }
}