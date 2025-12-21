package heparansulfate;

import java.util.Random;
import java.util.Vector;

/**
 * Homogeneous Markov model of species (chain sequence) abundances;
 * field this.g contains the expected fragment length distribution after digestion
 * by one heparinase.
 */
public class HCModel {
    int m = 0;
    int n = 0;
    double[][] P = null;
    double[][] pF = null;
    double[][][] pi = null;
    double c = 0.;
    BBSet bbs = null;
    CSpec cs = null;
    double[] g = null;
    double[] h = null;
    int lm = 0;
    double[] rhoF = null;

    public HCModel(int n, BBSet bbs, CSpec cs, double[][] tp, int lm) {
        this.n = n;
        this.bbs = bbs;
        this.m = bbs.m;
        this.cs = cs;
        this.lm = lm;
        P = new double[m][m];
        pF = new double[m][m];
        for (int i = 0; i < m; i++) {
            double s = 0.;
            for (int j = 0; j < m; j++) {
                s += tp[i][j];
            }
            for (int j = 0; j < m; j++) {
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

    void makeH() {
        h = new double[n];
        for (int i = 0; i < lm - 1; i++) {
            h[i] = g[i];
        }
        // computes expected distribution cumulative for length >= lm
        for (int i = lm - 1; i < n; i++) {
            h[lm - 1] += g[i];
        }
    }

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
        // length = 1
        g[0] = 0.;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                g[0] += bbs.rho[i] * cs.c[i] * cs.c[j] * P[i][j];
            }
        }
        g[0] /= c;
        g[0] *= (double) (n - 2);
        g[0] += 1.;
        g[0] /= (double) (n - 1);

        // 2 <= length <= n-1
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
                    tv *= (double) (n - ll - 1);
                    tv += 1.;
                    tv *= pi[l - 1][u][v];
                    tu += tv;
                }
                g[l] += bbs.rho[u] * cs.c[u] * tu;
            }
            g[l] /= (double) (n - 1);
            g[l] /= c;
        }
    }

    void makeC() {
        c = 0.;
        for (int i = 0; i < m; i++) {
            c += bbs.rho[i] * cs.c[i];
        }
    }

    int[] getSequence(Random rand) {
        int[] res = new int[n];
        res[0] = Utils.getRandIndexInF(rhoF, rand);
        for (int i = 1; i < n; i++) {
            res[i] = Utils.getRandIndexInF(pF[res[i - 1]], rand);
        }
        return res;
    }

    int[] getCuts(int[] seq, Random rand) {
        Vector<Integer> cuts = new Vector<Integer>();
        cuts.add(Integer.valueOf(0)); // Start cut at 0
        for (int i = 1; i < n; i++) {
            if (rand.nextDouble() < cs.c[seq[i]]) {
                cuts.add(Integer.valueOf(i));
            }
        }
        cuts.add(Integer.valueOf(n));
        int[] res = new int[cuts.size()];
        for (int i = 0; i < cuts.size(); i++) {
            res[i] = cuts.elementAt(i).intValue();
        }
        return res;
    }

    CleavedSequence getCutSequence(Random rand) {
        int[] seq = null;
        int[] cuts = new int[1];
        while (cuts.length < 2) {
            seq = getSequence(rand);
            cuts = getCuts(seq, rand);
        }
        return new CleavedSequence(seq, cuts);
    }

    void checkGL(Random rand) {
        int nsim = 1000000;
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

    public static double[] toVector(double[][] P) {
        double[] res = new double[P.length * P.length];
        int k = 0;
        for (int i = 0; i < P.length; i++) {
            for (int j = 0; j < P.length; j++) {
                res[k++] = P[i][j];
            }
        }
        return res;
    }

    public static double[][] toMatrix(int m, double[] p) {
        double[][] res = new double[m][m];
        int k = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                res[i][j] = p[k++];
            }
        }
        return res;
    }

    public static void main(String[] args) {
        String inDir = "input/";
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        CSpec cs = new CSpec(inDir + "US.hepI.txt", bbs);
        int m = cs.m;
        int n = 16;
        double[][] P = new double[m][m];
        Random rand = new Random();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = rand.nextDouble();
            }
        }
        // Note: HCModelLIC and ProjOnPolyHSet are required dependencies
        HCModelLIC mlic = new HCModelLIC(bbs);
        ProjOnPolyHSet proj = new ProjOnPolyHSet(toVector(P), mlic.A, mlic.b, rand);
        P = toMatrix(m, proj.optimum);
        int lm = 11;
        HCModel mm = new HCModel(n, bbs, cs, P, lm);
        mm.checkGL(rand);
    }
}