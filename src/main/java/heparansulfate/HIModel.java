package heparansulfate;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Homogeneity and Independence (H&amp;I) model: species abundances are defined by
 * the overall disaccharide composition.
 */
public class HIModel {
    /**
     * Number of building blocks.
     */
    int m = 0;
    /**
     * Chain length.
     */
    int n = 0;
    /**
     * Overall cleavage probability.
     */
    double c = 0;
    /**
     * Set of building blocks.
     */
    BBSet bbs = null;
    /**
     * Cleavage specificities.
     */
    CSpec cs = null;
    /**
     * Distribution of fragment length.
     */
    double[] g = null;
    /**
     * Distribution of fragment length with upper bound {@code lm}.
     */
    double[] h = null;
    /**
     * Upper bound for fragment length.
     */
    int lm = 0;
    /**
     * Cumulative version of overall building-block proportions, utilized for MC simulation.
     */
    double[] rhoF = null;

    /**
     * Computes heparinase fragment length distribution under H&amp;I ({@code this.g}).
     * @param n BKHS chain length.
     * @param bbs Set of building blocks (disaccharides and their overall proportions).
     * @param cs Cleavage specificity/yield (one heparinase only).
     * @param lm Maximum fragment length (cumulative for {@code h(l)} when {@code l+1 >= lm}).
     */
    public HIModel(int n, BBSet bbs, CSpec cs, int lm) {
        this.lm = lm;
        this.bbs = bbs;
        this.cs = cs;
        m = bbs.m;
        rhoF = new double[m];
        rhoF[0] = bbs.rho[0];
        for (int i = 1; i < m; i++) {
            rhoF[i] = rhoF[i - 1] + bbs.rho[i];
        }
        this.n = n;
        makeC();
        makeG();
        makeH();
    }

    /**
     * Numerical check via simulations.
     * @param rand Random number generator.
     * @param file Output file.
     */
    void checkGL(Random rand, String file) {
        int nsim = 10000000;
        int totfrag = 0;
        double[] g = new double[n];
        for (int sim = 0; sim < nsim; sim++) {
            CleavedSequence cseq = getCleavedSequence(rand);
            int[][] frag = cseq.getFragments();
            totfrag += frag.length;
            for (int i = 0; i < frag.length; i++) {
                g[frag[i].length - 1] += 1.;
            }
        }
        for (int i = 0; i < g.length; i++) {
            g[i] /= (double) totfrag;
        }
        double[] h = new double[n];
        for (int i = 0; i < lm; i++) {
            h[i] = g[i];
        }
        for (int i = lm; i < n; i++) {
            h[lm - 1] += g[i];
        }
        List<String> v = new ArrayList<>();
        v.add("l\tgl\tsimgl\thl\tsimhl");
        for (int i = 0; i < n; i++) {
            int L = i + 1;
            String outS = String.valueOf(L);
            outS += "\t" + this.g[i];
            outS += "\t" + g[i];
            outS += "\t" + this.h[i];
            outS += "\t" + h[i];
            v.add(outS);
        }
        Utils.saveFile(v, file);
    }

    /**
     * Returns a random chain sequence before cleavage.
     * @param rand Random number generator.
     * @return A random chain sequence before cleavage.
     */
    int[] getSequence(Random rand) {
        int[] res = new int[n];
        for (int i = 0; i < n; i++) {
            res[i] = Utils.getRandIndexInF(rhoF, rand);
        }
        return res;
    }

    /**
     * Generates a random set of cleavage positions in sequence {@code seq} based
     * on cleavage specificities (always a cut at position {@code n}).
     * @param seq BKHS sequence.
     * @param rand Random number generator.
     * @return A random set of cleavage positions in sequence {@code seq} based on
     * cleavage specificities.
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
     * Generates a random {@code CleavedSequence} (generates a BKHS chain and randomly
     * cleaves it based on cleavage specificities). See methods {@code getFragments()}
     * of class {@code CleavedSequence} to access resulting fragments.
     * @param rand Random number generator.
     * @return A random {@code CleavedSequence} (generates a BKHS chain and randomly cleaves
     * it based on cleavage specificities).
     */
    CleavedSequence getCleavedSequence(Random rand) {
        int[] seq = null;
        int[] cuts = new int[1];
        while (cuts.length < 2) { // only fragments obtained by cleavage are visible
            seq = getSequence(rand);
            cuts = getCuts(seq, rand);
        }
        return new CleavedSequence(seq, cuts);
    }

    /**
     * Computes version {@code this.h} of {@code this.g} (fragment length distribution).
     * Abundances for {@code ll >= lm} are summed.
     */
    void makeH() {
        h = new double[n];
        for (int l = 0; l < lm; l++) {
            h[l] = g[l];
        }
        for (int l = lm; l < n - 1; l++) {
            h[lm - 1] += g[l];
        }
    }

    /**
     * Computes the distribution of fragment length ({@code this.g}) expected under H&amp;I.
     */
    void makeG() {
        g = new double[n];
        for (int ll = 1; ll <= n - 1; ll++) {
            int l = ll - 1;
            g[l] = 1. + c * (double)(n - ll - 1);
            g[l] /= (double)(n - 1);
            g[l] *= Math.pow((1. - c), (double) l);
        }
    }

    /**
     * Computes parameter {@code this.c} used to compute {@code this.g}.
     */
    void makeC() {
        c = 0.;
        for (int i = 0; i < m; i++) {
            c += bbs.rho[i] * cs.c[i];
        }
    }

    /**
     * Saves {@code this.h} to a file.
     * @param file Output file.
     */
    void saveH(String file) {
        try {
            FileWriter fw = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("l\th");
            bw.newLine();
            for (int i = 0; i < lm; i++) {
                bw.write((i + 1) + "\t" + h[i]);
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Saves fragment length distribution expected under H&amp;I for HepI
     * and performs numerical check.
     * @param inDir Input directory (ends with "/").
     * @param outDir Output directory (ends with "/").
     */
    public static void hepI(String inDir, String outDir) {
        Random rand = new Random();
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        file = inDir + "US.hepI.txt";
        CSpec cs = new CSpec(file, bbs);
        int n = 16;
        HIModel pm = new HIModel(n, bbs, cs, 11);
        pm.checkGL(rand, outDir + "HIM.hepIgl.check.res");
        pm.saveH(outDir + "HIM.hepIhl.res");
    }

    /**
     * Saves fragment length distribution expected under H&amp;I for HepIII
     * and performs numerical check.
     * @param inDir Input directory (ends with "/").
     * @param outDir Output directory (ends with "/").
     */
    public static void hepIII(String inDir, String outDir) {
        Random rand = new Random();
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        file = inDir + "US.hepIII.txt";
        CSpec cs = new CSpec(file, bbs);
        int n = 16;
        HIModel pm = new HIModel(n, bbs, cs, 6);
        pm.checkGL(rand, outDir + "HIM.hepIIIgl.check.res");
        pm.saveH(outDir + "HIM.hepIIIhl.res");
    }

    /**
     * Main entry point.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        hepI(inDir, outDir);
        hepIII(inDir, outDir);
    }
}