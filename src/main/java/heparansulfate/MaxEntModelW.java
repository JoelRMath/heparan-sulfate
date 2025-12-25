package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * MaxEnt model of profiles of S/U composition and transition probabilities
 * along BKHS chains when BKHS is represented with a mixture of chain lengths
 */
public class MaxEntModelW {
    /**
     * solution (individual species abundances p) via geometric programming
     */
    MaxEntOptim opt = null;
    /**
     * number of disaccharides
     */
    int m = 2;
    /**
     * transition probabilities averaged over all positions, P[u][s]: from u to s
     */
    double[][] P = null;
    /**
     * disaccharides and overall proportions
     */
    BBSet bbs = null;
    /**
     * mixture model of chain lengths and sequence enumeration
     */
    MixSpecies msp = null;
    /**
     * prefix for output files
     */
    String outPref = null;

    /**
     * MaxEnt model of profiles of S/U composition and transition probabilities
     * along BKHS chains when BKHS is represented with a mixture of chain lengths;
     * saves MaxEnt individual species abundances (*.pind.res file)
     * @param outF prefix for output files
     * @param lmin smallest chain length
     * @param lmax largest chain length
     * @param sigma spread of chain lengths
     * @param mu average chain length
     * @param inDir input directory
     */
    public MaxEntModelW(String outF, int lmin, int lmax, double sigma, double mu, String inDir) {
        outPref = outF;
        LinEqCons lec = getWLEC(lmin, lmax, sigma, mu, inDir);
        opt = new MaxEntOptim(lec.A, lec.b);
        try {
            FileWriter fw = new FileWriter(outPref + ".pind.res");
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < msp.N; i++) {
                String s = String.valueOf(opt.p[i]);
                for (int j = 0; j < msp.seq[i].length; j++) {
                    s += "\t" + msp.seq[i][j];
                }
                bw.write(s);
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * default linear equality constraint with chain length distribution
     * = composition-1 + 2(digest-1) + chain length distribution -1
     * @param lmin smallest chain length
     * @param lmax largest chain length
     * @param sigma chain length spread
     * @param mu average chain length
     * @param inDir input directory
     * @return default linear equality constraint with chain length distribution
     */
    LinEqCons getWLEC(int lmin, int lmax, double sigma, double mu, String inDir) {
        msp = new MixSpecies(lmin, lmax, sigma, mu);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        List<LinEqCons> v = new ArrayList<>();
        v.add(LinEqCons.removeLastRow(msp.getRhoLEC(bbs)));
        v.add(LinEqCons.removeLastRow(msp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
                new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(msp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
                new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(msp.getWLEC()));
        return new LinEqCons(v);
    }

    /**
     * reads a file containing MaxEnt abundances of individual species
     * (pref.pind.res) and save files containing:
     * NRE and RE profiles of disaccharide composition (pref.prof.res)
     * disaccharide composition as a function of chain length (pref.rhol.res)
     * profile of transition probabilities when aligning by NRE (pref.ptnre.res)
     * profile of transition probabilities when aligning by RE (pref.ptre.res)
     * @param inDir input directory (ends with "/")
     * @param pref prefix for output files
     */
    public static void makeProf(String inDir, String pref) {
        double[] p = null;
        int[][] seq = null;
        int N = 0;
        try {
            FileReader fr = new FileReader(pref + ".pind.res");
            BufferedReader br = new BufferedReader(fr);
            while (br.readLine() != null) {
                N++;
            }
            br.close();
            fr.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        int lmax = 0;
        p = new double[N];
        seq = new int[N][];
        try {
            FileReader fr = new FileReader(pref + ".pind.res");
            BufferedReader br = new BufferedReader(fr);
            String line = null;
            int i = -1;
            while ((line = br.readLine()) != null) {
                i++;
                StringTokenizer st = new StringTokenizer(line, "\t");
                int n = st.countTokens() - 1;
                if (n > lmax) {
                    lmax = n;
                }
                seq[i] = new int[n];
                p[i] = Double.parseDouble(st.nextToken());
                int j = -1;
                while (st.hasMoreTokens()) {
                    j++;
                    seq[i][j] = Integer.parseInt(st.nextToken());
                }
            }
            br.close();
            fr.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        int lmin = lmax;
        int m = 2;
        String[] lab = new String[m];
        lab[0] = "U";
        lab[1] = "S";
        double[] lengthAtLeast = new double[lmax];
        double[][] nreProf = new double[lmax][m];
        for (int s = 0; s < N; s++) {
            for (int i = 0; i < seq[s].length; i++) {
                nreProf[i][seq[s][i]] += p[s];
                lengthAtLeast[i] += p[s];
            }
            if (seq[s].length < lmin) {
                lmin = seq[s].length;
            }
        }
        for (int i = 0; i < lmax; i++) {
            for (int j = 0; j < m; j++) {
                nreProf[i][j] /= lengthAtLeast[i];
            }
        }
        double[][] reProf = new double[lmax][m];
        for (int s = 0; s < N; s++) {
            for (int i = 0; i < seq[s].length; i++) {
                reProf[i][seq[s][seq[s].length - 1 - i]] += p[s];
            }
        }
        for (int i = 0; i < lmax; i++) {
            for (int j = 0; j < m; j++) {
                reProf[i][j] /= lengthAtLeast[i];
            }
        }
        List<String> v = new ArrayList<>();
        String sr = "nrePos";
        for (int i = 0; i < m; i++) {
            sr += "\tpnre" + lab[i];
        }
        sr += "\trePos";
        for (int i = 0; i < m; i++) {
            sr += "\tpre" + lab[i];
        }
        v.add(sr);
        for (int i = 0; i < lmax; i++) {
            sr = String.valueOf(i + 1);
            for (int j = 0; j < m; j++) {
                sr += "\t" + nreProf[i][j];
            }
            sr += "\t" + (lmax - i);
            for (int j = 0; j < m; j++) {
                sr += "\t" + reProf[i][j];
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".prof.res");
        int[] len = new int[lmax - lmin + 1];
        for (int i = 0; i < len.length; i++) {
            len[i] = lmin + i;
        }
        double[] w = new double[len.length];
        for (int s = 0; s < N; s++) {
            int l = seq[s].length - lmin;
            w[l] += p[s];
        }
        double[][] lenRho = new double[len.length][m];
        for (int s = 0; s < N; s++) {
            for (int i = 1; i < seq[s].length; i++) {
                int l = seq[s].length - lmin;
                lenRho[l][seq[s][i]] += p[s];
            }
        }
        for (int l = 0; l < len.length; l++) {
            for (int j = 0; j < m; j++) {
                lenRho[l][j] /= (double)(len[l] - 1);
                lenRho[l][j] /= w[l];
            }
        }
        v = new ArrayList<>();
        sr = "l\tw\trhoU\trhoS";
        v.add(sr);
        for (int l = 0; l < w.length; l++) {
            sr = String.valueOf(len[l]);
            sr += "\t" + w[l];
            sr += "\t" + lenRho[l][0];
            sr += "\t" + lenRho[l][1];
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".rhol.res");
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        double[][][] ptnre = new double[lmax][m][m];
        for (int s = 0; s < N; s++) {
            for (int i = 0; i < m; i++) {
                ptnre[0][i][seq[s][0]] += p[s];
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.;
            for (int j = 0; j < m; j++) {
                sum += ptnre[0][i][j];
            }
            for (int j = 0; j < m; j++) {
                ptnre[0][i][j] /= sum;
            }
        }
        for (int pos = 1; pos < lmax; pos++) {
            for (int s = 0; s < N; s++) {
                if (seq[s].length > pos) {
                    ptnre[pos][seq[s][pos - 1]][seq[s][pos]] += p[s];
                }
            }
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) {
                    sum += ptnre[pos][i][j];
                }
                for (int j = 0; j < m; j++) {
                    ptnre[pos][i][j] /= sum;
                }
            }
        }
        v = new ArrayList<>();
        sr = "pos";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                sr += "\tp" + bbs.name[i] + bbs.name[j];
            }
        }
        v.add(sr);
        for (int pos = 0; pos < lmax; pos++) {
            sr = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    sr += "\t" + ptnre[pos][i][j];
                }
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".ptnre.res");
        double[][][] ptre = new double[lmax][m][m];
        for (int s = 0; s < N; s++) {
            int pos = seq[s].length - lmax;
            if (pos >= 0) {
                for (int i = 0; i < m; i++) {
                    ptre[0][i][seq[s][pos]] += p[s];
                }
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.;
            for (int j = 0; j < m; j++) {
                sum += ptre[0][i][j];
            }
            for (int j = 0; j < m; j++) {
                ptre[0][i][j] /= sum;
            }
        }
        for (int pos = 1; pos < lmax; pos++) {
            for (int s = 0; s < N; s++) {
                int pos2 = seq[s].length - lmax + pos;
                if (pos2 >= 1) {
                    ptnre[pos][seq[s][pos2 - 1]][seq[s][pos2]] += p[s];
                }
            }
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) {
                    sum += ptnre[pos][i][j];
                }
                for (int j = 0; j < m; j++) {
                    ptnre[pos][i][j] /= sum;
                }
            }
        }
        v = new ArrayList<>();
        sr = "pos";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                sr += "\tp" + bbs.name[i] + bbs.name[j];
            }
        }
        v.add(sr);
        for (int pos = 0; pos < lmax; pos++) {
            sr = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    sr += "\t" + ptnre[pos][i][j];
                }
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".ptre.res");
    }

    /**
     * makes and saves individual species abundances for sigma = 3.5
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
     */
    public static void makeSig35(String inDir, String outDir) {
        new MaxEntModelW(outDir + "MEMW", 10, 20, 3.5, 16., inDir);
    }

    /**
     * makes and saves individual species abundances for sigma = 1.5
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
     */
    public static void makeSig15(String inDir, String outDir) {
        new MaxEntModelW(outDir + "MEMWSig1.5", 10, 20, 1.5, 16., inDir);
    }

    /**
     * saves distributions of chain lengths, composition as a function of chain length
     * and profiles (composition and transition probabilities)
     */
    public static void glAndRholAndProfExamples(String inDir, String outDir) {
        makeProf(inDir, outDir + "MEMW");
        makeProf(inDir, outDir + "MEMWSig1.5");
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        makeSig35(inDir, outDir);
        makeSig15(inDir, outDir);
        glAndRholAndProfExamples(inDir, outDir);
    }
}