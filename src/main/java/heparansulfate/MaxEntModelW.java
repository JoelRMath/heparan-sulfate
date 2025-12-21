package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * MaxEnt model of profiles of S/U composition and transition probabilities
 * along BKHS chains when BKHS is represented with a mixture of chain lengths.
 */
public class MaxEntModelW {
    MaxEntOptim opt = null;
    int m = 2;
    double[][] P = null;
    BBSet bbs = null;
    MixSpecies msp = null;
    String outPref = null;

    public MaxEntModelW(String outF, int lmin, int lmax, double sigma, double mu, String inDir) {
        outPref = outF;
        LinEqCons lec = getWLEC(lmin, lmax, sigma, mu, inDir);
        opt = new MaxEntOptim(lec.A, lec.b);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outPref + ".pind.res"))) {
            for (int i = 0; i < msp.N; i++) {
                StringBuilder s = new StringBuilder(String.valueOf(opt.p[i]));
                for (int j = 0; j < msp.seq[i].length; j++) {
                    s.append("\t").append(msp.seq[i][j]);
                }
                bw.write(s.toString());
                bw.newLine();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    LinEqCons getWLEC(int lmin, int lmax, double sigma, double mu, String inDir) {
        msp = new MixSpecies(lmin, lmax, sigma, mu);
        bbs = new BBSet(inDir + "US.ab.txt");
        Vector<LinEqCons> v = new Vector<>();
        v.add(LinEqCons.removeLastRow(msp.getRhoLEC(bbs)));
        v.add(LinEqCons.removeLastRow(msp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
                new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(msp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
                new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(msp.getWLEC()));
        return new LinEqCons(v);
    }

    public static void makeProf(String inDir, String pref) {
        double[] p = null;
        int[][] seq = null;
        int N = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(pref + ".pind.res"))) {
            while (br.readLine() != null) N++;
        } catch (Exception ex) { ex.printStackTrace(); }

        int lmax = 0;
        p = new double[N];
        seq = new int[N][];
        try (BufferedReader br = new BufferedReader(new FileReader(pref + ".pind.res"))) {
            String line;
            int i = 0;
            while ((line = br.readLine()) != null) {
                StringTokenizer st = new StringTokenizer(line, "\t");
                int n = st.countTokens() - 1;
                if (n > lmax) lmax = n;
                seq[i] = new int[n];
                p[i] = Double.parseDouble(st.nextToken());
                int j = 0;
                while (st.hasMoreTokens()) {
                    seq[i][j++] = Integer.parseInt(st.nextToken());
                }
                i++;
            }
        } catch (Exception ex) { ex.printStackTrace(); }

        int lmin = lmax;
        int m = 2;
        String[] lab = {"U", "S"};
        double[] lengthAtLeast = new double[lmax];
        double[][] nreProf = new double[lmax][m];
        
        for (int s = 0; s < N; s++) {
            for (int i = 0; i < seq[s].length; i++) {
                nreProf[i][seq[s][i]] += p[s];
                lengthAtLeast[i] += p[s];
            }
            if (seq[s].length < lmin) lmin = seq[s].length;
        }

        for (int i = 0; i < lmax; i++) {
            for (int j = 0; j < m; j++) if (lengthAtLeast[i] > 0) nreProf[i][j] /= lengthAtLeast[i];
        }

        double[][] reProf = new double[lmax][m];
        for (int s = 0; s < N; s++) {
            for (int i = 0; i < seq[s].length; i++) {
                reProf[i][seq[s][seq[s].length - 1 - i]] += p[s];
            }
        }
        for (int i = 0; i < lmax; i++) {
            for (int j = 0; j < m; j++) if (lengthAtLeast[i] > 0) reProf[i][j] /= lengthAtLeast[i];
        }

        Vector<String> v = new Vector<>();
        StringBuilder header = new StringBuilder("nrePos");
        for (String l : lab) header.append("\tpnre").append(l);
        header.append("\trePos");
        for (String l : lab) header.append("\tpre").append(l);
        v.add(header.toString());

        for (int i = 0; i < lmax; i++) {
            StringBuilder sr = new StringBuilder(String.valueOf(i + 1));
            for (int j = 0; j < m; j++) sr.append("\t").append(nreProf[i][j]);
            sr.append("\t").append(lmax - i);
            for (int j = 0; j < m; j++) sr.append("\t").append(reProf[i][j]);
            v.add(sr.toString());
        }
        Utils.saveFile(v, pref + ".prof.res");

        // Chain Length Distribution logic
        double[] w = new double[lmax - lmin + 1];
        for (int s = 0; s < N; s++) w[seq[s].length - lmin] += p[s];

        v = new Vector<>();
        v.add("l\tw");
        for (int i = 0; i < w.length; i++) v.add((lmin + i) + "\t" + w[i]);
        Utils.saveFile(v, pref + ".rhol.res");

        // Transition probabilities NRE
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        double[][][] ptnre = new double[lmax][m][m];
        for (int s = 0; s < N; s++) {
            for (int pos = 1; pos < seq[s].length; pos++) {
                ptnre[pos][seq[s][pos - 1]][seq[s][pos]] += p[s];
            }
        }
        // ... (Normalization and Saving ptnre)
    }

    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        new MaxEntModelW(outDir + "MEMW", 10, 20, 3.5, 16., inDir);
        makeProf(inDir, outDir + "MEMW");
    }
}