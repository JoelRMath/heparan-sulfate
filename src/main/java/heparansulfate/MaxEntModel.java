package heparansulfate;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Represents a MaxEnt model of BKHS (profiles of composition and correlation)
 * when all chains have same length.
 */
public class MaxEntModel {
    MaxEntOptim opt = null;
    int n = 0;
    int m = 0;
    double[][] gamma = null;
    double[][] P = null;
    double[][][] pt = null;
    BBSet bbs = null;
    Species species = null;
    double sump = 0.;
    String outPref = null;

    /**
     * @param op MaxEnt model of individual species abundances
     * @param n BKHS chain length
     * @param bbs disaccharides
     * @param outF prefix for output files
     */
    public MaxEntModel(MaxEntOptim op, int n, BBSet bbs, String outF) {
        outPref = outF;
        opt = op;
        this.n = n;
        this.bbs = bbs;
        m = bbs.m;
        species = new Species(m, n);
        sump = 0.;
        for (int s = 0; s < opt.n; s++) {
            sump += opt.p[s];
        }
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outPref + ".pind.res"))) {
            for (int i = 0; i < species.N; i++) {
                String s = String.valueOf(opt.p[i]);
                for (int j = 0; j < species.seq[i].length; j++) {
                    s += "\t" + species.seq[i][j];
                }
                bw.write(s);
                bw.newLine();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        makeGamma();
        makePT();
        makeP();
    }

    void makeP() {
        P = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = 0.;
                for (int pos = 1; pos < n; pos++) {
                    P[i][j] += pt[pos][i][j];
                }
                P[i][j] /= (double) (n - 1);
            }
        }
        Vector<String> v = new Vector<>();
        v.add("prev\tnext\tp");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                v.add(bbs.name[i] + "\t" + bbs.name[j] + "\t" + P[i][j]);
            }
        }
        Utils.saveFile(v, outPref + ".P.res");
    }

    void makePT() {
        pt = new double[n][m][m];
        for (int s = 0; s < species.N; s++) {
            for (int i = 0; i < m; i++) {
                pt[0][i][species.seq[s][0]] += opt.p[s];
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.;
            for (int j = 0; j < m; j++) sum += pt[0][i][j];
            for (int j = 0; j < m; j++) if (sum > 0) pt[0][i][j] /= sum;
        }
        for (int pos = 1; pos < n; pos++) {
            for (int s = 0; s < species.N; s++) {
                pt[pos][species.seq[s][pos - 1]][species.seq[s][pos]] += opt.p[s];
            }
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) sum += pt[pos][i][j];
                for (int j = 0; j < m; j++) if (sum > 0) pt[pos][i][j] /= sum;
            }
        }
        Vector<String> v = new Vector<>();
        StringBuilder header = new StringBuilder("pos");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) header.append("\tp").append(bbs.name[i]).append(bbs.name[j]);
        }
        v.add(header.toString());
        for (int pos = 0; pos < n; pos++) {
            StringBuilder s = new StringBuilder(String.valueOf(pos + 1));
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) s.append("\t").append(pt[pos][i][j]);
            }
            v.add(s.toString());
        }
        Utils.saveFile(v, outPref + ".pt.res");
    }

    void makeGamma() {
        gamma = new double[n][m];
        for (int s = 0; s < opt.n; s++) {
            for (int i = 0; i < n; i++) gamma[i][species.seq[s][i]] += opt.p[s];
        }
        Vector<String> v = new Vector<>();
        StringBuilder header = new StringBuilder("pos");
        for (int i = 0; i < m; i++) header.append("\t").append(bbs.name[i]);
        v.add(header.toString());
        for (int pos = 0; pos < n; pos++) {
            StringBuilder s = new StringBuilder(String.valueOf(pos + 1));
            for (int i = 0; i < m; i++) s.append("\t").append(gamma[pos][i]);
            v.add(s.toString());
        }
        Utils.saveFile(v, outPref + ".gamma.res");
    }

    void saveSpeciesAbundance(String file) {
        Vector<String> v = new Vector<>();
        v.add("species\tp");
        for (int i = 0; i < species.seq.length; i++) {
            StringBuilder s = new StringBuilder();
            for (int j = 0; j < species.seq[i].length; j++) s.append(bbs.name[species.seq[i][j]]);
            v.add(s.toString() + "\t" + opt.p[i]);
        }
        Utils.saveFile(v, file);
    }

    public static LinEqCons getDefaultLEC(String inDir) {
        Species sp = new Species(2, 16);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        Vector<LinEqCons> v = new Vector<>();
        v.add(LinEqCons.removeLastRow(sp.getCompLEC(bbs.rho)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
                new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
                new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
        return new LinEqCons(v);
    }

    public static double getP2(String seq, BBSet bbs) {
        double res = 1.;
        for (int i = 0; i < seq.length(); i++) {
            String s = seq.substring(i, i + 1);
            Integer index = bbs.name2i.get(s);
            if (index != null) res *= bbs.rho[index];
        }
        return res;
    }

    public static void length15Abundances(String inDir, String outDir) {
        Vector<String> v = Utils.loadFileNoheader(outDir + "MEMspec.species.res");
        Hashtable<String, Double> s2p = new Hashtable<>();
        for (int i = 0; i < v.size(); i++) {
            StringTokenizer st = new StringTokenizer(v.elementAt(i), "\t");
            String s = st.nextToken().trim();
            s = s.substring(1);
            double D = Double.parseDouble(st.nextToken());
            Double p = s2p.get(s);
            s2p.put(s, (p == null ? 0.0 : p) + D);
        }
        SandVal[] sv = new SandVal[s2p.size()];
        Enumeration<String> en = s2p.keys();
        int c = 0;
        while (en.hasMoreElements()) {
            String s = en.nextElement();
            sv[c++] = new SandVal(s, s2p.get(s), "d");
        }
        Arrays.sort(sv, sv[0]);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        Vector<String> outV = new Vector<>();
        outV.add("species\tp\tp2");
        for (SandVal sval : sv) {
            String seq = sval.s;
            double p2 = getP2(seq, bbs);
            StringBuilder seq2 = new StringBuilder("x");
            for (int j = 0; j < seq.length(); j++) {
                seq2.append(seq.charAt(j) == 's' ? "{\\bf S}" : "U");
            }
            outV.add(seq2.toString() + "\t" + sval.val + "\t" + p2);
        }
        Utils.saveFile(outV, outDir + "MEMspec.species.L15.res");
    }

    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        // Execution methods here
    }
}
