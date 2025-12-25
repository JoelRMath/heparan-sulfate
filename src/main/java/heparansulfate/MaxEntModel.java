package heparansulfate;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * represents a MaxEnt model of BKHS (profiles of composition and correlation)
 * when all chains have same length
 */
public class MaxEntModel {
    /**
     * solution (species abundances p) via geometric programming
     */
    MaxEntOptim opt = null;
    /**
     * chain length
     */
    int n = 0;
    /**
     * number of disaccharides
     */
    int m = 0;
    /**
     * composition profile
     */
    double[][] gamma = null;
    /**
     * transition probabilities averaged over all positions, P[u][s]: from u to s
     */
    double[][] P = null;
    /**
     * transition probabilities, pt[i][u][s]: from u at position i to s at i+1
     */
    double[][][] pt = null;
    /**
     * disaccharides and overall proportions
     */
    BBSet bbs = null;
    /**
     * sequences of species
     */
    Species species = null;
    /**
     * numerical check
     */
    double sump = 0.;
    /**
     * prefix for output files
     */
    String outPref = null;

    /**
     * represents a MaxEnt model of BKHS (profiles of composition and correlation)
     * when all chains have same length
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
        try {
            FileWriter fw = new FileWriter(outPref + ".pind.res");
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < species.N; i++) {
                String s = String.valueOf(opt.p[i]);
                for (int j = 0; j < species.seq[i].length; j++) {
                    s += "\t" + species.seq[i][j];
                }
                bw.write(s);
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        makeGamma();
        makePT();
        makeP();
    }

    /**
     * estimates transition probabilities averaged over all positions
     */
    void makeP() {
        P = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = 0.;
                for (int pos = 1; pos < n; pos++) {
                    P[i][j] += pt[pos][i][j];
                }
                P[i][j] /= (double)(n - 1);
            }
        }
        List<String> v = new ArrayList<>();
        String s = "prev\tnext\tp";
        v.add(s);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                s = bbs.name[i] + "\t" + bbs.name[j] + "\t" + P[i][j];
                v.add(s);
            }
        }
        Utils.saveFile(v, outPref + ".P.res");
    }

    /**
     * estimates transition probabilities at each position
     */
    void makePT() {
        pt = new double[n][m][m];
        for (int s = 0; s < species.N; s++) {
            for (int i = 0; i < m; i++) {
                pt[0][i][species.seq[s][0]] += opt.p[s];
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.;
            for (int j = 0; j < m; j++) {
                sum += pt[0][i][j];
            }
            for (int j = 0; j < m; j++) {
                pt[0][i][j] /= sum;
            }
        }
        for (int pos = 1; pos < n; pos++) {
            for (int s = 0; s < species.N; s++) {
                pt[pos][species.seq[s][pos - 1]][species.seq[s][pos]] += opt.p[s];
            }
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) {
                    sum += pt[pos][i][j];
                }
                for (int j = 0; j < m; j++) {
                    pt[pos][i][j] /= sum;
                }
            }
        }
        List<String> v = new ArrayList<>();
        String s = "pos";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                s += "\tp" + bbs.name[i] + bbs.name[j];
            }
        }
        v.add(s);
        for (int pos = 0; pos < n; pos++) {
            s = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    s += "\t" + pt[pos][i][j];
                }
            }
            v.add(s);
        }
        Utils.saveFile(v, outPref + ".pt.res");
    }

    /**
     * estimates profile of S/U composition
     */
    void makeGamma() {
        gamma = new double[n][m];
        for (int s = 0; s < opt.n; s++) {
            for (int i = 0; i < n; i++) {
                gamma[i][species.seq[s][i]] += opt.p[s];
            }
        }
        List<String> v = new ArrayList<>();
        String s = "pos";
        for (int i = 0; i < m; i++) {
            s += "\t" + bbs.name[i];
        }
        v.add(s);
        for (int pos = 0; pos < n; pos++) {
            s = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                s += "\t" + gamma[pos][i];
            }
            v.add(s);
        }
        Utils.saveFile(v, outPref + ".gamma.res");
    }

    /**
     * saves all individual species abundances
     * @param file output file
     */
    void saveSpeciesAbundance(String file) {
        List<String> v = new ArrayList<>();
        v.add("species\tp");
        for (int i = 0; i < species.seq.length; i++) {
            String s = "";
            for (int j = 0; j < species.seq[i].length; j++) {
                s += bbs.name[species.seq[i][j]];
            }
            v.add(s + "\t" + opt.p[i]);
        }
        Utils.saveFile(v, file);
    }

    /**
     * default linear equality constraint = composition-1 + 2(digest-1)
     * @return default linear equality constraint = composition-1 + 2(digest-1)
     */
    public static LinEqCons getDefaultLEC(String inDir) {
        Species sp = new Species(2, 16);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        List<LinEqCons> v = new ArrayList<>();
        v.add(LinEqCons.removeLastRow(sp.getCompLEC(bbs.rho)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
                new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
                new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
        return new LinEqCons(v);
    }

    /**
     * MaxEnt models for chain length = 13, 14, 15, 16, 17 or 18
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
     */
    public static void defaultModelAndN(String inDir, String outDir) {
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        for (int n = 13; n <= 18; n++) {
            Species sp = new Species(2, n);
            List<LinEqCons> v = new ArrayList<>();
            v.add(LinEqCons.removeLastRow(sp.getCompLEC(bbs.rho)));
            v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
                    new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
            v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
                    new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
            LinEqCons lec = new LinEqCons(v);
            MaxEntOptim opt = new MaxEntOptim(lec.A, lec.b);
            new MaxEntModel(opt, n, bbs, outDir + "MEMN" + n);
        }
    }

    /**
     * saves maxent S/U composition at each position when all chains have same length n = 16
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
     */
    public static void defaultModel(String inDir, String outDir) {
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        LinEqCons lec = getDefaultLEC(inDir);
        MaxEntOptim opt = new MaxEntOptim(lec.A, lec.b);
        new MaxEntModel(opt, 16, bbs, outDir + "MEM");
    }

    /**
     * saves individual species abundances when all chains have same length n = 16
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
     */
    public static void defaultModelSpecies(String inDir, String outDir) {
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        LinEqCons lec = getDefaultLEC(inDir);
        MaxEntOptim opt = new MaxEntOptim(lec.A, lec.b);
        MaxEntModel mem = new MaxEntModel(opt, 16, bbs, outDir + "MEMspec");
        mem.saveSpeciesAbundance(outDir + "MEMspec.species.res");
    }

    /**
     * saves species abundances after removing the non-reducing end
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
     */
    public static void length15Abundances(String inDir, String outDir) {
        List<String> v = Utils.loadFileNoheader(outDir + "MEMspec.species.res");
        Map<String, Double> s2p = new HashMap<>();
        for (int i = 0; i < v.size(); i++) {
            StringTokenizer st = new StringTokenizer(v.get(i), "\t");
            String s = st.nextToken().trim();
            s = s.substring(1, s.length());
            double D = Double.parseDouble(st.nextToken());
            Double p = s2p.get(s);
            if (p == null) {
                p = 0.;
            }
            p += D;
            s2p.put(s, p);
        }
        SandVal[] sv = new SandVal[s2p.size()];
        int c = -1;
        for (String s : s2p.keySet()) {
            c++;
            Double p = s2p.get(s);
            sv[c] = new SandVal(s, p, "d");
        }
        Arrays.sort(sv, sv[0]);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        v = new ArrayList<>();
        v.add("species\tp\tp2");
        for (int i = 0; i < sv.length; i++) {
            double D = sv[i].val;
            String seq = sv[i].s;
            double p2 = getP2(seq, bbs);
            String seq2 = "x";
            for (int j = 0; j < seq.length(); j++) {
                String t = seq.substring(j, j + 1);
                if (t.equals("s")) {
                    seq2 += "{\\bf S}";
                } else {
                    seq2 += "U";
                }
            }
            v.add(seq2 + "\t" + D + "\t" + p2);
        }
        Utils.saveFile(v, outDir + "MEMspec.species.L15.res");
    }

    /**
     * relative abundance of a species (seq) under the MaxEnt model defined by only
     * overall disaccharide composition constraint (model H&I)
     * @param seq species sequence
     * @param bbs disaccharide overall abundances
     * @return relative abundance of a species (seq) under the MaxEnt model defined
     * by only overall disaccharide composition constraint (model H&I)
     */
    public static double getP2(String seq, BBSet bbs) {
        double res = 1.;
        for (int i = 0; i < seq.length(); i++) {
            String s = seq.substring(i, i + 1);
            Integer I = bbs.name2i.get(s);
            res *= bbs.rho[I];
        }
        return res;
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        defaultModel(inDir, outDir);
        defaultModelAndN(inDir, outDir);
        defaultModelSpecies(inDir, outDir);
        length15Abundances(inDir, outDir);
    }
}