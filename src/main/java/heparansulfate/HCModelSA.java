package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;

/**
 * Optimization of a homogeneous Markov model (model H&C) via simulated
 * annealing to fit two heparinase digest constraint sets
 */
public class HCModelSA {
    /**
     * Markov models (e.g. one for hepI and one for hepIII)
     */
    HCModel[] mm = null;
    /**
     * f[i][] is the distribution of fragment length (experimental data)
     * for digest model i
     */
    double[][] f = null;
    /**
     * for perturbations
     */
    Random rand = null;
    /**
     * chain length
     */
    int n = 0;
    /**
     * number of building blocks
     */
    int m = 0;
    /**
     * number of digests
     */
    int nd = 0;
    /**
     * enzyme specificities
     */
    CSpec[] cs = null;
    /**
     * set of building blocks
     */
    BBSet bbs = null;
    /**
     * Markov model parameters
     */
    double[][] P = null;
    /**
     * best encountered model
     */
    double[][] bestP = null;

    /**
     * buffer for this.P
     */
    double[][] buffP = null;
    /**
     * linear inequality constraints (inequality version of nonnegativity,
     * sum to 1 and balance equations)
     */
    HCModelLIC mlic = null;
    /**
     * number of perturbations at each temperature
     */
    int nPert = 200;
    /**
     * annealing schedule (temperature multiplied by alpha when decreased)
     */
    double alpha = 0.999;
    /**
     * best objective function value
     */
    double bestE = 0.;

    /**
     * Optimization of a homogeneous Markov model (model H&C) via simulated annealing
     * to fit two heparinase digest constraint sets
     * @param n BKHS chain length
     * @param bbs disaccharides and their overall proportions
     * @param specFile cleavage specificities
     * @param consFile experimental constraints (files with experimentally measured
     * distributions of fragment length)
     * @param rand
     */
    public HCModelSA(int n, BBSet bbs, String[] specFile, String[] consFile,
                     String[] outFile, String modFile, Random rand) {
        this.bbs = bbs;
        m = bbs.m;
        this.n = n;
        this.rand = rand;
        nd = specFile.length;
        cs = new CSpec[nd];
        for (int i = 0; i < nd; i++) {
            cs[i] = new CSpec(specFile[i], bbs);
        }
        loadConstraints(consFile);
        initModel();
        double t = initT();
        int nup = 1;
        double E = bestE;
        while (nup > 0) {
            t *= alpha;
            System.err.println(t + "\t" + E + "\t" + bestE + "\t" + nup);
            nup = 0;
            for (int i = 0; i < nPert; i++) {
                perturb();
                double E2 = getE();
                boolean accept = false;
                if (E2 < E) {
                    accept = true;
                } else {
                    double p = Math.exp(-Math.abs(E2 - E) / t);
                    if (rand.nextDouble() < p) {
                        accept = true;
                        nup++;
                    }
                }
                if (accept) {
                    E = E2;
                    if (E < bestE) {
                        bestE = E;
                        saveP();
                    }
                } else {
                    restore();
                }
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = bestP[i][j];
            }
        }
        mm = new HCModel[nd];
        for (int i = 0; i < nd; i++) {
            mm[i] = new HCModel(n, bbs, cs[i], P, f[i].length);
            saveFit(outFile[i], i);
        }
        saveMM(modFile);
    }

    /**
     * saves the optimized transition probability matrix
     * @param file output file for matrix of transition probabilities
     */
    void saveMM(String file) {
        try {
            FileWriter fw = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("prev\tnext\tp");
            bw.newLine();
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    bw.write(bbs.name[i] + "\t" + bbs.name[j] + "\t" + P[i][j]);
                    bw.newLine();
                }
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * saves experimental digest data and modeled digest data
     * @param file
     * @param i digest number (heparinase number)
     */
    void saveFit(String file, int i) {
        try {
            FileWriter fw = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("l\tfexp\thmod");
            bw.newLine();
            for (int l = 0; l < mm[i].lm; l++) {
                bw.write((l + 1) + "\t" + f[i][l] + "\t" + mm[i].h[l]);
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * buffers encountered optimal P
     */
    void saveP() {
        bestP = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                bestP[i][j] = P[i][j];
            }
        }
    }

    /**
     * estimation of initial temperature (pr(upward jump) = 0.6)
     * @return estimation of initial temperature (pr(upward jump) = 0.6)
     */
    double initT() {
        double res = 0.;
        bestE = getE();
        saveP();
        double dE = 0.;
        for (int i = 0; i < nPert; i++) {
            perturb();
            dE += Math.abs(bestE - getE());
            restore();
        }
        dE /= (double) nPert;
        res = -dE / Math.log(0.6);
        return res;
    }

    /**
     * perturbation of transition probability matrix P (random perturbation of one
     * entry followed by projection) after buffering the current matrix in case
     * the perturbation is not accepted
     */
    void perturb() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                buffP[i][j] = P[i][j];
            }
        }
        P[rand.nextInt(m)][rand.nextInt(m)] = rand.nextDouble();
        ProjOnPolyHSet proj = new ProjOnPolyHSet(HCModel.toVector(P), mlic.A, mlic.b, rand);
        P = HCModel.toMatrix(m, proj.optimum);
        mm = new HCModel[nd];
        for (int i = 0; i < nd; i++) {
            mm[i] = new HCModel(n, bbs, cs[i], P, f[i].length);
        }
    }

    /**
     * restoration of previous P when a perturbation was not accepted
     */
    void restore() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = buffP[i][j];
            }
        }
        mm = new HCModel[nd];
        for (int i = 0; i < nd; i++) {
            mm[i] = new HCModel(n, bbs, cs[i], P, f[i].length);
        }
    }

    /**
     * initialization of global variables and random initialization of this.P
     */
    void initModel() {
        P = new double[m][m];
        buffP = new double[m][m];
        mlic = new HCModelLIC(bbs);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                P[i][j] = rand.nextDouble();
            }
        }
        ProjOnPolyHSet proj = new ProjOnPolyHSet(HCModel.toVector(P), mlic.A, mlic.b, rand);
        P = HCModel.toMatrix(m, proj.optimum);
        mm = new HCModel[nd];
        for (int i = 0; i < nd; i++) {
            mm[i] = new HCModel(n, bbs, cs[i], P, f[i].length);
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                buffP[i][j] = P[i][j];
            }
        }
    }

    /**
     * objective function (L1 distance between modeled and experimental fragment
     * length distributions)
     * @return objective function (L1 distance between modeled and experimental
     * fragment length distributions)
     */
    double getE() {
        double res = 0.;
        for (int i = 0; i < nd; i++) {
            double d = 0.;
            for (int j = 0; j < mm[i].lm; j++) {
                d += Math.abs(mm[i].h[j] - f[i][j]);
            }
            res += d;
        }
        return res;
    }

    /**
     * reads experimental fragment length distributions from two files
     * @param file digest constraint files
     */
    void loadConstraints(String[] file) {
        f = new double[nd][];
        for (int i = 0; i < nd; i++) {
            List<String> v = new ArrayList<>();
            try {
                FileReader fr = new FileReader(file[i]);
                BufferedReader br = new BufferedReader(fr);
                String line = br.readLine();
                while ((line = br.readLine()) != null) {
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    if (st.countTokens() == 2) {
                        v.add(line);
                    }
                }
                br.close();
                fr.close();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            f[i] = new double[v.size()];
            for (int j = 0; j < v.size(); j++) {
                String line = v.get(j);
                StringTokenizer st = new StringTokenizer(line, "\t");
                int L = Integer.parseInt(st.nextToken());
                double D = Double.parseDouble(st.nextToken());
                f[i][L - 1] = D;
            }
        }
    }

    /**
     * For supplementary material
     * @param inDir input directory (ends with"\\")
     * @param outDir output directory (ends with"\\")
     */
    public static void varyingN(String inDir, String outDir) {
        Random rand = new Random(1);
        int ns = 20;
        List<String> v = new ArrayList<>();
        v.add("N\tminE\tavE\tmaxE");
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        for (int n = 12; n <= 20; n++) {
            double avE = 0.;
            double minE = 10.;
            double maxE = 0.;
            for (int s = 1; s <= ns; s++) {
                String pref = "N" + n + "s" + s;
                String[] specFile = new String[2];
                specFile[0] = inDir + "US.hepI.txt";
                specFile[1] = inDir + "US.hepIII.txt";
                String[] consFile = new String[2];
                consFile[0] = inDir + "hepI.f.txt";
                consFile[1] = inDir + "hepIII.f.txt";
                String[] outFile = new String[2];
                outFile[0] = outDir + "MM." + pref + ".hepI.fit.res";
                outFile[1] = outDir + "MM." + pref + ".hepIII.fit.res";
                String modFile = outDir + "MM." + pref + ".P.res";
                HCModelSA sa = new HCModelSA(n, bbs, specFile, consFile, outFile, modFile, rand);
                avE += sa.bestE;
                if (sa.bestE < minE) {
                    minE = sa.bestE;
                }
                if (sa.bestE > maxE) {
                    maxE = sa.bestE;
                }
            }
            avE /= (double) ns;
            String t = n + "\t" + minE + "\t" + avE + "\t" + maxE;
            v.add(t);
        }
        Utils.saveFile(v, outDir + "HCMandN.res");
    }

    /**
     * generates data used for figure of H&C fit
     * @param inDir input directory (ends with "\\")
     * @param outDir output directory (ends with "\\")
     */
    public static void defaultN(String inDir, String outDir) {
        Random rand = new Random(1);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        int n = 16;
        for (int s = 1; s <= 100; s++) {
            String pref = "s" + s;
            String[] specFile = new String[2];
            specFile[0] = inDir + "US.hepI.txt";
            specFile[1] = inDir + "US.hepIII.txt";
            String[] consFile = new String[2];
            consFile[0] = inDir + "hepI.f.txt";
            consFile[1] = inDir + "hepIII.f.txt";
            String[] outFile = new String[2];
            outFile[0] = outDir + "HC." + pref + ".hepI.fit.res";
            outFile[1] = outDir + "HC." + pref + ".hepIII.fit.res";
            String modFile = outDir + "MM." + pref + ".P.res";
            new HCModelSA(n, bbs, specFile, consFile, outFile, modFile, rand);
        }
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        String outDir = "output\\";
        defaultN(inDir, outDir);
        varyingN(inDir, outDir);
    }
}