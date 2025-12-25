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
 * Fitting of N&I model parameters with heparinase digest constraints via simulated annealing
 */
public class NIModelSA {
    /**
     * position models (e.g. one for hepI and one for hepIII)
     */
    NIModel[] pm = null;
    /**
     * f[i][] is the distribution of fragment length (experimental data) for digest model i
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
     * enzyme specificities/yields
     */
    CSpec[] cs = null;
    /**
     * set of building blocks
     */
    BBSet bbs = null;
    /**
     * position model parameters
     */
    double[][] gamma = null;
    /**
     * best encountered model
     */
    double[][] gamBest = null;
    /**
     * buffer for this.gamma
     */
    double[][] gamBuff = null;
    /**
     * linear inequality constraints provided by overall building block proportions,
     * sum to 1 at each position and nonnegativity
     */
    NIModelLIC plic = null;
    /**
     * number of perturbations at each temperature
     */
    int nPert = 400;
    /**
     * annealing schedule (geometric temperature decrease)
     */
    double alpha = 0.999;
    /**
     * best encountered value of the objective function
     */
    double bestE = 0.;

    /**
     * Fitting of N&I model parameters with heparinase digest constraints via simulated annealing;
     * special constructor which adds constraints
     * @param n BKHS chain length
     * @param bbs disachharides and their overall proportions
     * @param specFile cleavage specificities/yields
     * @param consFile constraints (experimental fragment length distributions)
     * @param outFile output file for modeled fragment length distributions
     * @param modFile output file for optimized matrix Gamma of N&I model
     * @param zeta specified compositions at positions in pos (additional constraints)
     * @param pos positions at which composition is further specified (additional constraints)
     * @param rand
     */
    public NIModelSA(int n, BBSet bbs, String[] specFile, String[] consFile, String[] outFile,
                     String modFile, double[][] zeta, int[] pos, Random rand) {
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
        initModel(zeta, pos);
        double t = initT();
        double t0 = t;
        int nup = 1;
        double E = bestE;
        while (nup > 0) {
            t *= alpha;
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
                        saveGamma();
                    }
                } else {
                    restore();
                }
            }
            if (t < t0 / 100000.) {
                break;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamma[i][j] = gamBest[i][j];
            }
        }
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) {
            pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
            saveFit(outFile[i], i);
        }
        savePM(modFile);
    }

    /**
     * initialization when additional constraints
     * (specified compositions zeta in positions pos) are given
     * @param zeta disaccharide composition
     * @param pos position in chains
     */
    void initModel(double[][] zeta, int[] pos) {
        gamma = new double[n][m];
        gamBuff = new double[n][m];
        plic = new NIModelLIC(n, bbs, zeta, pos);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamma[i][j] = rand.nextDouble();
            }
        }
        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gamma), plic.A, plic.b, rand);
        gamma = NIModel.toMatrix(proj.optimum, n);
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) {
            pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamBuff[i][j] = gamma[i][j];
            }
        }
    }

    /**
     * Fitting of parameters of a N&I model with heparinase digest constraints
     * via simulated annealing
     * @param n BKHS chain length
     * @param bbs disachharides and their overall proportions
     * @param specFile cleavage specificities/yields
     * @param consFile constraints (experimental fragment length distributions)
     * @param outFile output file for modeled fragment length distributions
     * @param modFile output file for optimized matrix Gamma of N&I model
     * @param rand
     */
    public NIModelSA(int n, BBSet bbs, String[] specFile, String[] consFile,
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
                        saveGamma();
                    }
                } else {
                    restore();
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamma[i][j] = gamBest[i][j];
            }
        }
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) {
            pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
            saveFit(outFile[i], i);
        }
        savePM(modFile);
    }

    /**
     * saves matrix Gamma of the optimized N&I model
     * @param file
     */
    void savePM(String file) {
        try {
            FileWriter fw = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("position");
            for (int j = 0; j < m; j++) {
                bw.write("\t" + bbs.name[j]);
            }
            bw.newLine();
            for (int i = 0; i < n; i++) {
                bw.write(String.valueOf(i + 1));
                for (int j = 0; j < m; j++) {
                    bw.write("\t" + pm[0].gamma[i][j]);
                }
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * saves experimental and modeled (optimized) distribution of fragment length
     * after digestion by heparinase i
     * @param file output file
     * @param i heparinase number
     */
    void saveFit(String file, int i) {
        try {
            FileWriter fw = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("l\tfexp\thmod");
            bw.newLine();
            for (int l = 0; l < pm[i].lm; l++) {
                bw.write((l + 1) + "\t" + f[i][l] + "\t" + pm[i].h[l]);
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * buffers gamma into gamBest
     */
    void saveGamma() {
        gamBest = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamBest[i][j] = gamma[i][j];
            }
        }
    }

    /**
     * estimates initial temperature (pr(upward move) = 0.5)
     * @return initial temperature
     */
    double initT() {
        double res = 0.;
        bestE = getE();
        saveGamma();
        double dE = 0.;
        for (int i = 0; i < nPert; i++) {
            perturb();
            dE += Math.abs(bestE - getE());
            restore();
        }
        dE /= (double) nPert;
        res = -dE / Math.log(0.5);
        return res;
    }

    /**
     * one random perturbation of one gamma entry and projection,
     * after buffering current gamma
     */
    void perturb() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamBuff[i][j] = gamma[i][j];
            }
        }
        gamma[1 + rand.nextInt(n - 1)][rand.nextInt(m)] = rand.nextDouble();
        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gamma), plic.A, plic.b, rand);
        gamma = NIModel.toMatrix(proj.optimum, n);
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) {
            pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
        }
    }

    /**
     * restores previous gamma when a perturbation was not accepted
     */
    void restore() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamma[i][j] = gamBuff[i][j];
            }
        }
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) {
            pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
        }
    }

    /**
     * initialization of variables; gamma is set to a random matrix and projected
     */
    void initModel() {
        gamma = new double[n][m];
        gamBuff = new double[n][m];
        plic = new NIModelLIC(n, bbs);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamma[i][j] = rand.nextDouble();
            }
        }
        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gamma), plic.A, plic.b, rand);
        gamma = NIModel.toMatrix(proj.optimum, n);
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) {
            pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                gamBuff[i][j] = gamma[i][j];
            }
        }
    }

    /**
     * returns the objective function: L1 distance between experimental
     * and modeled fragment length distributions
     * @return L1 distance between experimental and modeled fragment
     * length distributions
     */
    double getE() {
        double res = 0.;
        for (int i = 0; i < nd; i++) {
            double d = 0.;
            for (int j = 0; j < pm[i].lm; j++) {
                d += Math.abs(pm[i].h[j] - f[i][j]);
            }
            res += d;
        }
        return res;
    }

    /**
     * loads experimental fragment length distributions
     * @param file input file
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
     * fits N&I model without added constraint at the reducing end;
     * for supplementary material
     * @param inDir input directory (ends with "\\")
     * @param outDir output directory (ends with "\\")
     */
    public static void withoutREC(String inDir, String outDir) {
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
            outFile[0] = outDir + "PM." + pref + ".hepI.fit.res";
            outFile[1] = outDir + "PM." + pref + ".hepIII.fit.res";
            String modFile = outDir + "PM." + pref + ".gamma.res";
            new NIModelSA(n, bbs, specFile, consFile, outFile, modFile, rand);
        }
    }

    /**
     * fits N&I model without added constraint at the reducing end
     * and for different values of BKHS chain length n; for supplementary material
     * @param inDir input directory (ends with "\\")
     * @param outDir output directory (ends with "\\")
     */
    public static void withoutRECandN(String inDir, String outDir) {
        Random rand = new Random(1);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        for (int n = 12; n <= 20; n++) {
            for (int s = 1; s <= 10; s++) {
                String pref = "N" + n + "s" + s;
                String[] specFile = new String[2];
                specFile[0] = inDir + "US.hepI.txt";
                specFile[1] = inDir + "US.hepIII.txt";
                String[] consFile = new String[2];
                consFile[0] = inDir + "hepI.f.txt";
                consFile[1] = inDir + "hepIII.f.txt";
                String[] outFile = new String[2];
                outFile[0] = outDir + "PM." + pref + ".hepI.fit.res";
                outFile[1] = outDir + "PM." + pref + ".hepIII.fit.res";
                String modFile = outDir + "PM." + pref + ".gamma.res";
                new NIModelSA(n, bbs, specFile, consFile, outFile, modFile, rand);
            }
        }
    }

    /**
     * fits N&I model with added constraint at the reducing end
     * and for different values of BKHS chain length n; for supplementary material
     * @param inDir input directory (ends with "\\")
     * @param outDir output directory (ends with "\\")
     */
    public static void withREC(String inDir, String outDir) {
        Random rand = new Random(1);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        int n = 16;
        int[] pos = new int[1];
        pos[0] = n - 1;
        int m = bbs.m;
        double[][] zeta = new double[pos.length][m];
        zeta[0][0] = 0.44;
        zeta[0][1] = 0.56;
        for (int s = 1; s <= 50; s++) {
            String pref = "s" + s;
            String[] specFile = new String[2];
            specFile[0] = inDir + "US.hepI.txt";
            specFile[1] = inDir + "US.hepIII.txt";
            String[] consFile = new String[2];
            consFile[0] = inDir + "hepI.f.txt";
            consFile[1] = inDir + "hepIII.f.txt";
            String[] outFile = new String[2];
            outFile[0] = outDir + "PMRE." + pref + ".hepI.fit.res";
            outFile[1] = outDir + "PMRE." + pref + ".hepIII.fit.res";
            String modFile = outDir + "PMRE." + pref + ".gamma.res";
            new NIModelSA(n, bbs, specFile, consFile, outFile, modFile, zeta, pos, rand);
        }
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        String outDir = "output\\";
        withoutREC(inDir, outDir);
        withoutRECandN(inDir, outDir);
        withREC(inDir, outDir);
    }
}