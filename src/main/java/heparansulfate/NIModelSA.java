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
 * Fits parameters of the N&amp;I model using heparinase digest constraints via simulated annealing.
 */
public class NIModelSA {
    /**
     * Position models (e.g., one for HepI and one for HepIII).
     */
    NIModel[] pm = null;
    /**
     * Experimental fragment length distribution {@code f[i][]} for digest model {@code i}.
     */
    double[][] f = null;
    /**
     * Random number generator for perturbations.
     */
    Random rand = null;
    /**
     * BKHS chain length.
     */
    int n = 0;
    /**
     * Number of building block types.
     */
    int m = 0;
    /**
     * Number of enzymatic digests.
     */
    int nd = 0;
    /**
     * Enzyme specificities and yields.
     */
    CSpec[] cs = null;
    /**
     * Set of building blocks.
     */
    BBSet bbs = null;
    /**
     * Position model parameters (Gamma matrix).
     */
    double[][] gamma = null;
    /**
     * Best model encountered during the annealing process.
     */
    double[][] gamBest = null;
    /**
     * Buffer for storing {@code this.gamma}.
     */
    double[][] gamBuff = null;
    /**
     * Linear inequality constraints including building block proportions, sum-to-1, and nonnegativity.
     */
    NIModelLIC plic = null;
    /**
     * Number of perturbations performed at each temperature step.
     */
    int nPert = 100;
    /**
     * Annealing schedule factor for geometric temperature decrease.
     */
    double alpha = 0.99;
    /**
     * Best encountered value of the objective function (energy).
     */
    double bestE = 0.;

    /**
     * Fits parameters of the N&amp;I model using simulated annealing with additional constraints.
     * 
     * @param n BKHS chain length.
     * @param bbs Disaccharides and their overall proportions.
     * @param specFile Files containing cleavage specificities.
     * @param consFile Files containing experimental fragment length distributions.
     * @param outFile Output files for modeled distributions.
     * @param modFile Output file for the optimized Gamma matrix.
     * @param zeta Specified compositions at specific positions (additional constraints).
     * @param pos Indices of positions with specified compositions.
     * @param rand Random number generator.
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
     * Initializes the model when additional positional constraints are provided.
     * @param zeta Positional disaccharide compositions.
     * @param pos Indices of the specified positions.
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
     * Fits parameters of a N&amp;I model using simulated annealing.
     * @param n BKHS chain length.
     * @param bbs Disaccharides and their overall proportions.
     * @param specFile Cleavage specificity files.
     * @param consFile Experimental constraint files.
     * @param outFile Output files for fit results.
     * @param modFile Output file for the optimized Gamma matrix.
     * @param rand Random number generator.
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
     * Saves the optimized Gamma matrix of the N&amp;I model to a file.
     * @param file Output file path.
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
     * Saves experimental and modeled fragment length distributions for a heparinase index.
     * @param file Output file path.
     * @param i Heparinase index.
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
     * Buffers the current {@code gamma} matrix into {@code gamBest}.
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
     * Estimates initial temperature for annealing (based on a 0.5 upward move probability).
     * @return The estimated initial temperature.
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
     * Performs a random perturbation of one {@code gamma} entry followed by a projection.
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
     * Restores the previous {@code gamma} matrix when a perturbation is rejected.
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
     * Initializes variables and sets {@code gamma} to a projected random matrix.
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
     * Computes the objective function: the L1 distance between experimental and modeled distributions.
     * @return The total L1 error value.
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
     * Loads experimental fragment length distributions from files.
     * @param file Array of input file paths.
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
     * Fits N&amp;I model without added constraint at the reducing end.
     * @param inDir Input directory.
     * @param outDir Output directory.
     */
    public static void withoutREC(String inDir, String outDir) {
        Random rand = new Random(1);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        int n = 16;
        for (int s = 1; s <= 100; s++) {
            System.out.println(s);
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
     * Fits N&amp;I model without reducing end constraints for different chain lengths.
     * @param inDir Input directory.
     * @param outDir Output directory.
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
     * Fits N&amp;I model with added reducing end constraints.
     * @param inDir Input directory.
     * @param outDir Output directory.
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
     * Main entry point for simulated annealing optimization.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/NI/";
        if (args.length >= 1){
            inDir = args[0];
        }
        if (args.length >= 2){
            outDir = args[1];
        }
        withoutREC(inDir, outDir);
    }
}