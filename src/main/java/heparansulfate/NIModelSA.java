package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Fitting of N&I model parameters with heparinase digest constraints via simulated annealing.
 */
public class NIModelSA {
    NIModel[] pm = null;
    double[][] f = null;
    Random rand = null;
    int n = 0;
    int m = 0;
    int nd = 0;
    CSpec[] cs = null;
    BBSet bbs = null;
    double[][] gamma = null;
    double[][] gamBest = null;
    double[][] gamBuff = null;
    NIModelLIC plic = null;
    int nPert = 400;
    double alpha = 0.999;
    double bestE = 0.;

    /**
     * Constructor which adds specific constraints at positions (e.g., Reducing End).
     */
    public NIModelSA(int n, BBSet bbs, String[] specFile, String[] consFile, String[] outFile,
                     String modFile, double[][] zeta, int[] pos, Random rand) {
        this.bbs = bbs;
        this.m = bbs.m;
        this.n = n;
        this.rand = rand;
        this.nd = specFile.length;
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
                if (t < t0 / 100000.) {
                    break;
                }
            }
        }
        finalizeModel(outFile, modFile);
    }

    /**
     * Standard constructor for fitting N&I model parameters.
     */
    public NIModelSA(int n, BBSet bbs, String[] specFile, String[] consFile,
                     String[] outFile, String modFile, Random rand) {
        this.bbs = bbs;
        this.m = bbs.m;
        this.n = n;
        this.rand = rand;
        this.nd = specFile.length;
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
        finalizeModel(outFile, modFile);
    }

    private void finalizeModel(String[] outFile, String modFile) {
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

    void initModel(double[][] zeta, int[] pos) {
        gamma = new double[n][m];
        gamBuff = new double[n][m];
        plic = new NIModelLIC(n, bbs, zeta, pos);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) gamma[i][j] = rand.nextDouble();
        }
        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gamma), plic.A, plic.b, rand);
        gamma = NIModel.toMatrix(proj.optimum, n);
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
        saveGamma();
    }

    void initModel() {
        gamma = new double[n][m];
        gamBuff = new double[n][m];
        plic = new NIModelLIC(n, bbs);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) gamma[i][j] = rand.nextDouble();
        }
        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gamma), plic.A, plic.b, rand);
        gamma = NIModel.toMatrix(proj.optimum, n);
        pm = new NIModel[nd];
        for (int i = 0; i < nd; i++) pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
        saveGamma();
    }

    void savePM(String file) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            bw.write("position");
            for (int j = 0; j < m; j++) bw.write("\t" + bbs.name[j]);
            bw.newLine();
            for (int i = 0; i < n; i++) {
                bw.write(String.valueOf(i + 1));
                for (int j = 0; j < m; j++) bw.write("\t" + pm[0].gamma[i][j]);
                bw.newLine();
            }
        } catch (Exception ex) { ex.printStackTrace(); }
    }

    void saveFit(String file, int i) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            bw.write("l\tfexp\thmod");
            bw.newLine();
            for (int l = 0; l < pm[i].lm; l++) {
                bw.write((l + 1) + "\t" + f[i][l] + "\t" + pm[i].h[l]);
                bw.newLine();
            }
        } catch (Exception ex) { ex.printStackTrace(); }
    }

    void saveGamma() {
        gamBest = new double[n][m];
        for (int i = 0; i < n; i++) {
            System.arraycopy(gamma[i], 0, gamBest[i], 0, m);
        }
    }

    double initT() {
        bestE = getE();
        saveGamma();
        double dE = 0.;
        for (int i = 0; i < nPert; i++) {
            perturb();
            dE += Math.abs(bestE - getE());
            restore();
        }
        return - (dE / (double) nPert) / Math.log(0.5);
    }

    void perturb() {
        for (int i = 0; i < n; i++) System.arraycopy(gamma[i], 0, gamBuff[i], 0, m);
        gamma[1 + rand.nextInt(n - 1)][rand.nextInt(m)] = rand.nextDouble();
        ProjOnPolyHSet proj = new ProjOnPolyHSet(NIModel.toVector(gamma), plic.A, plic.b, rand);
        gamma = NIModel.toMatrix(proj.optimum, n);
        for (int i = 0; i < nd; i++) pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
    }

    void restore() {
        for (int i = 0; i < n; i++) System.arraycopy(gamBuff[i], 0, gamma[i], 0, m);
        for (int i = 0; i < nd; i++) pm[i] = new NIModel(bbs, cs[i], gamma, f[i].length);
    }

    double getE() {
        double res = 0.;
        for (int i = 0; i < nd; i++) {
            for (int j = 0; j < pm[i].lm; j++) {
                res += Math.abs(pm[i].h[j] - f[i][j]);
            }
        }
        return res;
    }

    void loadConstraints(String[] file) {
        f = new double[nd][];
        for (int i = 0; i < nd; i++) {
            Vector<String> v = new Vector<>();
            try (BufferedReader br = new BufferedReader(new FileReader(file[i]))) {
                String line = br.readLine(); // skip header
                while ((line = br.readLine()) != null) {
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    if (st.countTokens() == 2) v.add(line);
                }
            } catch (Exception ex) { ex.printStackTrace(); }
            f[i] = new double[v.size()];
            for (int j = 0; j < v.size(); j++) {
                StringTokenizer st = new StringTokenizer(v.elementAt(j), "\t");
                int L = Integer.parseInt(st.nextToken());
                f[i][L - 1] = Double.parseDouble(st.nextToken());
            }
        }
    }
}
