package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Optimization of a homogeneous Markov model (model H&C) via simulated
 * annealing to fit two heparinase digest constraint sets.
 */
public class HCModelSA {
    HCModel[] mm = null;
    double[][] f = null;
    Random rand = null;
    int n = 0;
    int m = 0;
    int nd = 0;
    CSpec[] cs = null;
    BBSet bbs = null;
    double[][] P = null;
    double[][] bestP = null;
    double[][] buffP = null;
    HCModelLIC mlic = null;
    int nPert = 200;
    double alpha = 0.999;
    double bestE = 0.;

    public HCModelSA(int n, BBSet bbs, String[] specFile, String[] consFile,
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
            System.err.println("Temp: " + t + " Current E: " + E + " Best E: " + bestE + " Upward: " + nup);
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
        
        // Finalize with best found P
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

    void saveMM(String file) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            bw.write("prev\tnext\tp");
            bw.newLine();
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    bw.write(bbs.name[i] + "\t" + bbs.name[j] + "\t" + P[i][j]);
                    bw.newLine();
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void saveFit(String file, int i) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            bw.write("l\tfexp\thmod");
            bw.newLine();
            for (int l = 0; l < mm[i].lm; l++) {
                bw.write((l + 1) + "\t" + f[i][l] + "\t" + mm[i].h[l]);
                bw.newLine();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void saveP() {
        bestP = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                bestP[i][j] = P[i][j];
            }
        }
    }

    double initT() {
        bestE = getE();
        saveP();
        double dE = 0.;
        for (int i = 0; i < nPert; i++) {
            perturb();
            dE += Math.abs(bestE - getE());
            restore();
        }
        dE /= (double) nPert;
        return -dE / Math.log(0.6);
    }

    void perturb() {
        for (int i = 0; i < m; i++) {
            System.arraycopy(P[i], 0, buffP[i], 0, m);
        }
        P[rand.nextInt(m)][rand.nextInt(m)] = rand.nextDouble();
        ProjOnPolyHSet proj = new ProjOnPolyHSet(HCModel.toVector(P), mlic.A, mlic.b, rand);
        P = HCModel.toMatrix(m, proj.optimum);
        for (int i = 0; i < nd; i++) {
            mm[i] = new HCModel(n, bbs, cs[i], P, f[i].length);
        }
    }

    void restore() {
        for (int i = 0; i < m; i++) {
            System.arraycopy(buffP[i], 0, P[i], 0, m);
        }
        for (int i = 0; i < nd; i++) {
            mm[i] = new HCModel(n, bbs, cs[i], P, f[i].length);
        }
    }

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
        saveP();
    }

    double getE() {
        double res = 0.;
        for (int i = 0; i < nd; i++) {
            for (int j = 0; j < mm[i].lm; j++) {
                res += Math.abs(mm[i].h[j] - f[i][j]);
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
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            f[i] = new double[v.size()];
            for (int j = 0; j < v.size(); j++) {
                StringTokenizer st = new StringTokenizer(v.elementAt(j), "\t");
                int L = Integer.parseInt(st.nextToken());
                double D = Double.parseDouble(st.nextToken());
                f[i][L - 1] = D;
            }
        }
    }
}