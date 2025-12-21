package heparansulfate;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * enumeration of all possible sequences of length n with m disaccharides
 */
public class Species {
    int n = 0;
    int m = 0;
    public int N = 0;
    public int[][] seq = null;

    /**
     * enumeration of all possible sequences of length n with m disaccharides
     * @param m number of disaccharides
     * @param n BKHS chain length
     */
    public Species(int m, int n) {
        this.m = m;
        this.n = n;
        N = (int) Math.pow((double) m, (double) n);
        seq = new int[N][n];
        int[] buff = new int[n];
        int pos = 0;
        for (int s = 1; s < N; s++) {
            pos = 0;
            boolean kg = true;
            while (kg) {
                buff[pos]++;
                if (buff[pos] < m) {
                    kg = false;
                } else {
                    buff[pos] = 0;
                    pos++;
                }
            }
            for (int i = 0; i < n; i++) {
                seq[s][i] = buff[i];
            }
        }
    }

    public LinEqCons getHomogeneityLEC(BBSet bbs) {
        double[][] A = new double[(n - 1) * (m - 1)][N];
        double[] b = new double[(n - 1) * (m - 1)];
        String[] tp = new String[(n - 1) * (m - 1)];
        int r = -1;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < m - 1; j++) {
                r++;
                b[r] = bbs.rho[j];
                tp[r] = "hom";
                for (int s = 0; s < N; s++) {
                    if (seq[s][i] == j) {
                        A[r][s] = 1.;
                    }
                }
            }
        }
        return new LinEqCons(A, b, tp);
    }

    public LinEqCons getFragLEC(double[] f, CSpec cs, BBSet bbs) {
        double[][] A = new double[f.length][N];
        double[] b = new double[f.length];
        String[] lab = new String[f.length];
        for (int i = 0; i < b.length; i++) {
            b[i] = f[i];
            lab[i] = "frag";
        }
        int lm = f.length;
        for (int l = 0; l < lm - 1; l++) {
            for (int s = 0; s < N; s++) {
                A[l][s] = getFragContrib(s, l + 1, cs);
            }
        }
        for (int s = 0; s < N; s++) {
            A[lm - 1][s] = 0.;
            for (int l = lm - 1; l < n - 1; l++) {
                A[lm - 1][s] += getFragContrib(s, l + 1, cs);
            }
        }
        double denom = 0.;
        for (int i = 0; i < bbs.m; i++) {
            denom += bbs.rho[i] * cs.c[i];
        }
        denom *= (double) (n - 1);
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                A[i][j] /= denom;
            }
        }
        return new LinEqCons(A, b, lab);
    }

    double getFragContrib(int s, int ll, CSpec cs) {
        double res = 0.;
        if (ll > n) return res;
        for (int i = 1; i < n - ll; i++) {
            double t = cs.c[seq[s][i]];
            for (int j = 1; j <= ll - 1; j++) {
                t *= (1. - cs.c[seq[s][i + j]]);
            }
            t *= cs.c[seq[s][i + ll]];
            res += t;
        }
        double t = cs.c[seq[s][n - ll]];
        for (int j = n - ll + 1; j < n; j++) {
            t *= (1. - cs.c[seq[s][j]]);
        }
        res += t;
        return res;
    }

    public LinEqCons getCompLEC(double[] rho) {
        int m_local = rho.length;
        double[][] A = new double[m_local][N];
        double[] b = new double[m_local];
        String[] lab = new String[m_local];
        for (int i = 0; i < m_local; i++) {
            b[i] = rho[i];
            lab[i] = "comp";
            for (int s = 0; s < N; s++) {
                for (int j = 1; j < n; j++) {
                    if (seq[s][j] == i) A[i][s] += 1.;
                }
                A[i][s] /= (double) (n - 1);
            }
        }
        return new LinEqCons(A, b, lab);
    }

    public LinEqCons getNormLEC() {
        double[][] A = new double[1][N];
        double[] b = {1.};
        String[] lab = {"norm"};
        for (int i = 0; i < N; i++) A[0][i] = 1.;
        return new LinEqCons(A, b, lab);
    }

    public double[] getP(double[][] gamma) {
        double[] res = new double[N];
        for (int s = 0; s < N; s++) {
            res[s] = 1.;
            for (int i = 0; i < n; i++) {
                res[s] *= gamma[i][seq[s][i]];
            }
        }
        return res;
    }

    /**
     * Fixed: Updated Vector to List to match updated LinEqCons constructor
     */
    public LinEqCons getCompleteLEC(BBSet bbs, String[] csFile, String[] fragFile) {
        List<LinEqCons> v = new ArrayList<LinEqCons>();
        v.add(getNormLEC());
        v.add(LinEqCons.removeLastRow(getCompLEC(bbs.rho)));
        for (int i = 0; i < csFile.length; i++) {
            CSpec cs = new CSpec(csFile[i], bbs);
            double[] f = loadFragAbund(fragFile[i]);
            v.add(LinEqCons.removeLastRow(getFragLEC(f, cs, bbs)));
        }
        return new LinEqCons(v);
    }

    /**
     * Fixed: Updated Vector to List to match updated Utils return type
     */
    public static double[][] loadGamma(String file) {
        List<String> v = Utils.loadFileNoheader(file);
        StringTokenizer st = new StringTokenizer(v.get(0), "\t");
        int m_local = st.countTokens() - 1;
        int n_local = v.size();
        double[][] res = new double[n_local][m_local];
        for (int i = 0; i < n_local; i++) {
            st = new StringTokenizer(v.get(i), "\t");
            st.nextToken(); // skip position
            int j = 0;
            while (st.hasMoreTokens()) {
                res[i][j++] = Double.parseDouble(st.nextToken());
            }
        }
        return res;
    }

    /**
     * Fixed: Updated Vector to List to match updated Utils return type
     */
    public static double[] loadFragAbund(String file) {
        List<String> v = Utils.loadFileNoheader(file);
        int n_local = v.size();
        double[] res = new double[n_local];
        for (int i = 0; i < n_local; i++) {
            StringTokenizer st = new StringTokenizer(v.get(i), "\t");
            st.nextToken(); // skip length label
            res[i] = Double.parseDouble(st.nextToken());
        }
        return res;
    }
}