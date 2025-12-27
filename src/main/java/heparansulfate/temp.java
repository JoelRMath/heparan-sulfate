package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.StringTokenizer;
import java.util.Vector;

public class temp {
    public static void makeProf(String inDir, String pref){
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
                p[i] = Double.valueOf(st.nextToken());
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
        lab[0] = "S";
        lab[1] = "U";
        
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
            int sLen = seq[s].length;
            for (int i = 0; i < sLen; i++) {
                // Indexing from the end: RE is at index 0 in the reProf array
                reProf[i][seq[s][sLen - 1 - i]] += p[s];
            }
        }
        for (int i = 0; i < lmax; i++) {
            for (int j = 0; j < m; j++) {
                reProf[i][j] /= lengthAtLeast[i];
            }
        }
        
        Vector<String> v = new Vector<String>();
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
            sr += "\t" + String.valueOf(lmax - i);
            for (int j = 0; j < m; j++) {
                sr += "\t" + reProf[i][j];
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".prof.res");
        
        int[] lenArray = new int[lmax - lmin + 1];
        for (int i = 0; i < lenArray.length; i++) {
            lenArray[i] = lmin + i;
        }
        double[] w = new double[lenArray.length];
        for (int s = 0; s < N; s++) {
            int lIdx = seq[s].length - lmin;
            w[lIdx] += p[s];
        }
        double[][] lenRho = new double[lenArray.length][m];
        for (int s = 0; s < N; s++) {
            int lIdx = seq[s].length - lmin;
            for (int i = 0; i < seq[s].length; i++) {
                lenRho[lIdx][seq[s][i]] += p[s];
            }
        }
        for (int l = 0; l < lenArray.length; l++) {
            for (int j = 0; j < m; j++) {
                lenRho[l][j] /= (double) (lenArray[l]);
                lenRho[l][j] /= w[l];
            }
        }
        v = new Vector<String>();
        sr = "l\tw\trhoS\trhoU";
        v.add(sr);
        for (int l = 0; l < w.length; l++) {
            sr = String.valueOf(lenArray[l]);
            sr += "\t" + String.valueOf(w[l]);
            sr += "\t" + String.valueOf(lenRho[l][0]); // rhoS
            sr += "\t" + String.valueOf(lenRho[l][1]); // rhoU
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".rhol.res");
        
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        
        // --- NRE Transition Block ---
        double[][][] ptnre = new double[lmax][m][m];
        for (int s = 0; s < N; s++) {
            for (int i = 0; i < m; i++) {
                ptnre[0][i][seq[s][0]] += p[s];
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.;
            for (int j = 0; j < m; j++) sum += ptnre[0][i][j];
            for (int j = 0; j < m; j++) if (sum > 0) ptnre[0][i][j] /= sum;
        }
        for (int pos = 1; pos < lmax; pos++) {
            for (int s = 0; s < N; s++) {
                if (seq[s].length > pos) {
                    ptnre[pos][seq[s][pos - 1]][seq[s][pos]] += p[s];
                }
            }
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) sum += ptnre[pos][i][j];
                for (int j = 0; j < m; j++) if (sum > 0) ptnre[pos][i][j] /= sum;
            }
        }
        // Save PTNRE block (standard forward alignment)
        v = new Vector<String>();
        sr = "pos";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) sr += "\tp" + bbs.name[i] + bbs.name[j];
        }
        v.add(sr);
        for (int pos = 0; pos < lmax; pos++) {
            sr = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) sr += "\t" + String.valueOf(ptnre[pos][i][j]);
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".ptnre.res");

        // --- RE Transition Block ---
        double[][][] ptre = new double[lmax][m][m];
        for (int s = 0; s < N; s++) {
            int[] cseq = completeSeq(seq[s], lmax);
            if (cseq[0] >= 0){
                for (int i = 0; i < m; i++) {
                    ptre[0][i][cseq[0]] += p[s];
                }
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.;
            for (int j = 0; j < m; j++) sum += ptre[0][i][j];
            for (int j = 0; j < m; j++) if (sum > 0) ptre[0][i][j] /= sum;
        }
        for (int pos = 1; pos < lmax; pos++) {
            for (int s = 0; s < N; s++) {
                int[] cseq = completeSeq(seq[s], lmax);
                if (cseq[pos-1] >= 0) {
                    ptre[pos][cseq[pos - 1]][cseq[pos]] += p[s];
                }
            }
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) sum += ptre[pos][i][j];
                for (int j = 0; j < m; j++) if (sum > 0) ptre[pos][i][j] /= sum;
            }
        }
        // Save PTRE block (standard forward alignment)
        v = new Vector<String>();
        sr = "pos";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) sr += "\tp" + bbs.name[i] + bbs.name[j];
        }
        v.add(sr);
        for (int pos = 0; pos < lmax; pos++) {
            sr = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) sr += "\t" + String.valueOf(ptre[pos][i][j]);
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".ptre.res");
        
/*         // --- RE Transition Block (Distance-from-RE Strategy) ---
        double[][][] ptre = new double[lmax][m][m];
        
        for (int s = 0; s < N; s++) {
            int sLen = seq[s].length;
            
            // 1. Initial state at the RE (The Seed)
            if (sLen > 0) {
                int reVal = seq[s][sLen - 1];
                for (int i = 0; i < m; i++) {
                    ptre[0][i][reVal] += p[s];
                }
            }
            
            // 2. Transitions: p(current | previous) indexed by distance from RE
            for (int i = 1; i < sLen; i++) {
                // Distance of the current unit 'i' from the Reducing End
                int distFromRE = (sLen - 1) - i; 
                
                // We want distFromRE=0 (the transition ending at RE) to be pos=1
                int posIdx = distFromRE + 1;
                
                if (posIdx < lmax) {
                    // We keep the forward transition identity: p(seq[i] | seq[i-1])
                    ptre[posIdx][seq[s][i-1]][seq[s][i]] += p[s];
                }
            }
        }
        
        // Normalize RE positions
        for (int pos = 0; pos < lmax; pos++) {
            for (int i = 0; i < m; i++) {
                double sum = 0.;
                for (int j = 0; j < m; j++) sum += ptre[pos][i][j];
                for (int j = 0; j < m; j++) {
                    if (sum > 0) ptre[pos][i][j] /= sum;
                }
            }
        }
 */        
        // Save PTRE
        v = new Vector<String>();
        sr = "pos";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) sr += "\tp" + bbs.name[i] + bbs.name[j];
        }
        v.add(sr);
        for (int pos = 0; pos < lmax; pos++) {
            sr = String.valueOf(pos + 1);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    sr += "\t" + String.valueOf(ptre[pos][i][j]);
                }
            }
            v.add(sr);
        }
        Utils.saveFile(v, pref + ".ptre.res");
    }

    /**
    * Elongate/complete an int sequence to length lmax by filling the first positions with -1
    * @param seq int sequence to be elongated
    * @return elongated sequence of seq
    */
    public static int[] completeSeq(int[] seq, int lmax){
        
        int[] cseq = new int[lmax];
        int len = seq.length;
        int c = -1;
        for (int i = 0; i <= (lmax-len); i++){
            cseq[i] = -1;
            c++;
        }
        for (int i = c ; i < lmax; i++){
            cseq[i] = seq[i-(lmax-len)];
        }
        return(cseq);
    }

    public static void main(String[] args){

        int[] seq = {1,2,3};
        System.out.println(Arrays.toString(seq));
        int lmax = 5;
        int[] cseq = completeSeq(seq, lmax);
        System.out.println(Arrays.toString(cseq));
    }
}
