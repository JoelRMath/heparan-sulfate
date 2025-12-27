package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

/**
* MaxEnt model of profiles of S/U composition and transition probabilities
* along BKHS chains when BKHS is represented with a mixture of chain lengths
*/
public class MaxEntModelW {
    /**
    * solution (individual species abundances p) via geometric programming
    */
    MaxEntOptim opt = null;
    /**
    * number of disaccharides
    */
    int m = 2;
    /**
    * transition probabilities averaged over all positions, P[u][s]: from u to s
    */
    double[][] P = null;
    /**
    * disaccharides and overall proportions
    */
    BBSet bbs = null;
    /**
    * mixture model of chain lengths and sequence enumeration
    */
    MixSpecies msp = null;
    /**
    * prefix for output files
    */
    String outPref = null;
    
    /**
    * Helper class for sorting species by abundance (New Implementation)
    */
    class SpeciesEntry implements Comparable<SpeciesEntry> {
        double p;
        int[] seq;
        
        public SpeciesEntry(double p, int[] seq) {
            this.p = p;
            this.seq = seq;
        }
        
        @Override
        public int compareTo(SpeciesEntry o) {
            // Sort descending (largest probability first)
            return Double.compare(o.p, this.p);
        }
    }
    
    /**
    * MaxEnt model of profiles of S/U composition and transition probabilities
    * along BKHS chains when BKHS is represented with a mixture of chain lengths;
    * saves MaxEnt individual species abundances (*.pind.res file)
    * @param outF prefix for output files
    * @param lmin smallest chain length
    * @param lmax largest chain length
    * @param sigma spread of chain lengths
    * @param mu average chain length
    * @param inDir input directory
    */
    public MaxEntModelW(String outF, int lmin, int lmax, double sigma, double mu, String inDir) {
        outPref = outF;
        LinEqCons lec = getWLEC(lmin, lmax, sigma, mu, inDir);
        opt = new MaxEntOptim(lec.A, lec.b);
        try {
            FileWriter fw = new FileWriter(outPref + ".pind.res");
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < msp.N; i++) {
                String s = String.valueOf(opt.p[i]);
                for (int j = 0; j < msp.seq[i].length; j++) {
                    s += "\t" + msp.seq[i][j];
                }
                bw.write(s);
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    /**
    * POST-PUBLICATION CONVENIENCE CONSTRUCTOR.
    * Optimized to filter species by cumulative mass threshold to reduce file size for Git/GitHub.
    * MaxEnt model of profiles of S/U composition and transition probabilities
    * along BKHS chains when BKHS is represented with a mixture of chain lengths;
    * saves filtered MaxEnt individual species abundances (*.pind.res file)
    * @param outF prefix for output files
    * @param lmin smallest chain length
    * @param lmax largest chain length
    * @param sigma spread of chain lengths
    * @param mu average chain length
    * @param inDir input directory
    * @param massThreshold cumulative abundance threshold (e.g. 0.999)
    */
    public MaxEntModelW(String outF, int lmin, int lmax, double sigma, double mu, String inDir, double massThreshold) {
        outPref = outF;
        LinEqCons lec = getWLEC(lmin, lmax, sigma, mu, inDir);
        opt = new MaxEntOptim(lec.A, lec.b);
        
        // 1. Collect all species in memory
        List<SpeciesEntry> entries = new ArrayList<>();
        for (int i = 0; i < msp.N; i++) {
            entries.add(new SpeciesEntry(opt.p[i], msp.seq[i]));
        }
        
        // 2. Sort by abundance (descending)
        Collections.sort(entries);
        
        // 3. Write only until threshold is reached
        try {
            FileWriter fw = new FileWriter(outPref + ".pind.res");
            BufferedWriter bw = new BufferedWriter(fw);
            double currentMass = 0.0;
            int count = 0;
            
            for (SpeciesEntry e : entries) {
                String s = String.valueOf(e.p);
                for (int val : e.seq) {
                    s += "\t" + val;
                }
                bw.write(s);
                bw.newLine();
                
                currentMass += e.p;
                count++;
                if (currentMass >= massThreshold) {
                    break;
                }
            }
            bw.close();
            fw.close();
            System.out.println("Saved " + outPref + " with " + count + " species (Mass: " + currentMass + ")");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    /**
    * default linear equality constraint with chain length distribution
    * = composition-1 + 2(digest-1) + chain length distribution -1
    * @param lmin smallest chain length
    * @param lmax largest chain length
    * @param sigma chain length spread
    * @param mu average chain length
    * @param inDir input directory
    * @return default linear equality constraint with chain length distribution
    */
    LinEqCons getWLEC(int lmin, int lmax, double sigma, double mu, String inDir) {
        msp = new MixSpecies(lmin, lmax, sigma, mu);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        List<LinEqCons> v = new ArrayList<>();
        v.add(LinEqCons.removeLastRow(msp.getRhoLEC(bbs)));
        v.add(LinEqCons.removeLastRow(msp.getFragLEC(Species.loadFragAbund(inDir + "hepI.f.txt"),
        new CSpec(inDir + "US.hepI.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(msp.getFragLEC(Species.loadFragAbund(inDir + "hepIII.f.txt"),
        new CSpec(inDir + "US.hepIII.txt", bbs), bbs)));
        v.add(LinEqCons.removeLastRow(msp.getWLEC()));
        return new LinEqCons(v);
    }
    
    /**
    * reads a file containing MaxEnt abundances of individual species
    * (pref.pind.res) and save files containing:
    * NRE and RE profiles of disaccharide composition (pref.prof.res)
    * disaccharide composition as a function of chain length (pref.rhol.res)
    * profile of transition probabilities when aligning by NRE (pref.ptnre.res)
    * profile of transition probabilities when aligning by RE (pref.ptre.res)
    * @param inDir input directory (ends with "/")
    * @param pref prefix for output files
    */
    public static void makeProf(String inDir, String pref) {
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
        
        // --- RE Transition Block (Distance-from-RE Strategy) ---
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
    
    /**
    * invert an int sequence
    * @param seq int sequence to be inverted
    * @return inverse sequence of seq
    */
    public static int[] revertSeq(int[] seq){
        
        int[] rseq = new int[seq.length];
        for (int k = 0; k < rseq.length; k++) {
            rseq[k] = seq[rseq.length - 1 - k];
        }
        return(rseq);
    }
    /**
    * makes and saves individual species abundances for sigma = 3.5
    * @param inDir input directory (ends with "/")
    * @param outDir output directory (ends with "/")
    */
    public static void makeSig35(String inDir, String outDir) {
        new MaxEntModelW(outDir + "MEMW", 10, 20, 3.5, 16., inDir);
    }
    
    /**
    * POST-PUBLICATION CONVENIENCE METHOD.
    * Makes and saves individual species abundances for sigma = 3.5
    * filtering species to keep 0.999 cumulative mass.
    * Use this to reduce file sizes for version control (GitHub).
    * Output Prefix: MEMWF
    * @param inDir input directory (ends with "/")
    * @param outDir output directory (ends with "/")
    */
    public static void makeSig35Filtered(String inDir, String outDir) {
        new MaxEntModelW(outDir + "MEMWF", 10, 20, 3.5, 16., inDir, 0.9999);
    }
    
    /**
    * makes and saves individual species abundances for sigma = 1.5
    * @param inDir input directory (ends with "/")
    * @param outDir output directory (ends with "/")
    */
    public static void makeSig15(String inDir, String outDir) {
        new MaxEntModelW(outDir + "MEMWSig1.5", 10, 20, 1.5, 16., inDir);
    }
    
    /**
    * POST-PUBLICATION CONVENIENCE METHOD.
    * Makes and saves individual species abundances for sigma = 1.5
    * filtering species to keep 0.999 cumulative mass.
    * Use this to reduce file sizes for version control (GitHub).
    * Output Prefix: MEMWSig1.5F
    * @param inDir input directory (ends with "/")
    * @param outDir output directory (ends with "/")
    */
    public static void makeSig15Filtered(String inDir, String outDir) {
        new MaxEntModelW(outDir + "MEMWSig1.5F", 10, 20, 1.5, 16., inDir, 0.9999);
    }
    
    /**
    * saves distributions of chain lengths, composition as a function of chain length
    * and profiles (composition and transition probabilities) for ORIGINAL FULL FILES
    * @param inDir input directory (ends with "/")
    * @param outDir output directory (ends with "/")
    */
    public static void glAndRholAndProfExamples(String inDir, String outDir) {
        makeProf(inDir, outDir + "MEMW");
        makeProf(inDir, outDir + "MEMWSig1.5");
    }
    
    /**
    * POST-PUBLICATION CONVENIENCE METHOD.
    * Saves distributions and profiles for FILTERED FILES (suffix "F").
    * @param inDir input directory (ends with "/")
    * @param outDir output directory (ends with "/")
    */
    public static void glAndRholAndProfExamplesFiltered(String inDir, String outDir) {
        makeProf(inDir, outDir + "MEMWF");
        makeProf(inDir, outDir + "MEMWSig1.5F");
    }
    
    /**
    * Main entry point
    * @param args command line arguments
    */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        
        // Original (Full) Run - Comment out if disk space is critical
        makeSig35(inDir, outDir);
        makeSig15(inDir, outDir);
        glAndRholAndProfExamples(inDir, outDir);
        
        // Filtered Run (0.999 mass) - Optimized for GitHub storage
        makeSig35Filtered(inDir, outDir);
        makeSig15Filtered(inDir, outDir);
        glAndRholAndProfExamplesFiltered(inDir, outDir);
    }
}