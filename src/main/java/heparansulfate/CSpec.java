package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Cleavage specificities, based on a BBSet object (set of building blocks)
 * and a file giving specificity for each building block name.
 */
public class CSpec {
    /**
     * number of building blocks
     */
    int m = 0;
    /**
     * cleavage specificities
     */
    double[] c = null;

    /**
     * Creates cleavage specificities, based on a BBSet object (set of building blocks) and
     * a file giving specificity for each building block name.
     * @param file ascii tab-delimited file with one header row: name \t specificity
     * @param bbset set of building blocks
     */
    public CSpec(String file, BBSet bbset) {
        Vector<String> v = new Vector<String>();
        try {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine(); // Skip header
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
        
        m = v.size();
        c = new double[m];
        for (int i = 0; i < m; i++) {
            String s = v.elementAt(i);
            StringTokenizer st = new StringTokenizer(s, "\t");
            String name = st.nextToken().trim().toLowerCase();
            
            Integer I = bbset.name2i.get(name);
            if (I != null) {
                // Modern Java: Double.parseDouble is preferred over new Double()
                c[I] = Double.parseDouble(st.nextToken().trim());
            } else {
                System.err.println("Warning: Building block '" + name + "' not found in BBSet.");
            }
        }
    }

    /**
     * for testing
     * @param args
     */
    public static void main(String[] args) {
        // Updated paths for Mac compatibility (Forward slashes)
        String inDir = "input/";
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        
        file = inDir + "US.hepI.txt";
        CSpec cs = new CSpec(file, bbs);
        System.out.println("Hep I Specificities:");
        for (int i = 0; i < cs.m; i++) {
            System.out.println(cs.c[i]);
        }
        
        file = inDir + "US.hepIII.txt";
        cs = new CSpec(file, bbs);
        System.out.println("Hep III Specificities:");
        for (int i = 0; i < cs.m; i++) {
            System.out.println(cs.c[i]);
        }
    }
}