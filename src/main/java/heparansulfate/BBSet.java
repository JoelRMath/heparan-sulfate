package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Set of building blocks (e.g. S and U), including names and relative abundances
 */
public class BBSet {
    /**
     * number of building blocks
     */
    int m = 0;
    /**
     * building block labels, lower case
     */
    String[] name = null;
    /**
     * building block relative abundances, they must sum to 1
     */
    double[] rho = null;
    /**
     * maps name (lower case) to index in this.name
     */
    Hashtable<String, Integer> name2i = null;

    /**
     * creates a set of building blocks (names and relative abundances) from a file
     * * @param file ascii tab-delimited with one header row: name \t abundance
     */
    public BBSet(String file) {
        Vector<String> v = new Vector<String>();
        try {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine(); // skip header
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
        name = new String[m];
        rho = new double[m];
        name2i = new Hashtable<String, Integer>();
        for (int i = 0; i < m; i++) {
            String s = v.elementAt(i);
            StringTokenizer st = new StringTokenizer(s, "\t");
            name[i] = st.nextToken().trim().toLowerCase();
            // Modern Java prefers valueOf over 'new Integer'
            name2i.put(name[i], Integer.valueOf(i));
            // Modern Java prefers Double.parseDouble over 'new Double'
            rho[i] = Double.parseDouble(st.nextToken().trim());
        }
        double sum = 0.;
        for (int i = 0; i < m; i++) {
            sum += rho[i];
        }
        if (sum > 0) {
            for (int i = 0; i < m; i++) {
                rho[i] /= sum;
            }
        }
    }

    /**
     * for testing
     * * @param args
     */
    public static void main(String[] args) {
        // Mac uses forward slashes. Ensure this file exists in your project root!
        String file = "input/US.ab.txt"; 
        BBSet bs = new BBSet(file);
        for (int i = 0; i < bs.m; i++) {
            System.out.println(bs.name[i] + "\t" + bs.rho[i]);
        }
    }
}