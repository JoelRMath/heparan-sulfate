package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

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
    Map<String, Integer> name2i = null;

    /**
     * creates a set of building blocks (names and relative abundances) from a file
     * @param file ascii tab-delimited with one header row: name \t abundance
     */
    public BBSet(String file) {
        List<String> v = new ArrayList<>();
        try {
            FileReader fr = new FileReader(file);
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
        m = v.size();
        name = new String[m];
        rho = new double[m];
        name2i = new HashMap<>();
        for (int i = 0; i < m; i++) {
            String s = v.get(i);
            StringTokenizer st = new StringTokenizer(s, "\t");
            name[i] = st.nextToken().trim().toLowerCase();
            name2i.put(name[i], i);
            rho[i] = Double.parseDouble(st.nextToken().trim());
        }
        double sum = 0.;
        for (int i = 0; i < m; i++) {
            sum += rho[i];
        }
        for (int i = 0; i < m; i++) {
            rho[i] /= sum;
        }
    }

    /**
     * for testing
     * @param args command line arguments
     */
    public static void main(String[] args) {
        String file = "input/US.ab.txt";
        BBSet bs = new BBSet(file);
        for (int i = 0; i < bs.m; i++) {
            System.out.println(bs.name[i] + "\t" + bs.rho[i]);
        }
    }
}