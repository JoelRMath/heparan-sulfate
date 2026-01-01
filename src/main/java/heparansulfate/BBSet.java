package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * Represents a set of chemical building blocks (e.g., disaccharides S and U),
 * storing their unique labels and their relative abundances. 
 * <p>
 * This class handles the ingestion of building block data and ensures that 
 * abundances are normalized to sum to 1.0.
 */
public class BBSet {
    /**
     * Number of building blocks
     */
    int m = 0;
    /** * Array of building block labels (e.g., "s", "u"), stored in lower case. 
     */
    String[] name = null;
    /**
     * Building block relative abundances. Their sum is normalized to 1 during construction.
     */
    double[] rho = null;
    /** * A lookup map connecting lower-case labels to their respective index 
     * in the {@code name} and {@code rho} arrays.
     */
    Map<String, Integer> name2i = null;

    /**
     * Constructs a set of building blocks from a data file.
     * <p>
     * The input file must be an ASCII tab-delimited file with a single header
     * row. Each subsequent row should follow the format: {@code name \t abundance}.
     * Labels are automatically converted to lower case, and abundances are 
     * normalized relative to the total sum of the input values.
     * * @param file The path to the tab-delimited file containing building block data.
     * @param file ASCII tab-delimited file containing building block names and proportions.
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
     * For testing
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