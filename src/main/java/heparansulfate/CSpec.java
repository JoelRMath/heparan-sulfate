package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * cleavage specificities, based on a BBSet object (set of building blocks)
 * and a file giving specificity for each building block name
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
     * a file giving specificity for each building block name
     * @param file ascii tab-delimited file with one header row: name \t specificity
     * @param bbset set of building blocks
     */
    public CSpec(String file, BBSet bbset) {
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
        c = new double[m];
        for (int i = 0; i < m; i++) {
            String s = v.get(i);
            StringTokenizer st = new StringTokenizer(s, "\t");
            String name = st.nextToken().trim().toLowerCase();
            Integer I = bbset.name2i.get(name);
            c[I] = Double.parseDouble(st.nextToken().trim());
        }
    }

    /**
     * for testing
     * @param args command line arguments
     */
    public static void main(String[] args) {
        String inDir = "input\\";
        String file = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(file);
        file = inDir + "US.hepI.txt";
        CSpec cs = new CSpec(file, bbs);
        for (int i = 0; i < cs.m; i++) {
            System.out.println(cs.c[i]);
        }
        file = inDir + "US.hepIII.txt";
        cs = new CSpec(file, bbs);
        for (int i = 0; i < cs.m; i++) {
            System.out.println(cs.c[i]);
        }
    }
}