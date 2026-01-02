package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
* Cleavage specificities, based on a {@code BBSet} object (set of building blocks)
* and a file giving specificity for each building block name.
*/
public class CSpec {
    /**
    * Number of building blocks.
    */
    int m = 0;
    /**
    * Cleavage specificities.
    */
    double[] c = null;
    
    /**
    * Creates cleavage specificities, based on a {@code BBSet} object (set of building blocks) and
    * a file giving specificity for each building block name.
    * @param file ASCII tab-delimited file with one header row: {@code name \t specificity}.
    * @param bbset Set of building blocks.
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
    
    /** Returns a string representation of the cleavage specificities for 
    * documentation and R integration.
    * @param inDir The absolute path to the input directory.
    * @return A formatted string of specificity values.
    */
    public static String getHeparanSulfateSpecSummary(String inDir) {
        StringBuilder sb = new StringBuilder();
        
        String abFile = inDir + "US.ab.txt";
        BBSet bbs = new BBSet(abFile);
        String hepI = inDir + "US.hepI.txt";
        CSpec csI = new CSpec(hepI, bbs);
        sb.append("Heparinase I Specificities:\n");
        for (double val : csI.c) {
            sb.append(val).append("\n");
        }
        String hepIII = inDir + "US.hepIII.txt";
        CSpec csIII = new CSpec(hepIII, bbs);
        sb.append("\nHeparinase III Specificities:\n");
        for (double val : csIII.c) {
            sb.append(val).append("\n");
        }        
        return sb.toString();
    }
    /**
    * For testing.
    * @param args Command line arguments.
    */
    public static void main(String[] args) {
        String inDir = args[0];
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