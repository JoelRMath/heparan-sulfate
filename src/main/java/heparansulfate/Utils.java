package heparansulfate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;

/**
 * general methods (mostly I/O) used by several classes
 */
public class Utils {

    /**
     * Default constructor
     */
    public Utils() {
    }

    /**
     * random selection of index i in F based on probabilities F[i+1]-F[i]
     * @param F cumulative probabilities
     * @param rand random number generator
     * @return random index i in F based on probabilities F[i+1]-F[i]
     */
    public static int getRandIndexInF(double[] F, Random rand) {
        double d = rand.nextDouble() * 0.99999999;
        int res = 0;
        if (d <= F[res]) {
            return res;
        } else {
            while (d > F[res]) {
                res++;
            }
            return res;
        }
    }

    /**
     * minimum
     * @param x array of values
     * @return min of x
     */
    public static double getMin(double[] x) {
        double res = x[0];
        for (int i = 0; i < x.length; i++) {
            if (x[i] < res) {
                res = x[i];
            }
        }
        return res;
    }

    /**
     * maximum
     * @param x array of values
     * @return max of x
     */
    public static double getMax(double[] x) {
        double res = x[0];
        for (int i = 0; i < x.length; i++) {
            if (x[i] > res) {
                res = x[i];
            }
        }
        return res;
    }

    /**
     * average
     * @param x array of values
     * @return average of x
     */
    public static double getAverage(double[] x) {
        double res = 0.;
        for (int i = 0; i < x.length; i++) {
            res += x[i];
        }
        res /= (double) x.length;
        return res;
    }

    /**
     * standard deviation
     * @param x array of values
     * @return standard deviation of x
     */
    public static double getSD(double[] x) {
        double res = 0.;
        double av = getAverage(x);
        for (int i = 0; i < x.length; i++) {
            res += (x[i] - av) * (x[i] - av);
        }
        res /= (double)(x.length - 1);
        res = Math.sqrt(res);
        return res;
    }

    /**
     * loads an ASCII file into a vector, line by line
     * @param file ASCII file, tab-delimited
     * @return Vector representation (one line in each element) of file
     */
    public static List<String> loadFileNoheader(String file) {
        List<String> res = new ArrayList<>();
        try {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                StringTokenizer st = new StringTokenizer(line, "\t");
                if (st.countTokens() >= 2) {
                    res.add(line);
                }
            }
            br.close();
            fr.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return res;
    }

    /**
     * saves each element of Vector v into a line of ASCII file ’file’
     * @param v contains file lines, each one being tab-delimited
     * @param file output file
     */
    public static void saveFile(List<String> v, String file) {
        try {
            FileWriter fw = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < v.size(); i++) {
                bw.write(v.get(i));
                bw.newLine();
            }
            bw.close();
            fw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}