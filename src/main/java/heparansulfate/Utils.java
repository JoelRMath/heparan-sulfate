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
     * random selection of index i in F based on probabilities F[i+1]-F[i]
     * @param F cumulative probabilities
     * @param rand Random number generator
     * @return random index i in F based on probabilities
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
     * minimum of array x
     */
    public static double getMin(double[] x) {
        double res = x[0];
        for (double val : x) {
            if (val < res) {
                res = val;
            }
        }
        return res;
    }

    /**
     * maximum of array x
     */
    public static double getMax(double[] x) {
        double res = x[0];
        for (double val : x) {
            if (val > res) {
                res = val;
            }
        }
        return res;
    }

    /**
     * average of array x
     */
    public static double getAverage(double[] x) {
        double res = 0.;
        for (double val : x) {
            res += val;
        }
        res /= (double) x.length;
        return res;
    }

    /**
     * standard deviation of array x
     */
    public static double getSD(double[] x) {
        double res = 0.;
        double av = getAverage(x);
        for (double val : x) {
            res += (val - av) * (val - av);
        }
        res /= (double) (x.length - 1);
        res = Math.sqrt(res);
        return res;
    }

    /**
     * loads an ASCII file into a list, line by line, skipping header
     * @param file ASCII file, tab-delimited
     * @return List representation of file lines
     */
    public static List<String> loadFileNoheader(String file) {
        List<String> res = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            // Read and discard header
            br.readLine();
            String line;
            while ((line = br.readLine()) != null) {
                StringTokenizer st = new StringTokenizer(line, "\t");
                if (st.countTokens() >= 2) {
                    res.add(line);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return res;
    }

    /**
     * saves each element of List v into a line of ASCII file
     * @param v contains file lines, each one being tab-delimited
     * @param file output file path
     */
    public static void saveFile(List<String> v, String file) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            for (String line : v) {
                bw.write(line);
                bw.newLine();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}