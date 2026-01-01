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
 * General utility methods, primarily handling I/O operations and basic statistics, 
 * used across several classes in the package.
 */
public class Utils {

    /**
     * Default constructor.
     */
    public Utils() {
    }

    /**
     * Performs a random selection of an index based on cumulative probabilities.
     * @param F Array of cumulative probabilities.
     * @param rand Random number generator.
     * @return A random index {@code i} such that the selection is based on the 
     * interval probabilities defined by {@code F}.
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
     * Finds the minimum value in an array.
     * @param x Array of double values.
     * @return The minimum value in {@code x}.
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
     * Finds the maximum value in an array.
     * @param x Array of double values.
     * @return The maximum value in {@code x}.
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
     * Calculates the arithmetic average of an array.
     * @param x Array of double values.
     * @return The average of the values in {@code x}.
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
     * Calculates the sample standard deviation of an array.
     * @param x Array of double values.
     * @return The standard deviation of the values in {@code x}.
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
     * Loads a tab-delimited ASCII file into a list of strings, skipping the header.
     * @param file Path to the tab-delimited ASCII file.
     * @return A {@code List} of strings, where each element represents a line from the file.
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
     * Saves each element of a list into a line of an ASCII file.
     * @param v List containing the lines to be saved.
     * @param file Path to the output file.
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