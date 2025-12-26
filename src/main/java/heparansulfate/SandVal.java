package heparansulfate;

import java.util.Comparator;

/**
 * convenience class to sort pairs
 */
public class SandVal implements Comparator<SandVal> {
    /**
     * String field
     */
    public String s = null;
    /**
     * double field, used for sorting
     */
    double val = 0.;
    /**
     * 1 for sorting by ascending values and -1 by descending values
     */
    int sign = 1;

    /**
     * convenience class to sort two arrays, based on values in one,
     * while preserving pairs
     * @param s String (label)
     * @param val value to sort by
     * @param dir direction (ascending "a" or descending "d")
     */
    public SandVal(String s, double val, String dir) {
        this.s = s;
        this.val = val;
        if (dir.toLowerCase().startsWith("a"))
            sign = 1;
        else
            sign = -1;
    }

    /**
     * defines order
     * @param sv1 the first object to be compared.
     * @param sv2 the second object to be compared.
     * @return a negative integer, zero, or a positive integer as the first argument is less than, equal to, or greater than the second.
     */
    @Override
    public int compare(SandVal sv1, SandVal sv2) {
        if (sv1.val == sv2.val) {
            return 0;
        } else {
            if (sv1.val < sv2.val) {
                return -1 * sign;
            } else {
                return sign;
            }
        }
    }
}