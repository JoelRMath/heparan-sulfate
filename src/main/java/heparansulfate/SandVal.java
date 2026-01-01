package heparansulfate;

import java.util.Comparator;

/**
 * Convenience class to sort pairs of strings and double values.
 */
public class SandVal implements Comparator<SandVal> {
    /**
     * String field representing the label.
     */
    public String s = null;
    /**
     * Double field used as the sorting criteria.
     */
    double val = 0.;
    /**
     * Sorting direction: {@code 1} for ascending and {@code -1} for descending.
     */
    int sign = 1;

    /**
     * Constructs a pair containing a string label and a numerical value for sorting.
     * @param s String label.
     * @param val Value to sort by.
     * @param dir Direction: "a" for ascending or "d" for descending.
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
     * Compares two {@code SandVal} objects based on their numerical values and the specified direction.
     * @param sv1 The first object to be compared.
     * @param sv2 The second object to be compared.
     * @return A negative integer, zero, or a positive integer as the first argument 
     * is less than, equal to, or greater than the second.
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