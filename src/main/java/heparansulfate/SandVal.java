package heparansulfate;

import java.util.Comparator;

/**
 * Convenience class to sort pairs of Strings and doubles.
 */
public class SandVal implements Comparator<SandVal> {
    /**
     * String field (usually a sequence label)
     */
    public String s = null;
    /**
     * double field, used for sorting
     */
    public double val = 0.;
    /**
     * 1 for sorting by ascending values and -1 by descending values
     */
    int sign = 1;

    /**
     * Required for Comparator initialization in some contexts.
     */
    public SandVal() {}

    /**
     * Convenience class to sort two arrays, based on values in one,
     * while preserving pairs.
     * @param s String (label)
     * @param val value to sort by
     * @param dir direction (ascending "a" or descending "d")
     */
    public SandVal(String s, double val, String dir) {
        this.s = s;
        this.val = val;
        if (dir.toLowerCase().startsWith("a")) {
            sign = 1;
        } else {
            sign = -1;
        }
    }

    /**
     * Defines order for sorting.
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
