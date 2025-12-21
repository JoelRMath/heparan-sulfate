package heparansulfate;

/**
 * Represents a chain that was cut (cleaved). This includes chain sequence,
 * the cut positions and methods to return the fragments. This class is
 * utilized for numerical check (cleavage simulations) of derived equations
 * for distributions of fragment lengths.
 */
public class CleavedSequence {
    /**
     * Sequence of the original chain.
     */
    int[] oriseq = null;
    /**
     * Positions of the cuts. By convention, cuts happen on the left side
     * and for a chain of length n there is always a cut at 0 and at n.
     */
    int[] cutPos = null;

    /**
     * Creates a new cut sequence of a chain.
     * * @param seq  chain sequence
     * @param cuts Positions of the cut
     */
    public CleavedSequence(int[] seq, int[] cuts) {
        oriseq = seq;
        cutPos = cuts;
    }

    /**
     * Returns the fragments as an array of sequences,
     * each sequence being represented by int[].
     * * @return fragments as an int[][]
     */
    public int[][] getFragments() {
        int[][] res = new int[cutPos.length - 1][];
        for (int i = 0; i < cutPos.length - 1; i++) {
            int l = cutPos[i + 1] - cutPos[i];
            res[i] = new int[l];
            int pos = -1;
            for (int j = cutPos[i]; j < cutPos[i + 1]; j++) {
                pos++;
                res[i][pos] = oriseq[j];
            }
        }
        return res;
    }

    /**
     * For testing
     * * @param args
     */
    public static void main(String[] args) {
        int[] seq = new int[10];
        for (int i = 0; i < 10; i++)
            seq[i] = i + 1;

        // Example: cuts at every position except the last jump which is 8 to 10
        int[] cuts = new int[10];
        cuts[0] = 0;
        cuts[1] = 1;
        cuts[2] = 2;
        cuts[3] = 3;
        cuts[4] = 4;
        cuts[5] = 5;
        cuts[6] = 6;
        cuts[7] = 7;
        cuts[8] = 8;
        cuts[9] = 10;

        CleavedSequence cs = new CleavedSequence(seq, cuts);
        int[][] test = cs.getFragments();

        System.out.println("Original Sequence:");
        for (int i = 0; i < cs.oriseq.length; i++)
            System.out.print(cs.oriseq[i] + " ");
        
        System.out.println("\nFragments:");
        for (int i = 0; i < test.length; i++) {
            for (int j = 0; j < test[i].length; j++) {
                System.out.print(test[i][j]);
            }
            System.out.println("");
        }
    }
}