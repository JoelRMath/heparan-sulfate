package heparansulfate;

/**
 * Represents a molecular chain that has been enzymatically or chemically cleaved.
 * <p>
 * This class stores the original sequence and the specific indices where cleavage 
 * occurred. It provides methods to extract the resulting fragments as individual 
 * sub-sequences. It is primarily used for numerical validation and cleavage 
 * simulations to check derived equations for fragment length distributions.
 */
public class CleavedSequence {

    /**
     * The integer-encoded sequence of the original parent chain before cleavage.
     */
    int[] oriseq = null;

    /**
     * The indices representing cut positions. 
     * <p>
     * <b>Convention:</b> Cuts occur on the left side of the index. For a chain of 
     * length {@code n}, the array must always include {@code 0} (start) and 
     * {@code n} (end) to correctly delineate the first and last fragments.
     */
    int[] cutPos = null;

    /**
     * Constructs a representation of a cleaved chain.
     * @param seq The original parent chain sequence.
     * @param cuts An ordered array of cut positions, including {@code 0} and {@code n}.
     */
    public CleavedSequence(int[] seq, int[] cuts) {
        oriseq = seq;
        cutPos = cuts;
    }

    /**
     * Extracts and returns the fragments resulting from the cleavage.
     * <p>
     * Each fragment is returned as an independent {@code int[]} array. The order of 
     * fragments in the resulting 2D array corresponds to their linear order 
     * (from start to end) in the original sequence.
     * * 
     * * @return A 2D array {@code int[fragment_index][sequence_data]} containing 
     * the sub-sequences of all resulting fragments.
     */
    int[][] getFragments() {
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
     * Main method for verifying the fragment extraction logic via simulation.
     * @param args Command line arguments (not used).
     */
    public static void main(String[] args) {
        int[] seq = new int[10];
        for (int i = 0; i < 10; i++)
            seq[i] = i + 1;
        
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
        
        for (int i = 0; i < cs.oriseq.length; i++)
            System.out.print(cs.oriseq[i]);
        System.out.println("");
        
        for (int i = 0; i < test.length; i++) {
            for (int j = 0; j < test[i].length; j++)
                System.out.print(test[i][j]);
            System.out.println("");
        }
    }
}