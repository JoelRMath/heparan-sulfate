package heparansulfate;

import java.util.ArrayList;
import java.util.List;

/**
 * Calculates the infeasibility of a specific constraint set (overall disaccharide composition 
 * and heparinase digest) for BKHS chains of a given length {@code n}.
 */
public class ProfileFeasibility {
    /**
     * Enumerates molecular species (all possible sequences for the given length).
     */
    Species sp = null;
    /**
     * Matrix representing the linear equality constraints {@code Ax = b}.
     */
    double[][] A = null;
    /**
     * Constant vector for the linear equality constraints {@code Ax = b}.
     */
    double[] b = null;
    /**
     * Convenience wrapper containing all equality constraints.
     */
    LinEqCons lec = null;
    /**
     * Phase I of the Simplex algorithm used to determine feasibility.
     */
    SimplexPhaseI sp1 = null;
    /**
     * Calculated infeasibility value (residual cost from Phase I).
     */
    double infeasibility = 0.;
    /**
     * Labels for the disaccharides used in the model.
     */
    String[] lab = null;

    /**
     * Constructs an instance to determine the infeasibility of the constraint set for BKHS chains of length {@code n}.
     * 
     * @param m Number of disaccharides.
     * @param n BKHS chain length.
     * @param lab Disaccharide labels.
     * @param inDir Input directory path (must end with "/").
     * @param outDir Output directory path (must end with "/").
     */
    public ProfileFeasibility(int m, int n, String[] lab, String inDir, String outDir) {
        this.lab = lab;
        sp = new Species(m, n);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        String[] csFile = new String[2];
        String[] fragFile = new String[2];
        csFile[0] = inDir + "US.hepI.txt";
        csFile[1] = inDir + "US.hepIII.txt";
        fragFile[0] = inDir + "hepI.f.txt";
        fragFile[1] = inDir + "hepIII.f.txt";
        lec = sp.getCompleteLEC(bbs, csFile, fragFile);
        A = lec.A;
        b = lec.b;
        sp1 = new SimplexPhaseI(A, b);
        infeasibility = sp1.finalCost;
        System.err.println(n + "\t" + infeasibility);
    }

    /**
     * Generates a feasibility profile for BKHS chain lengths ranging from 5 to 20.
     * @param inDir Input directory path (must end with "/").
     * @param outDir Output directory path (must end with "/").
     */
    public static void makeProfile(String inDir, String outDir) {
        String[] lab = new String[2];
        lab[0] = "U";
        lab[1] = "S";
        List<String> v = new ArrayList<>();
        v.add("n\tinfeasibility");
        for (int n = 5; n <= 20; n++) {
            ProfileFeasibility hf = new ProfileFeasibility(2, n, lab, inDir, outDir);
            v.add(n + "\t" + hf.infeasibility);
        }
        Utils.saveFile(v, outDir + "Feasibility.res");
    }

    /**
     * Main entry point for generating the feasibility profile.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/LP/";
        makeProfile(inDir, outDir);
    }
}