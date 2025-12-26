package heparansulfate;

import java.util.ArrayList;
import java.util.List;

/**
 * each instance gives infeasibility of the constraint set
 * (overall disaccharide composition and heparinase digest)
 * for BKHS chains of length n
 */
public class ProfileFeasibility {
    /**
     * enumerates molecular species (all possible sequences)
     */
    Species sp = null;
    /**
     * matrix for constraints
     */
    double[][] A = null;
    /**
     * vector for constraints
     */
    double[] b = null;
    /**
     * wrapper for all constraints in equality form
     */
    LinEqCons lec = null;
    /**
     * phase I of the simplex
     */
    SimplexPhaseI sp1 = null;
    /**
     * infeasibility (produced by this.sp1)
     */
    double infeasibility = 0.;
    /**
     * disaccaride labels
     */
    String[] lab = null;

    /**
     * gives infeasibility of the constraint set (overall disaccharide composition and
     * heparinase digest) for BKHS chains of length n
     * @param m number of disaccharides
     * @param n BKHS chain length
     * @param lab disaccharide labels
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
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
     * makes profile of infeasibility for BKHS chain length from 5 to 20
     * @param inDir input directory (ends with "/")
     * @param outDir output directory (ends with "/")
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
     * Main entry point
     * @param args command line arguments
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        makeProfile(inDir, outDir);
    }
}