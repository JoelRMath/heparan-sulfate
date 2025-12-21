package heparansulfate;

import java.util.ArrayList;
import java.util.List;

/**
 * each instance gives infeasibility of the constraint set
 * (overall disaccharide composition and heparinase digest)
 * for BKHS chains of length n
 */
public class ProfileFeasibility {
    Species sp = null;
    double[][] A = null;
    double[] b = null;
    LinEqCons lec = null;
    SimplexPhaseI sp1 = null;
    public double infeasibility = 0.;
    String[] lab = null;

    /**
     * Constructor for a single chain length n
     */
    public ProfileFeasibility(int m, int n, String[] lab, String inDir, String outDir) {
        this.lab = lab;
        sp = new Species(m, n);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        String[] csFile = {inDir + "US.hepI.txt", inDir + "US.hepIII.txt"};
        String[] fragFile = {inDir + "hepI.f.txt", inDir + "hepIII.f.txt"};
        
        // Complete LEC combines the normalization, composition, and digest constraints
        lec = sp.getCompleteLEC(bbs, csFile, fragFile);
        A = lec.A;
        b = lec.b;
        
        // Simplex Phase I finds the minimum artificial variable cost (infeasibility)
        sp1 = new SimplexPhaseI(A, b);
        infeasibility = sp1.finalCost;
        System.err.println("Length " + n + " -> Infeasibility: " + infeasibility);
    }

    /**
     * The method called by ReproduceFigures.java
     * @param inDir Input directory path
     * @param outDir Output directory path
     * @param nMin Start of chain length range
     * @param nMax End of chain length range
     */
    public static void run(String inDir, String outDir, int nMin, int nMax) {
        String[] labels = {"U", "S"};
        List<String> results = new ArrayList<>();
        results.add("n\tinfeasibility");

        for (int n = nMin; n <= nMax; n++) {
            ProfileFeasibility pf = new ProfileFeasibility(2, n, labels, inDir, outDir);
            results.add(n + "\t" + pf.infeasibility);
            
            // Save inside the loop so you can see progress in the file
            Utils.saveFile(results, outDir + "Feasibility.res");
        }
    }

    /**
     * Legacy method preserved for compatibility with existing code
     */
    public static void makeProfile(String inDir, String outDir) {
        run(inDir, outDir, 5, 20);
    }

    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        run(inDir, outDir, 5, 20);
    }
}