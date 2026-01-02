package heparansulfate;

import java.util.ArrayList;
import java.util.List;

/**
 * Feasibility problem solver (Phase I via Simplex) to test for the possibility
 * of homogeneous disaccharide composition along BKHS chains when combined
 * with heparinase digest constraints.
 */
public class HomogeneityFeasibility {
    /**
     * Enumeration of BKHS chain sequences.
     */
    Species sp = null;
    /**
     * Matrix in the constraint {@code Ap = b}.
     */
    double[][] A = null;
    /**
     * Vector in the constraint {@code Ap = b}.
     */
    double[] b = null;
    /**
     * Convenience wrapper for the linear equality constraint {@code Ap = b}.
     */
    LinEqCons lec = null;
    /**
     * Solution to the feasibility problem (Phase I) via linear programming
     * using the Simplex method.
     */
    SimplexPhaseI sp1 = null;
    /**
     * Output of the feasibility problem: a value close to 0 indicates that 
     * homogeneity is feasible under the given constraints.
     */
    double infeasibility = 0.;
    /**
     * Disaccharide labels.
     */
    String[] lab = null;

    /**
     * Constructs a feasibility problem to test the compatibility of homogeneous 
     * disaccharide composition with heparinase digest constraints.
     * @param m Number of disaccharides.
     * @param n BKHD chain length.
     * @param lab Disaccharide labels.
     * @param inDir Input directory path (must end with "/").
     */
    public HomogeneityFeasibility(int m, int n, String[] lab, String inDir) {
        this.lab = lab;
        sp = new Species(m, n);
        BBSet bbs = new BBSet(inDir + "US.ab.txt");
        String[] csFile = new String[2];
        String[] fragFile = new String[2];
        csFile[0] = inDir + "US.hepI.txt";
        csFile[1] = inDir + "US.hepIII.txt";
        fragFile[0] = inDir + "hepI.f.txt";
        fragFile[1] = inDir + "hepIII.f.txt";
        
        List<LinEqCons> v = new ArrayList<>();
        v.add(sp.getNormLEC());
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(fragFile[0]),
                new CSpec(csFile[0], bbs), bbs)));
        v.add(LinEqCons.removeLastRow(sp.getFragLEC(Species.loadFragAbund(fragFile[1]),
                new CSpec(csFile[1], bbs), bbs)));
        v.add(sp.getHomogeneityLEC(bbs));
        
        lec = new LinEqCons(v);
        A = lec.A;
        b = lec.b;
        sp1 = new SimplexPhaseI(A, b);
        infeasibility = sp1.finalCost;
        System.err.println(n + "\t" + infeasibility);
    }

    /**
     * Computes the feasibility of combining homogeneity and heparinase digest constraints
     * across a range of BKHS chain lengths (from 5 to 20).
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
            HomogeneityFeasibility hf = new HomogeneityFeasibility(2, n, lab, inDir);
            v.add(n + "\t" + hf.infeasibility);
        }
        Utils.saveFile(v, outDir + "HomogeneityFeasibility.res");
    }

    /**
     * Main entry point for the homogeneity feasibility analysis.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/LP/";
        if (args.length >= 1){
            inDir = args[0];
        }
        if (args.length >= 2){
            outDir = args[1];
        }
        makeProfile(inDir, outDir);
    }
}