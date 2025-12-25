package heparansulfate;

import java.util.ArrayList;
import java.util.List;

/**
 * Feasibility problem (phase I via simplex) to test for the possibility
 * of homogeneous disaccharide composition along BKHS chains when combined
 * with heparinase digest constraints
 */
public class HomogeneityFeasibility {
    /**
     * enumeration of BKHS chain sequences
     */
    Species sp = null;
    /**
     * matrix in constraing Ap = b
     */
    double[][] A = null;
    /**
     * vector in constraint Ap = b
     */
    double[] b = null;
    /**
     * convenience wrapper for constraint Ap = b
     */
    LinEqCons lec = null;
    /**
     * solution to feasibility problem (phase I) via linear programming
     * (simplex method)
     */
    SimplexPhaseI sp1 = null;
    /**
     * output of the feasibility problem: close to 0 if homogeneity is feasible
     */
    double infeasibility = 0.;
    /**
     * disaccharide labels
     */
    String[] lab = null;

    /**
     * Feasibility problem (phase I via simplex) to test for the possibility
     * of homogeneous disaccharide composition along BKHS chains when combined with
     * heparinase digest constraints
     * @param m number of disaccharides
     * @param n BKHD chain length
     * @param lab disaccharide labels
     * @param inDir input directory (ends with "/")
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
     * computes feasibility of combining homogeneity and heparinase digest constraints
     * for BKHS chain lengths between 5 and 20
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
            HomogeneityFeasibility hf = new HomogeneityFeasibility(2, n, lab, inDir);
            v.add(n + "\t" + hf.infeasibility);
            Utils.saveFile(v, outDir + "HomogeneityFeasibility.res");
        }
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/";
        makeProfile(inDir, outDir);
    }
}