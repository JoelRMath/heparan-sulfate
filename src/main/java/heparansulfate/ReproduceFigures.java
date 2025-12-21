package heparansulfate;

/**
 * Master class to generate all computation outputs for paper reproduction.
 */
public class ReproduceFigures {

    public static void main(String[] args) {
        // Paths relative to the project root (where the pom.xml is)
        String inDir = "input/";
        String outDir = "output/";

        System.out.println("=== Starting HS Sequencing Reproduction ===");

        // --- FIGURE 2: Feasibility Profile ---
        // Range: n = 5 to 20 as per Pradines et al. (2016)
        System.out.println("Step 1: Generating Feasibility Profile (Figure 2)...");
        ProfileFeasibility.run(inDir, outDir, 5, 20);

        // --- FIGURE 3: Model Fits (H&I and H&C) ---
        // Requires n=16 (determined from feasibility)
        System.out.println("Step 2: Generating Model Fits (Figure 3)...");
        // We call the existing static methods in your model classes
        HIModel.hepI(inDir, outDir);
        HIModel.hepIII(inDir, outDir);
        // Note: HCModelSA.run() would be called here once you have it ready

        System.out.println("=== Data Generation Complete. Check /output folder ===");
    }
}