package heparansulfate;

/**
 * Master class to generate all computation outputs for reproduction of figures.
 */
public class ReproduceFigures {

    public static void main(String[] args) {
        String inDir = "input/";
        String outDir = "output/HIM/";

        System.out.println("Figuer 1: Generating results for HIM...");
        HIModel.hepI(inDir, outDir); 
        HIModel.hepIII(inDir, outDir);
        System.out.println("... done");

    }
}