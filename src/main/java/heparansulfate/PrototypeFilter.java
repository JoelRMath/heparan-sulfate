package heparansulfate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This class is a post-publication class implemented in order prototype reduction of the size of species abundance files (from 126MB to less than 20MB)
 */
public class PrototypeFilter {

    /**
     * Inner class to hold line data for sorting
     */
    static class Entry implements Comparable<Entry> {
        String line;
        double probability;

        public Entry(String line, double probability) {
            this.line = line;
            this.probability = probability;
        }

        // Sort descending (largest probability first)
        @Override
        public int compareTo(Entry o) {
            return Double.compare(o.probability, this.probability);
        }
    }

    public static void analyzeReduction(String filePath, double targetMass) {
        System.out.println("Reading file: " + filePath + " ...");
        List<Entry> entries = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            int lineCount = 0;
            
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                lineCount++;

                // Split by tab
                String[] parts = line.split("\t");
                
                // Format: First column is abundance
                try {
                    double p = Double.parseDouble(parts[0]);
                    entries.add(new Entry(line, p));
                } catch (NumberFormatException | ArrayIndexOutOfBoundsException e) {
                    System.err.println("Skipping malformed line " + lineCount + ": " + line.substring(0, Math.min(line.length(), 20)) + "...");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        System.out.println("Sorting " + entries.size() + " entries...");
        
        // 1. Sort Descending
        Collections.sort(entries);

        // 2. Accumulate Mass
        double currentMass = 0.0;
        int linesKept = 0;
        int totalLines = entries.size();

        for (Entry e : entries) {
            currentMass += e.probability;
            linesKept++;
            // Stop once we cross the threshold
            if (currentMass >= targetMass) {
                break;
            }
        }

        // 3. Report Stats
        double percentKept = (double) linesKept / totalLines * 100.0;
        
        System.out.println("------------------------------------------------");
        System.out.println("File Analyzed               : " + filePath);
        System.out.println("Target Cumulative Abundance : " + targetMass);
        System.out.println("Total Species (Lines)       : " + totalLines);
        System.out.println("Species Kept                : " + linesKept);
        System.out.println("Reduction                   : " + String.format("%.2f", (100 - percentKept)) + "% of lines removed");
        System.out.println("Actual Cumulative Mass      : " + currentMass);
        System.out.println("------------------------------------------------");
    }

    public static void main(String[] args) {
        // Adjust filenames to match your output folder structure
        String file1 = "output/MEW/MEMW.pind.res";
        // String file2 = "output/MEW/MEMW.Sig15.res"; 

        // Analyze standard model
        analyzeReduction(file1, 0.999);
        
        System.out.println();
        
        // Analyze Sigma 1.5 model (if it exists)
        // analyzeReduction(file2, 0.9999);
    }
}