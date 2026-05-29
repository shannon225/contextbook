package edu.washington.gs.maccoss.encyclopedia.context;

import java.awt.Color;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.swing.JComponent;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;

public class ContextMProphetPlotter {

    private static class MProphetRow {
        double score;
        double posteriorErrorProb;

        MProphetRow(double score, double posteriorErrorProb) {
            this.score = score;
            this.posteriorErrorProb = posteriorErrorProb;
        }
    }

    public static void plotContextMProphetResults(
            File backgroundTargetFile,
            File backgroundDecoyFile,
            File referenceTargetFile,
            File referenceDecoyFile,
            File outputDirectory
    ) throws Exception {

        if (!outputDirectory.exists()) {
            outputDirectory.mkdirs();
        }

        ArrayList<MProphetRow> backgroundTargets = readMProphetRows(backgroundTargetFile, true);
        ArrayList<MProphetRow> backgroundDecoys = readMProphetRows(backgroundDecoyFile, true);
        ArrayList<MProphetRow> referenceTargets = readMProphetRows(referenceTargetFile, false);
        ArrayList<MProphetRow> referenceDecoys = readMProphetRows(referenceDecoyFile, false);
       
        
        writeDensityPlot(
                backgroundTargets,
                backgroundDecoys,
                new File(outputDirectory, "background_target_decoy_score_density.pdf"),
                "Background peptides: primary score density"
        );

        writeDensityPlot(
                referenceTargets,
                referenceDecoys,
                new File(outputDirectory, "reference_target_decoy_score_density.pdf"),
                "Reference peptides: primary score density"
        );

        writeScoreVsPosteriorErrorPlot(
                referenceTargets,
                new File(outputDirectory, "reference_targets_score_vs_posterior_error_prob.pdf"),
                "Reference targets: score vs posterior error probability"
        );
    }

    private static ArrayList<MProphetRow> readMProphetRows(File file, boolean keepOnlyTargetScores) throws Exception {
        ArrayList<MProphetRow> rows = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String header = br.readLine();

            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("pi_0=") || line.trim().isEmpty()) {
                    continue;
                }

                String[] columns = line.split("\t", -1);

                if (keepOnlyTargetScores && !isTargetScoreRow(columns)) {
                    continue;
                }

                double score = Double.parseDouble(columns[1]);
                double posteriorErrorProb = Double.parseDouble(columns[3]);

                rows.add(new MProphetRow(score, posteriorErrorProb));
            }
        }

        return rows;
    }

    private static void writeDensityPlot(
            ArrayList<MProphetRow> targets,
            ArrayList<MProphetRow> decoys,
            File outputFile,
            String title
    ) {

        XYTraceInterface targetTrace = makeDensityTrace(targets, "Targets", Color.RED);
        XYTraceInterface decoyTrace = makeDensityTrace(decoys, "Decoys", Color.BLUE);

        XYTraceInterface[] traces;

        if (decoys.isEmpty()) {
            traces = new XYTraceInterface[] { targetTrace };
        } else {
            traces = new XYTraceInterface[] { targetTrace, decoyTrace };
        }

        Charter.writeAsPDF(
                outputFile,
                "Primary score",
                "Density",
                true,
                traces
        );
    }
    
    private static boolean isTargetScoreRow(String[] columns) {
        String id = columns[0].toLowerCase();

        return id.contains("target");
    }
    
    private static XYTraceInterface makeDensityTrace(
            ArrayList<MProphetRow> rows,
            String name,
            Color color
    ) {

        if (rows.isEmpty()) {
            return new XYTrace(new double[] {}, new double[] {}, GraphType.line, name, color, 2.0f);
        }

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;

        for (MProphetRow row : rows) {
            min = Math.min(min, row.score);
            max = Math.max(max, row.score);
        }

        double padding = (max - min) * 0.1;
        min -= padding;
        max += padding;

        int nPoints = 200;
        double[] x = new double[nPoints];
        double[] y = new double[nPoints];

        double bandwidth = estimateBandwidth(rows);

        for (int i = 0; i < nPoints; i++) {
            x[i] = min + i * (max - min) / (nPoints - 1);
            y[i] = gaussianKDE(rows, x[i], bandwidth);
        }

        return new XYTrace(x, y, GraphType.line, name, color, 2.0f);
    }

    private static double estimateBandwidth(ArrayList<MProphetRow> rows) {
        if (rows.size() < 2) {
            return 1.0;
        }

        double mean = 0.0;
        for (MProphetRow row : rows) {
            mean += row.score;
        }
        mean /= rows.size();

        double variance = 0.0;
        for (MProphetRow row : rows) {
            double diff = row.score - mean;
            variance += diff * diff;
        }
        variance /= (rows.size() - 1);

        double sd = Math.sqrt(variance);

        double bandwidth = 1.06 * sd * Math.pow(rows.size(), -0.2);

        if (bandwidth <= 0.0 || Double.isNaN(bandwidth)) {
            bandwidth = 1.0;
        }

        return bandwidth;
    }
    
   
    private static int getMaxCount(int[] targetCounts, int[] decoyCounts) {
        int max = 0;

        for (int count : targetCounts) {
            max = Math.max(max, count);
        }

        for (int count : decoyCounts) {
            max = Math.max(max, count);
        }

        return max;
    }

    private static double getMinScore(
            ArrayList<MProphetRow> targets,
            ArrayList<MProphetRow> decoys
    ) {
        double min = Double.POSITIVE_INFINITY;

        for (MProphetRow row : targets) {
            min = Math.min(min, row.score);
        }

        for (MProphetRow row : decoys) {
            min = Math.min(min, row.score);
        }

        return min;
    }

    private static double getMaxScore(
            ArrayList<MProphetRow> targets,
            ArrayList<MProphetRow> decoys
    ) {
        double max = Double.NEGATIVE_INFINITY;

        for (MProphetRow row : targets) {
            max = Math.max(max, row.score);
        }

        for (MProphetRow row : decoys) {
            max = Math.max(max, row.score);
        }

        return max;
    }

    private static double gaussianKDE(ArrayList<MProphetRow> rows, double x, double bandwidth) {
        double sum = 0.0;

        for (MProphetRow row : rows) {
            double z = (x - row.score) / bandwidth;
            sum += Math.exp(-0.5 * z * z);
        }

        return sum / (rows.size() * bandwidth * Math.sqrt(2.0 * Math.PI));
    }

    private static void writeScoreVsPosteriorErrorPlot(
            ArrayList<MProphetRow> referenceTargets,
            File outputFile,
            String title
    ) {

        double[] x = new double[referenceTargets.size()];
        double[] y = new double[referenceTargets.size()];

        for (int i = 0; i < referenceTargets.size(); i++) {
            x[i] = referenceTargets.get(i).score;
            y[i] = referenceTargets.get(i).posteriorErrorProb;
        }

        XYTraceInterface scatter = new XYTrace(
                x,
                y,
                GraphType.point,
                "Reference targets",
                Color.BLACK,
                2.0f
        );

        Charter.writeAsPDF(
                outputFile,
                "Primary score",
                "Posterior error probability",
                false,
                scatter
        );
    }
    
}