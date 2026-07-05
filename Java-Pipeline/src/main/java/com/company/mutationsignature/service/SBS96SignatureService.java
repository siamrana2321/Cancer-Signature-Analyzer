package com.company.mutationsignature.service;

import com.company.mutationsignature.model.SBS96Response;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.springframework.stereotype.Service;

import java.io.*;
import java.util.*;

@Service
public class SBS96SignatureService {

    private final Map<String, String> etiologyMap = new HashMap<>();

    private static final String GENOME_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\genome\\GRCh38.fa";
    private static final String COSMIC_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\cosmic\\COSMIC_v3.6_SBS_GRCh38.txt";
    private static final String SBS_Etiology_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\cosmic\\SBS_Etiology.csv";

    private static final Map<Character, Character> COMPLEMENT = new LinkedHashMap<>();
    static {
        COMPLEMENT.put('A', 'T');
        COMPLEMENT.put('T', 'A');
        COMPLEMENT.put('C', 'G');
        COMPLEMENT.put('G', 'C');
    }

    public SBS96SignatureService() throws Exception {
        loadEtiologyTable();
        System.out.println("SBS Etiology Loaded Successfully");
    }

    // ---------- Public analysis method ----------
    public SBS96Response analyze(String filePath) throws Exception {
        List<Mutation> mutations = parseVcf(filePath);
        if (mutations.isEmpty()) {
            throw new RuntimeException("No SNVs found in the uploaded VCF file.");
        }

        Map<String, Integer> sbs96 = buildSbs96(mutations);
        double[] tumorVector = buildNormalizedTumorVector(sbs96);
        CosmicMatrix cosmic = loadCosmicMatrix();
        NnlsResult fit = fitNnls(tumorVector, cosmic);
        List<SBS96Response.SignatureResult> results = buildNnlsResults(fit, cosmic);

        SBS96Response response = new SBS96Response();
        response.setReconstructionCosine(fit.getCosine());
        response.setPearson(fit.getPearson());
        response.setCosmicVersion("COSMIC v3.6 SBS GRCh38");
        response.setRmse(fit.getRmse());
        response.setTotalMutations(mutations.size());
        response.setSbs96Spectrum(sbs96);
        response.setSbs96Grouped(groupSbs96(sbs96));
        response.setSbs96Percentage(buildSbs96Percentages(sbs96));
        response.setCosmicContributions(results);
        response.setBiologicalAnnotations(buildAnnotations(results));
        response.setClinicalSummary(buildClinicalSummary(results));

        System.out.println("SBS96 Total = " + sbs96.values().stream().mapToInt(Integer::intValue).sum());
        System.out.println("Cosine = " + fit.getCosine());
        System.out.println("Pearson = " + fit.getPearson());
        System.out.println("RMSE = " + fit.getRmse());

        return response;
    }

    // ---------- Helper methods ----------
    private void loadEtiologyTable() throws Exception {
        File file = new File(SBS_Etiology_PATH);
        if (!file.exists()) {
            throw new RuntimeException("SBS etiology file not found.");
        }
        try (BufferedReader br = new BufferedReader(new FileReader(SBS_Etiology_PATH))) {
            String line;
            br.readLine(); // skip header
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(",", 2);
                if (parts.length >= 2) {
                    etiologyMap.put(parts[0].trim(), parts[1].trim());
                }
            }
        }
        System.out.println("Loaded SBS Etiology Entries = " + etiologyMap.size());
    }

    private Map<String, Map<String, Double>> buildSbs96Percentages(Map<String, Integer> sbs96) {
        double total = sbs96.values().stream().mapToInt(Integer::intValue).sum();
        if (total == 0) return new LinkedHashMap<>();

        Map<String, Map<String, Double>> grouped = new LinkedHashMap<>();
        for (Map.Entry<String, Integer> entry : sbs96.entrySet()) {
            String key = entry.getKey();
            String mutation = key.substring(key.indexOf('[') + 1, key.indexOf(']'));
            char left = key.charAt(0);
            char right = key.charAt(key.length() - 1);
            char ref = mutation.charAt(0);
            String context = "" + left + ref + right;
            double percentage = (entry.getValue() / total) * 100.0;
            grouped.computeIfAbsent(mutation, k -> new LinkedHashMap<>()).put(context, percentage);
        }
        return grouped;
    }

    private Map<String, Map<String, Integer>> groupSbs96(Map<String, Integer> sbs96) {
        Map<String, Map<String, Integer>> grouped = new LinkedHashMap<>();
        for (Map.Entry<String, Integer> entry : sbs96.entrySet()) {
            String key = entry.getKey();
            String mutation = key.substring(key.indexOf('[') + 1, key.indexOf(']'));
            char left = key.charAt(0);
            char right = key.charAt(key.length() - 1);
            char ref = mutation.charAt(0);
            String context = "" + left + ref + right;
            grouped.computeIfAbsent(mutation, k -> new LinkedHashMap<>()).put(context, entry.getValue());
        }
        return grouped;
    }

    private double calculatePearson(double[] a, double[] b) {
        return new PearsonsCorrelation().correlation(a, b);
    }

    private double[] buildNormalizedTumorVector(Map<String, Integer> sbs96) {
        double total = sbs96.values().stream().mapToInt(Integer::intValue).sum();
        if (total == 0.0) throw new RuntimeException("Empty SBS96 spectrum.");
        double[] vector = new double[96];
        int idx = 0;
        for (Integer value : sbs96.values()) {
            vector[idx++] = value / total;
        }
        return vector;
    }

    // ---------- NNLS fitting ----------
    private NnlsResult fitNnls(double[] tumor, CosmicMatrix cosmic) {
        int signatures = cosmic.getSignatureNames().size();
        double[][] A = cosmic.getMatrix();
        double[] weights = new double[signatures];

        for (int iter = 0; iter < 500; iter++) {
            boolean changed = false;
            for (int j = 0; j < signatures; j++) {
                double numerator = 0.0, denominator = 0.0;
                for (int row = 0; row < 96; row++) {
                    double prediction = 0.0;
                    for (int k = 0; k < signatures; k++) {
                        prediction += A[row][k] * weights[k];
                    }
                    numerator += A[row][j] * (tumor[row] - prediction + A[row][j] * weights[j]);
                    denominator += A[row][j] * A[row][j];
                }
                if (denominator == 0.0) continue;
                double updated = Math.max(0.0, numerator / denominator);
                if (Math.abs(updated - weights[j]) > 1e-10) changed = true;
                weights[j] = updated;
            }
            if (!changed) break;
        }

        boolean allZero = true;
        for (double w : weights) if (w > 0.0) { allZero = false; break; }
        if (allZero) throw new RuntimeException("Unable to fit COSMIC signatures.");

        double[] reconstructed = new double[96];
        for (int row = 0; row < 96; row++) {
            for (int col = 0; col < signatures; col++) {
                reconstructed[row] += A[row][col] * weights[col];
            }
        }

        NnlsResult result = new NnlsResult();
        result.setWeights(weights);
        result.setReconstructed(reconstructed);
        result.setCosine(cosineSimilarity(tumor, reconstructed));
        result.setPearson(calculatePearson(normalize(tumor), normalize(reconstructed)));
        result.setRmse(calculateRmse(normalize(tumor), normalize(reconstructed)));
        return result;
    }

    private double[] normalize(double[] vector) {
        double sum = 0.0;
        for (double v : vector) sum += v;
        if (sum == 0) throw new RuntimeException("Cannot normalize zero vector.");
        double[] norm = new double[vector.length];
        for (int i = 0; i < vector.length; i++) norm[i] = vector[i] / sum;
        return norm;
    }

    private double cosineSimilarity(double[] a, double[] b) {
        double dot = 0.0, normA = 0.0, normB = 0.0;
        for (int i = 0; i < a.length; i++) {
            dot += a[i] * b[i];
            normA += a[i] * a[i];
            normB += b[i] * b[i];
        }
        if (normA == 0 || normB == 0) return 0.0;
        return dot / (Math.sqrt(normA) * Math.sqrt(normB));
    }

    private double calculateRmse(double[] original, double[] reconstructed) {
        double sum = 0.0;
        for (int i = 0; i < original.length; i++) {
            double diff = original[i] - reconstructed[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum / original.length);
    }

    private List<SBS96Response.SignatureResult> buildNnlsResults(NnlsResult fit, CosmicMatrix cosmic) {
        List<SBS96Response.SignatureResult> results = new ArrayList<>();
        double total = 0.0;
        for (double w : fit.getWeights()) total += w;
        if (total == 0) return Collections.emptyList();

        for (int i = 0; i < fit.getWeights().length; i++) {
            double contribution = fit.getWeights()[i];
            if (contribution <= 0.0) continue;
            double percentage = (contribution / total) * 100.0;
            if (percentage < 1.0) continue;

            SBS96Response.SignatureResult result = new SBS96Response.SignatureResult();
            result.setSignature(cosmic.getSignatureNames().get(i));
            result.setSimilarity(percentage);
            result.setEtiology(getSignatureMeaning(cosmic.getSignatureNames().get(i)));
            results.add(result);
        }

        results.sort((a, b) -> Double.compare(b.getSimilarity(), a.getSimilarity()));
        return results.size() > 3 ? new ArrayList<>(results.subList(0, 3)) : results;
    }

    private String getSignatureMeaning(String signature) {
        return etiologyMap.getOrDefault(signature, "Unknown");
    }

    private List<SBS96Response.BiologicalAnnotation> buildAnnotations(List<SBS96Response.SignatureResult> signatures) {
        List<SBS96Response.BiologicalAnnotation> annotations = new ArrayList<>();
        for (SBS96Response.SignatureResult result : signatures) {
            SBS96Response.BiologicalAnnotation ann = new SBS96Response.BiologicalAnnotation();
            ann.setSignature(result.getSignature());
            ann.setEtiology(getSignatureMeaning(result.getSignature()));
            double sim = result.getSimilarity();
            if (sim >= 20) ann.setConfidence("High");
            else if (sim >= 5) ann.setConfidence("Moderate");
            else ann.setConfidence("Low");
            annotations.add(ann);
        }
        return annotations;
    }

    private String buildClinicalSummary(List<SBS96Response.SignatureResult> signatures) {
        StringBuilder sb = new StringBuilder();
        for (SBS96Response.SignatureResult result : signatures) {
            if (result.getSimilarity() < 5.0) continue;
            sb.append(result.getSignature())
                    .append(" (")
                    .append(String.format("%.1f", result.getSimilarity()))
                    .append("%) : ")
                    .append(result.getEtiology())
                    .append(". ");
        }
        return sb.toString();
    }

    // ---------- COSMIC matrix loading ----------
    private CosmicMatrix loadCosmicMatrix() throws Exception {
        List<String> contexts = new ArrayList<>();
        File file = new File(COSMIC_PATH);
        if (!file.exists()) throw new RuntimeException("COSMIC matrix file not found.");

        try (BufferedReader reader = new BufferedReader(new FileReader(COSMIC_PATH))) {
            String header = reader.readLine();
            String[] columns = header.split("\t");
            List<String> signatureNames = new ArrayList<>();
            for (int i = 1; i < columns.length; i++) signatureNames.add(columns[i]);

            int signatureCount = signatureNames.size();
            double[][] matrix = new double[96][signatureCount];
            String line;
            int row = 0;
            while ((line = reader.readLine()) != null && row < 96) {
                String[] tokens = line.split("\t");
                contexts.add(tokens[0]);
                for (int col = 1; col < tokens.length; col++) {
                    matrix[row][col - 1] = Double.parseDouble(tokens[col]);
                }
                row++;
            }
            return new CosmicMatrix(matrix, signatureNames, contexts);
        }
    }

    // ---------- VCF parsing ----------
    private List<Mutation> parseVcf(String filePath) throws Exception {
        List<Mutation> mutations = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] cols = line.split("\t");
                if (cols.length < 5) continue;
                String chromosome = cols[0];
                long position = Long.parseLong(cols[1]);
                String ref = cols[3];
                String alt = cols[4];
                if (ref.length() == 1 && alt.length() == 1) {
                    mutations.add(new Mutation(chromosome, position, ref, alt));
                }
            }
        }
        return mutations;
    }

    // ---------- SBS96 construction ----------
    private Map<String, Integer> buildSbs96(List<Mutation> mutations) throws Exception {
        Map<String, Integer> sbs96 = initializeSbs96();
        try (IndexedFastaSequenceFile genome = new IndexedFastaSequenceFile(new File(GENOME_PATH))) {
            for (Mutation mutation : mutations) {
                String trinucleotide;
                try {
                    ReferenceSequence seq = genome.getSubsequenceAt(
                            mutation.getChromosome(),
                            mutation.getPosition() - 1,
                            mutation.getPosition() + 1
                    );
                    trinucleotide = new String(seq.getBases()).toUpperCase();
                    if (trinucleotide.length() != 3) continue;
                } catch (Exception e) {
                    continue;
                }
                char left = trinucleotide.charAt(0);
                char ref = mutation.getRef().toUpperCase().charAt(0);
                char alt = mutation.getAlt().toUpperCase().charAt(0);
                char right = trinucleotide.charAt(2);
                String context = normalizeContext(left, ref, alt, right);
                sbs96.put(context, sbs96.getOrDefault(context, 0) + 1);
            }
        }
        return sbs96;
    }

    private Map<String, Integer> initializeSbs96() {
        Map<String, Integer> sbs96 = new LinkedHashMap<>();
        char[] bases = {'A', 'C', 'G', 'T'};
        String[] substitutions = {"C>A", "C>G", "C>T", "T>A", "T>C", "T>G"};
        for (char left : bases) {
            for (String sub : substitutions) {
                for (char right : bases) {
                    sbs96.put(left + "[" + sub + "]" + right, 0);
                }
            }
        }
        return sbs96;
    }

    private String normalizeContext(char left, char ref, char alt, char right) {
        if (ref == 'C' || ref == 'T') {
            return left + "[" + ref + ">" + alt + "]" + right;
        }
        String rc = reverseComplement("" + left + ref + right);
        char nLeft = rc.charAt(0);
        char nRight = rc.charAt(2);
        char nRef = COMPLEMENT.get(ref);
        char nAlt = COMPLEMENT.get(alt);
        return nLeft + "[" + nRef + ">" + nAlt + "]" + nRight;
    }

    private String reverseComplement(String seq) {
        StringBuilder sb = new StringBuilder();
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(COMPLEMENT.get(seq.charAt(i)));
        }
        return sb.toString();
    }

    // ---------- Inner helper classes ----------
    private static class Mutation {
        private final String chromosome;
        private final long position;
        private final String ref;
        private final String alt;

        Mutation(String chromosome, long position, String ref, String alt) {
            this.chromosome = chromosome;
            this.position = position;
            this.ref = ref;
            this.alt = alt;
        }
        String getChromosome() { return chromosome; }
        long getPosition() { return position; }
        String getRef() { return ref; }
        String getAlt() { return alt; }
    }

    private static class CosmicMatrix {
        private final double[][] matrix;
        private final List<String> signatureNames;
        private final List<String> contexts;

        CosmicMatrix(double[][] matrix, List<String> signatureNames, List<String> contexts) {
            this.matrix = matrix;
            this.signatureNames = signatureNames;
            this.contexts = contexts;
        }
        double[][] getMatrix() { return matrix; }
        List<String> getSignatureNames() { return signatureNames; }
        List<String> getContexts() { return contexts; }
    }

    private static class NnlsResult {
        private double[] weights;
        private double[] reconstructed;
        private double cosine;
        private double rmse;
        private double pearson;

        double[] getWeights() { return weights; }
        void setWeights(double[] weights) { this.weights = weights; }
        double[] getReconstructed() { return reconstructed; }
        void setReconstructed(double[] reconstructed) { this.reconstructed = reconstructed; }
        double getCosine() { return cosine; }
        void setCosine(double cosine) { this.cosine = cosine; }
        double getRmse() { return rmse; }
        void setRmse(double rmse) { this.rmse = rmse; }
        double getPearson() { return pearson; }
        void setPearson(double pearson) { this.pearson = pearson; }
    }
}