package com.company.mutationsignature.service;

import com.company.mutationsignature.model.ID83Response;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

@Service
public class ID83SignatureService {
    private static final Logger log = LoggerFactory.getLogger(ID83SignatureService.class);

    // Biological annotations (extend as needed)
    private static final Map<String, String> ETIOLOGY_MAP = new HashMap<>();
    static {
        ETIOLOGY_MAP.put("ID1", "Slippage during DNA replication of the replicated DNA strand");
        ETIOLOGY_MAP.put("ID2", "Slippage during DNA replication of the replicated DNAs strand");
        ETIOLOGY_MAP.put("ID3", "Tobacco smoking");
        ETIOLOGY_MAP.put("ID4", "Unknown");
        ETIOLOGY_MAP.put("ID5", "Unknown");
        ETIOLOGY_MAP.put("ID6", "Defective homologous recombination DNA damage repair");
        ETIOLOGY_MAP.put("ID7", "Defective DNA mismatch repair");
        ETIOLOGY_MAP.put("ID8", "Repair of DNA double strand breaks by NHEJ or mutations in TOP2A");
        ETIOLOGY_MAP.put("ID9", "Unknown");
        ETIOLOGY_MAP.put("ID10", "Unknown");
        ETIOLOGY_MAP.put("ID11", "Unknown");
        ETIOLOGY_MAP.put("ID12", "Unknown");
        ETIOLOGY_MAP.put("ID13", "Ultraviolet light exposure");
        ETIOLOGY_MAP.put("ID14", "Unknown");
        ETIOLOGY_MAP.put("ID15", "Unknown");
        ETIOLOGY_MAP.put("ID16", "Unknown");
        ETIOLOGY_MAP.put("ID17", "Mutations in topoisomerase TOP2A");
        ETIOLOGY_MAP.put("ID18", "Colibactin exposure (E.coli bacteria carrying pks pathogenicity island)");
        ETIOLOGY_MAP.put("ID19", "Unknown");
        ETIOLOGY_MAP.put("ID20", "Unknown");
        ETIOLOGY_MAP.put("ID21", "Unknown");
        ETIOLOGY_MAP.put("ID22", "Unknown");
        ETIOLOGY_MAP.put("ID23", "Aristolochic acid exposure");
        ETIOLOGY_MAP.put("ID24", "Unknown");
        ETIOLOGY_MAP.put("ID25", "Unknown");
    }
    private static final String GENOME_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\genome\\GRCh38.fa";
    private static final String COSMIC_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\cosmic\\COSMIC_Human_ID-83_GRCh37_v3.6.tsv";

    private static class Annotation {
        String etiology, mechanism, association;
        Annotation(String e, String m, String a) { etiology=e; mechanism=m; association=a; }
    }

    // ------------------- Main pipeline -------------------
    public ID83Response runPipeline(Path vcfPath) throws Exception {
        log.info("Starting pipeline: VCF={}", vcfPath);

        // 1. Load genome
        Map<String, String> genome = loadGenome(Paths.get(GENOME_PATH));
        log.info("Genome loaded: {} chromosomes", genome.size());

        // 2. Build ID83 matrix from VCF
        Map<String, Integer> matrix = buildMatrixFromVcf(vcfPath, genome);
        log.info("Matrix built with {} categories", matrix.size());

        // 3. Load COSMIC signatures
        CosmicData cosmic = loadCosmicSignatures(Paths.get(COSMIC_PATH));
        log.info("COSMIC loaded: {} categories, {} signatures", cosmic.categories.size(), cosmic.signatureNames.size());

        // 4. Align matrix with COSMIC order
        double[] sampleVector = alignToCosmic(matrix, cosmic.categories);
        long totalMutations = (long) Arrays.stream(sampleVector).sum();
        if (totalMutations == 0) {
            throw new IllegalArgumentException("Total mutations = 0. Check VCF and genome.");
        }

        // 5. NNLS fitting
        double[] exposures = solveNNLS(cosmic.matrix, sampleVector, 2000, 1e-8);

        // 6. Reconstruct and compute metrics
        double[] reconstructed = reconstruct(cosmic.matrix, exposures);
        double cosSim = cosineSimilarity(sampleVector, reconstructed);
        double pearson = computePearson(sampleVector, reconstructed);

        // L1, L2, KL divergence
        double l1Norm = computeL1Norm(sampleVector, reconstructed);
        double l2Norm = computeL2Norm(sampleVector, reconstructed);
        double l1NormPercent = (totalMutations > 0) ? 100.0 * l1Norm / totalMutations : 0.0;
        double l2NormOriginal = Math.sqrt(Arrays.stream(sampleVector).map(v -> v*v).sum());
        double l2NormPercent = (l2NormOriginal > 0) ? 100.0 * l2Norm / l2NormOriginal : 0.0;

        // KL divergence (normalize vectors)
        double[] p = normalize(sampleVector);
        double[] q = normalize(reconstructed);
        double klDiv = computeKLDivergence(p, q);

        // 7. Prepare exposure results, sort, take top 3
        List<ExposureResult> allExposures = new ArrayList<>();
        double totalExposure = Arrays.stream(exposures).sum();
        for (int i = 0; i < exposures.length; i++) {
            if (exposures[i] > 1e-6) {
                double percent = (totalExposure > 0) ? 100.0 * exposures[i] / totalExposure : 0.0;
                allExposures.add(new ExposureResult(cosmic.signatureNames.get(i), exposures[i], percent));
            }
        }
        allExposures.sort((a, b) -> Double.compare(b.percent, a.percent));
        List<ExposureResult> top3 = allExposures.stream().limit(3).collect(Collectors.toList());

        // 8. Build contributions, annotations, and clinical summary from top 3
        List<ID83Response.CosmicContribution> contributions = new ArrayList<>();
        List<ID83Response.BiologicalAnnotation> annotationsList = new ArrayList<>();
        StringBuilder clinicalSummary = new StringBuilder();

        for (ExposureResult er : top3) {
            String sig = er.signature;
            double percentRounded = round4(er.percent);
//            Annotation ann = ANNOTATIONS.getOrDefault(sig, new Annotation("Unknown", "Unknown", "Unknown"));
            String etiology = ETIOLOGY_MAP.getOrDefault(sig, "Unknown");

            // Confidence based on percentage threshold (as in SBS96)
            String confidence;
            if (percentRounded >= 20.0) {
                confidence = "High";
            } else if (percentRounded >= 5.0) {
                confidence = "Moderate";
            } else {
                confidence = "Low";
            }

            contributions.add(new ID83Response.CosmicContribution(sig, percentRounded, etiology));
            annotationsList.add(new ID83Response.BiologicalAnnotation(sig, etiology, confidence));

            if (clinicalSummary.length() > 0) clinicalSummary.append(". ");
            clinicalSummary.append(sig).append(" (").append(String.format("%.1f", percentRounded)).append("%) : ").append(etiology);
        }

        // 9. Build id83Percentage map with rounded values (4 decimals)
        Map<String, Double> percentageMap = new LinkedHashMap<>();
        for (Map.Entry<String, Integer> entry : matrix.entrySet()) {
            double pct = 100.0 * entry.getValue() / totalMutations;
            percentageMap.put(entry.getKey(), round4(pct));
        }
        Map<String, Map<String, Double>> grouped = buildGroupedPercentage(matrix, totalMutations);


        // 10. Build final response
        ID83Response response = new ID83Response();
        response.setCosmicVersion("COSMIC v3.6 ID83 GRCh38");
        response.setTotalMutations(totalMutations);
        response.setId83Spectrum(matrix);
        response.setId83Percentage(percentageMap);
        response.setId83Grouped(grouped);
        response.setCosmicContributions(contributions);
        response.setBiologicalAnnotations(annotationsList);
        response.setClinicalSummary(clinicalSummary.toString());
        response.setReconstructionCosine(cosSim);
        response.setPearson(pearson);

        // Set new metrics
        response.setL1Norm(l1Norm);
        response.setL1NormPercent(l1NormPercent);
        response.setL2Norm(l2Norm);
        response.setL2NormPercent(l2NormPercent);
        response.setKlDivergence(klDiv);

        log.info("Pipeline completed. Cosine = {}, Pearson = {}, L1_Norm% = {}", cosSim, pearson, l1NormPercent);
        return response;
    }

    // Helper to round to 4 decimal places
    private double round4(double value) {
        return Math.round(value * 10000.0) / 10000.0;
    }

    // Normalize vector to sum to 1
    private double[] normalize(double[] vector) {
        double sum = 0.0;
        for (double v : vector) sum += v;
        if (sum == 0) throw new IllegalArgumentException("Cannot normalize zero vector");
        double[] norm = new double[vector.length];
        for (int i = 0; i < vector.length; i++) norm[i] = vector[i] / sum;
        return norm;
    }

    // ------ Metrics computations ------
    private double computeL1Norm(double[] a, double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum += Math.abs(a[i] - b[i]);
        }
        return sum;
    }

    private double computeL2Norm(double[] a, double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    private double computeKLDivergence(double[] p, double[] q) {
        double kl = 0.0;
        for (int i = 0; i < p.length; i++) {
            if (p[i] > 0 && q[i] > 0) {
                kl += p[i] * Math.log(p[i] / q[i]);
            }
        }
        return kl;
    }

    // ------------------- Private helpers (unchanged from previous) -------------------
    private Map<String, String> loadGenome(Path fastaPath) throws Exception {
        Map<String, String> seqs = new LinkedHashMap<>();
        try (BufferedReader reader = Files.newBufferedReader(fastaPath)) {
            String line;
            String currentChrom = null;
            StringBuilder currentSeq = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith(">")) {
                    if (currentChrom != null) {
                        seqs.put(currentChrom, currentSeq.toString());
                    }
                    String[] parts = line.substring(1).trim().split("\\s+");
                    currentChrom = parts[0];
                    currentSeq = new StringBuilder();
                } else {
                    currentSeq.append(line);
                }
            }
            if (currentChrom != null) {
                seqs.put(currentChrom, currentSeq.toString());
            }
        }
        return seqs;
    }

    private Map<String, Integer> buildMatrixFromVcf(Path vcfPath, Map<String, String> genome) throws Exception {
        Map<String, Integer> matrix = initializeMatrix();
        try (BufferedReader reader = Files.newBufferedReader(vcfPath)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] cols = line.split("\t");
                if (cols.length < 5) continue;
                String chrom = cols[0];
                int pos = Integer.parseInt(cols[1]);
                String ref = cols[3];
                String alt = cols[4];
                String category = classifyVariant(genome, chrom, pos, ref, alt);
                if (category != null) {
                    matrix.merge(category, 1, Integer::sum);
                }
            }
        }
        return matrix;
    }

    private Map<String, Integer> initializeMatrix() {
        Map<String, Integer> matrix = new LinkedHashMap<>();
        for (String op : Arrays.asList("Del", "Ins")) {
            for (String base : Arrays.asList("C", "T")) {
                for (int r = 0; r <= 5; r++) {
                    matrix.put("1:" + op + ":" + base + ":" + r, 0);
                }
            }
        }
        for (int L : Arrays.asList(2, 3, 4, 5)) {
            for (String op : Arrays.asList("Del", "Ins")) {
                for (int r = 0; r <= 5; r++) {
                    matrix.put(L + ":" + op + ":R:" + r, 0);
                }
            }
        }
        matrix.put("2:Del:M:1", 0);
        matrix.put("3:Del:M:1", 0);
        matrix.put("3:Del:M:2", 0);
        matrix.put("4:Del:M:1", 0);
        matrix.put("4:Del:M:2", 0);
        matrix.put("4:Del:M:3", 0);
        for (int mh = 1; mh <= 5; mh++) {
            matrix.put("5:Del:M:" + mh, 0);
        }
        return matrix;
    }

    private String classifyVariant(Map<String, String> genome, String chrom, int pos, String ref, String alt) {
        if (!genome.containsKey(chrom)) {
            if (!chrom.startsWith("chr")) {
                chrom = "chr" + chrom;
            }
        }
        if (!genome.containsKey(chrom)) return null;
        ref = ref.toUpperCase();
        alt = alt.toUpperCase();
        String seq = genome.get(chrom);

        int lenRef = ref.length(), lenAlt = alt.length();
        String mutType;
        if (lenRef - lenAlt == lenRef - 1) mutType = "Del";
        else if (lenAlt - lenRef == lenAlt - 1) mutType = "Ins";
        else return null;

        if (mutType.equals("Del")) {
            String typeSeq = ref.substring(1);
            int typeLen = typeSeq.length();
            String sequence = typeSeq;

            int posRev = pos;
            String actualSeq = "";
            int start0 = posRev - typeLen;
            if (start0 >= 1) actualSeq = seq.substring(start0, posRev);
            while (posRev - typeLen > 0 && actualSeq.equals(typeSeq)) {
                sequence = actualSeq + sequence;
                posRev -= typeLen;
                int newStart = posRev - typeLen;
                if (newStart >= 1) actualSeq = seq.substring(newStart, posRev);
                else actualSeq = "";
            }

            int posForward = pos + typeLen;
            String newSeq = "";
            if (posForward + typeLen <= seq.length()) {
                newSeq = seq.substring(posForward, posForward + typeLen);
            }
            while (newSeq.equals(typeSeq)) {
                sequence += newSeq;
                posForward += typeLen;
                if (posForward + typeLen <= seq.length()) {
                    newSeq = seq.substring(posForward, posForward + typeLen);
                } else newSeq = "";
            }

            if (typeLen > 1 && sequence.length() == typeLen) {
                String forwardHom = ref.substring(1, ref.length() - 1);
                String reverseHom = ref.substring(2);
                int forHom = 0, revHom = 0;
                String forSeq = "", revSeq = "";

                int p = pos + typeLen;
                for (int i = forwardHom.length(); i >= 1; i--) {
                    if (p + i <= seq.length()) {
                        String sub = seq.substring(p, p + i);
                        if (sub.equals(forwardHom.substring(0, i))) {
                            forSeq = forwardHom.substring(0, i);
                            forHom = forSeq.length();
                            break;
                        }
                    }
                }
                p = pos;
                for (int i = reverseHom.length(); i >= 1; i--) {
                    if (p - i >= 0) {
                        String sub = seq.substring(p - i, p);
                        if (sub.equals(reverseHom.substring(reverseHom.length() - i))) {
                            revSeq = reverseHom.substring(reverseHom.length() - i);
                            revHom = revSeq.length();
                            break;
                        }
                    }
                }
                if (forHom > 0 || revHom > 0) {
                    if (forHom >= revHom) sequence += forSeq;
                    else sequence = revSeq + sequence;
                    int mhLen = Math.min(sequence.length() - typeLen, 5);
                    int lenClass = Math.min(typeLen, 5);
                    return lenClass + ":Del:M:" + mhLen;
                }
            }

            int lenClass = Math.min(typeLen, 5);
            int repeatClass = Math.min(sequence.length() / typeLen - 1, 5);
            if (typeLen == 1) {
                char base = typeSeq.charAt(0);
                return "1:Del:" + ((base == 'C' || base == 'G') ? "C" : "T") + ":" + repeatClass;
            } else {
                return lenClass + ":Del:R:" + repeatClass;
            }
        }

        else {
            String typeSeq = alt.substring(1);
            int typeLen = typeSeq.length();
            String sequence = typeSeq;

            int posRev = pos;
            String subSeq = "";
            if (posRev - typeLen >= 0) subSeq = seq.substring(posRev - typeLen, posRev);
            while (posRev - typeLen > 0 && subSeq.equals(typeSeq)) {
                sequence = subSeq + sequence;
                posRev -= typeLen;
                if (posRev - typeLen >= 0) subSeq = seq.substring(posRev - typeLen, posRev);
                else subSeq = "";
            }

            int posForward = pos;
            if (posForward + typeLen <= seq.length()) {
                subSeq = seq.substring(posForward, posForward + typeLen);
            } else subSeq = "";
            while (subSeq.equals(typeSeq)) {
                sequence += subSeq;
                posForward += typeLen;
                if (posForward + typeLen <= seq.length()) {
                    subSeq = seq.substring(posForward, posForward + typeLen);
                } else subSeq = "";
            }

            int repeatClass = Math.min(sequence.length() / typeLen - 1, 5);
            int lenClass = Math.min(typeLen, 5);
            if (typeLen == 1) {
                char base = typeSeq.charAt(0);
                return "1:Ins:" + ((base == 'C' || base == 'G') ? "C" : "T") + ":" + repeatClass;
            } else {
                return lenClass + ":Ins:R:" + repeatClass;
            }
        }
    }

    private CosmicData loadCosmicSignatures(Path cosmicPath) throws Exception {
        try (BufferedReader reader = Files.newBufferedReader(cosmicPath)) {
            String header = reader.readLine();
            if (header == null) throw new IllegalArgumentException("Empty COSMIC file");
            String[] headers = header.split("\t");
            List<String> sigNames = new ArrayList<>();
            for (int i = 1; i < headers.length; i++) sigNames.add(headers[i].trim());
            int nSig = sigNames.size();

            List<String> categories = new ArrayList<>();
            List<double[]> rows = new ArrayList<>();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] cols = line.split("\t");
                if (cols.length < 1 + nSig) continue;
                categories.add(cols[0].trim());
                double[] row = new double[nSig];
                for (int i = 0; i < nSig; i++) row[i] = Double.parseDouble(cols[i+1].trim());
                rows.add(row);
            }
            double[][] matrix = new double[rows.size()][nSig];
            for (int i = 0; i < rows.size(); i++) matrix[i] = rows.get(i);
            return new CosmicData(categories, sigNames, matrix);
        }
    }

    private static class CosmicData {
        List<String> categories, signatureNames;
        double[][] matrix;
        CosmicData(List<String> c, List<String> s, double[][] m) { categories=c; signatureNames=s; matrix=m; }
    }

    private double[] alignToCosmic(Map<String, Integer> matrix, List<String> cosmicCategories) {
        double[] vec = new double[cosmicCategories.size()];
        for (int i = 0; i < cosmicCategories.size(); i++) {
            Integer count = matrix.get(cosmicCategories.get(i));
            vec[i] = (count == null) ? 0.0 : count.doubleValue();
        }
        return vec;
    }

    private double[] solveNNLS(double[][] A, double[] b, int maxIter, double tol) {
        int m = A.length, n = A[0].length;
        double[] x = new double[n];
        Arrays.fill(x, 1.0);

        for (int iter = 0; iter < maxIter; iter++) {
            double[] Atb = new double[n];
            double[] AtAx = new double[n];
            for (int j = 0; j < n; j++) {
                double sumB = 0.0, sumAx = 0.0;
                for (int i = 0; i < m; i++) {
                    sumB += A[i][j] * b[i];
                    double ax = 0.0;
                    for (int k = 0; k < n; k++) ax += A[i][k] * x[k];
                    sumAx += A[i][j] * ax;
                }
                Atb[j] = sumB;
                AtAx[j] = sumAx;
            }
            double maxChange = 0.0;
            for (int j = 0; j < n; j++) {
                double old = x[j];
                if (AtAx[j] < 1e-15) x[j] = 0.0;
                else x[j] = old * Atb[j] / AtAx[j];
                if (Double.isNaN(x[j]) || Double.isInfinite(x[j])) x[j] = 0.0;
                maxChange = Math.max(maxChange, Math.abs(x[j] - old));
            }
            if (maxChange < tol) break;
        }
        return x;
    }

    private double[] reconstruct(double[][] A, double[] x) {
        double[] rec = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            double sum = 0.0;
            for (int j = 0; j < x.length; j++) sum += A[i][j] * x[j];
            rec[i] = sum;
        }
        return rec;
    }

    private double cosineSimilarity(double[] a, double[] b) {
        double dot = 0.0, na = 0.0, nb = 0.0;
        for (int i = 0; i < a.length; i++) {
            dot += a[i] * b[i];
            na += a[i] * a[i];
            nb += b[i] * b[i];
        }
        if (na == 0 || nb == 0) return 0.0;
        return dot / (Math.sqrt(na) * Math.sqrt(nb));
    }

//    private double computeRMSE(double[] observed, double[] reconstructed) {
//        double sum = 0.0;
//        for (int i = 0; i < observed.length; i++) {
//            double diff = observed[i] - reconstructed[i];
//            sum += diff * diff;
//        }
//        return Math.sqrt(sum / observed.length);
//    }

    private double computePearson(double[] x, double[] y) {
        double meanX = Arrays.stream(x).average().orElse(0.0);
        double meanY = Arrays.stream(y).average().orElse(0.0);
        double num = 0.0, denX = 0.0, denY = 0.0;
        for (int i = 0; i < x.length; i++) {
            double dx = x[i] - meanX;
            double dy = y[i] - meanY;
            num += dx * dy;
            denX += dx * dx;
            denY += dy * dy;
        }
        if (denX == 0 || denY == 0) return 0.0;
        return num / (Math.sqrt(denX) * Math.sqrt(denY));
    }

    private static class ExposureResult {
        String signature;
        double exposure;
        double percent;
        ExposureResult(String sig, double exp, double pct) { signature=sig; exposure=exp; percent=pct; }
    }

    private Map<String, Map<String, Double>> buildGroupedPercentage(
            Map<String, Integer> matrix,
            long totalMutations
    ) {
        // Define the group order (matches the Python plot)
        List<String> groupOrder = Arrays.asList(
                "1bp Deletion C", "1bp Deletion T",
                "1bp Insertion C", "1bp Insertion T",
                ">1bp Deletion 2", ">1bp Deletion 3", ">1bp Deletion 4", ">1bp Deletion 5+",
                ">1bp Insertion 2", ">1bp Insertion 3", ">1bp Insertion 4", ">1bp Insertion 5+",
                "MH 2", "MH 3", "MH 4", "MH 5+"
        );

        Map<String, Map<String, Double>> grouped = new LinkedHashMap<>();
        // Initialize empty maps for each group
        for (String group : groupOrder) {
            grouped.put(group, new LinkedHashMap<>());
        }

        // Helper to add a value to a group with a sub-key
        // sub-key is the repeat class or microhomology length (as string)
        for (Map.Entry<String, Integer> entry : matrix.entrySet()) {
            String key = entry.getKey();
            int count = entry.getValue();
            double pct = (totalMutations > 0) ? 100.0 * count / totalMutations : 0.0;
            pct = round4(pct); // round to 4 decimals

            String group = null;
            String subKey = null;

            // 1bp deletions
            if (key.matches("1:Del:[CT]:[0-5]")) {
                String base = key.split(":")[2];
                group = "1bp Deletion " + base;
                subKey = key.split(":")[3];
            }
            // 1bp insertions
            else if (key.matches("1:Ins:[CT]:[0-5]")) {
                String base = key.split(":")[2];
                group = "1bp Insertion " + base;
                subKey = key.split(":")[3];
            }
            // >1bp deletions at repeats
            else if (key.matches("[2-5]:Del:R:[0-5]")) {
                String[] parts = key.split(":");
                int length = Integer.parseInt(parts[0]);
                String groupLength = (length == 5) ? "5+" : String.valueOf(length);
                group = ">1bp Deletion " + groupLength;
                subKey = parts[3];
            }
            // >1bp insertions at repeats
            else if (key.matches("[2-5]:Ins:R:[0-5]")) {
                String[] parts = key.split(":");
                int length = Integer.parseInt(parts[0]);
                String groupLength = (length == 5) ? "5+" : String.valueOf(length);
                group = ">1bp Insertion " + groupLength;
                subKey = parts[3];
            }
            // microhomology deletions
            else if (key.matches("[2-5]:Del:M:[1-5]")) {
                String[] parts = key.split(":");
                int length = Integer.parseInt(parts[0]);
                String groupLength = (length == 5) ? "5+" : String.valueOf(length);
                group = "MH " + groupLength;
                subKey = parts[3];
            }

            if (group != null && subKey != null) {
                Map<String, Double> inner = grouped.get(group);
                if (inner != null) {
                    inner.put(subKey, pct);
                } else {
                    // If group not in order, still add it (shouldn't happen)
                    Map<String, Double> newInner = new LinkedHashMap<>();
                    newInner.put(subKey, pct);
                    grouped.put(group, newInner);
                }
            }
        }

        // Ensure all groups have at least an empty map (already done)
        return grouped;
    }
}