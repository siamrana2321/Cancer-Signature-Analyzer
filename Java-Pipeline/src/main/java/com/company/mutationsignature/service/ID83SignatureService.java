package com.company.mutationsignature.service;

import com.company.mutationsignature.model.ID83Response;
import org.apache.commons.math3.linear.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * ID83 mutational-signature service.
 */
@Service
public class ID83SignatureService {

    private static final Logger log = LoggerFactory.getLogger(ID83SignatureService.class);

    private static final String GENOME_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\genome\\GRCh38.fa";

    private static final String COSMIC_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\cosmic\\COSMIC_v3.4_ID_GRCh37.txt";

    private static final double NNLS_TOL = 1e-12;
    private static final int NNLS_MAX_ITER = 10000;
    private static final double NNLS_ADD_PENALTY = 0.05;
    private static final double NNLS_REMOVE_PENALTY = 0.01;
    private static final double INITIAL_REMOVE_PENALTY = 0.05;
    private static final double MIN_ACTIVITY = 1e-8;
    private static final double EPS = 1e-12;

    private static final Map<String, String> ETIOLOGY_MAP = new HashMap<>();

    static {
        ETIOLOGY_MAP.put("ID1", "Slippage during DNA replication of the replicated DNA strand");
        ETIOLOGY_MAP.put("ID2", "Slippage during DNA replication of the replicated DNA strand");
        ETIOLOGY_MAP.put("ID3", "Tobacco smoking");
        ETIOLOGY_MAP.put("ID4", "Unknown");
        ETIOLOGY_MAP.put("ID5", "Unknown");
        ETIOLOGY_MAP.put("ID6", "Defective homologous recombination DNA repair");
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

    // -------------------------------------------------------------------------
    // Public pipeline
    // -------------------------------------------------------------------------

    public ID83Response runPipeline(Path vcfPath) throws Exception {
        log.info("Starting ID83 pipeline: {}", vcfPath);

        Genome genome = loadGenome(Paths.get(GENOME_PATH));
        CosmicData cosmic = loadCosmicSignatures(Paths.get(COSMIC_PATH));

        validateCosmicChannels(cosmic.categories);

        MatrixBuildResult build = buildMatrixFromVcf(vcfPath, genome);
        double[] sample = alignToCosmic(build.matrix, cosmic.categories);
        long totalMutations = Math.round(sum(sample));

        if (totalMutations == 0L) {
            throw new IllegalArgumentException(
                    "No classifiable ID83 indels were found. " +
                            "Input=" + build.totalAlleles +
                            ", REF mismatch=" + build.skipReasons.get(SkipReason.REF_MISMATCH));
        }

        FitResult fit = sparseRefit(cosmic.matrix, sample);
        double[] reconstructed = multiply(cosmic.matrix, fit.exposures);

        double cosine = cosineSimilarity(sample, reconstructed);
        double pearson = pearson(sample, reconstructed);
        double l1 = l1(sample, reconstructed);
        double l2 = l2(sample, reconstructed);
        double l1Percent = 100.0 * l1 / Math.max(sum(sample), EPS);
        double sampleL2 = vectorNorm(sample);
        double l2Percent = 100.0 * l2 / Math.max(sampleL2, EPS);
        double kl = klDivergence(normalize(sample), normalize(reconstructed));

        List<ExposureResult> exposureResults = new ArrayList<>();
        double totalActivity = sum(fit.exposures);

        for (int i = 0; i < fit.exposures.length; i++) {
            if (fit.exposures[i] > MIN_ACTIVITY) {
                double pct = 100.0 * fit.exposures[i] / Math.max(totalActivity, EPS);
                exposureResults.add(new ExposureResult(
                        cosmic.signatureNames.get(i),
                        fit.exposures[i],
                        pct
                ));
            }
        }

        exposureResults.sort((a, b) -> Double.compare(b.activity, a.activity));

        List<ID83Response.CosmicContribution> contributions = new ArrayList<>();
        List<ID83Response.BiologicalAnnotation> annotations = new ArrayList<>();
        StringBuilder summary = new StringBuilder();

        for (ExposureResult er : exposureResults) {
            String etiology = ETIOLOGY_MAP.getOrDefault(er.signature, "Unknown");
            contributions.add(new ID83Response.CosmicContribution(
                    er.signature, round4(er.percentage), etiology));
            annotations.add(new ID83Response.BiologicalAnnotation(
                    er.signature, etiology, contributionLevel(er.percentage)));
            if (summary.length() > 0) summary.append(". ");
            summary.append(er.signature)
                    .append(" (")
                    .append(String.format(Locale.US, "%.1f", er.percentage))
                    .append("%): ")
                    .append(etiology);
        }

        Map<String, Double> percentage = new LinkedHashMap<>();
        for (String channel : canonicalId83Channels()) {
            int count = build.matrix.get(channel);
            percentage.put(channel, round4(100.0 * count / totalMutations));
        }

        ID83Response response = new ID83Response();
        response.setCosmicVersion("COSMIC v3.4 ID83 GRCh38");
        response.setTotalMutations(totalMutations);
        response.setId83Spectrum(build.matrix);
        response.setId83Percentage(percentage);
        response.setId83Grouped(buildGroupedPercentage(build.matrix, totalMutations));
        response.setCosmicContributions(contributions);
        response.setBiologicalAnnotations(annotations);
        response.setClinicalSummary(summary.toString());
        response.setReconstructionCosine(cosine);
        response.setPearson(pearson);
        response.setL1Norm(l1);
        response.setL1NormPercent(l1Percent);
        response.setL2Norm(l2);
        response.setL2NormPercent(l2Percent);
        response.setKlDivergence(kl);

        log.info(
                "ID83 complete: records={}, alleles={}, classified={}, skipped={}, cosine={}, signatures={}",
                build.totalRecords, build.totalAlleles, totalMutations,
                build.totalAlleles - totalMutations, cosine, exposureResults.size());

        for (Map.Entry<SkipReason, Long> e : build.skipReasons.entrySet()) {
            if (e.getValue() > 0L) {
                log.info("ID83 skip reason {} = {}", e.getKey(), e.getValue());
            }
        }

        return response;
    }

    // -------------------------------------------------------------------------
    // VCF -> normalized indels -> ID83 Profile
    // -------------------------------------------------------------------------

    private MatrixBuildResult buildMatrixFromVcf(Path vcfPath, Genome genome) throws Exception {
        LinkedHashMap<String, Integer> matrix = initializeMatrix();
        MatrixBuildResult result = new MatrixBuildResult(matrix);

        try (BufferedReader reader = Files.newBufferedReader(vcfPath)) {
            String line;

            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') {
                    continue;
                }

                result.totalRecords++;

                String[] cols = line.split("\t", -1);
                if (cols.length < 5) {
                    increment(result, SkipReason.MALFORMED_RECORD);
                    continue;
                }

                String chrom = genome.resolveChromosome(cols[0]);
                if (chrom == null) {
                    increment(result, SkipReason.CHROMOSOME_NOT_FOUND);
                    continue;
                }

                final int pos1;
                try {
                    pos1 = Integer.parseInt(cols[1]);
                } catch (NumberFormatException ex) {
                    increment(result, SkipReason.MALFORMED_RECORD);
                    continue;
                }

                String ref = cols[3].toUpperCase(Locale.ROOT);
                String altField = cols[4].toUpperCase(Locale.ROOT);

                if (altField.contains(",")) {
                    result.totalAlleles++;
                    increment(result, SkipReason.MULTIALLELIC);
                    continue;
                }

                String alt = altField;
                result.totalAlleles++;

                if (!isSequenceAllele(ref) || !isSequenceAllele(alt)) {
                    increment(result, SkipReason.SYMBOLIC_OR_INVALID_ALLELE);
                    continue;
                }

                if (ref.length() == alt.length()) {
                    increment(result, SkipReason.NOT_INDEL);
                    continue;
                }

                NormalizedIndel indel = normalizeIndel(
                        chrom, pos1, ref, alt, genome.sequence(chrom));

                if (indel == null || indel.deleted.length() == indel.inserted.length()) {
                    increment(result, SkipReason.COMPLEX_VARIANT);
                    continue;
                }

                if (!validateReference(indel, genome.sequence(chrom))) {
                    increment(result, SkipReason.REF_MISMATCH);
                    continue;
                }

                String channel = classifyId83(indel, genome.sequence(chrom));

                if (channel == null || !matrix.containsKey(channel)) {
                    increment(result, SkipReason.UNCLASSIFIED);
                    continue;
                }

                matrix.put(channel, matrix.get(channel) + 1);
            }
        }

        return result;
    }

    private NormalizedIndel normalizeIndel(
            String chrom, int pos1, String rawRef, String rawAlt, String chromosome) {

        int start0 = pos1 - 1;
        String ref = rawRef;
        String alt = rawAlt;

        while (ref.length() > 1 && alt.length() > 1
                && ref.charAt(ref.length() - 1) == alt.charAt(alt.length() - 1)) {
            ref = ref.substring(0, ref.length() - 1);
            alt = alt.substring(0, alt.length() - 1);
        }

        int prefix = 0;
        int min = Math.min(ref.length(), alt.length());
        while (prefix < min && ref.charAt(prefix) == alt.charAt(prefix)) {
            prefix++;
        }

        start0 += prefix;
        ref = ref.substring(prefix);
        alt = alt.substring(prefix);

        if (!ref.isEmpty() && !alt.isEmpty()) {
            return null;
        }

        String deleted = ref;
        String inserted = alt;

        if (!deleted.isEmpty()) {
            while (start0 > 0) {
                char previous = chromosome.charAt(start0 - 1);
                char last = deleted.charAt(deleted.length() - 1);
                if (previous != last) {
                    break;
                }
                deleted = previous + deleted.substring(0, deleted.length() - 1);
                start0--;
            }
        } else if (!inserted.isEmpty()) {
            while (start0 > 0) {
                char previous = chromosome.charAt(start0 - 1);
                char last = inserted.charAt(inserted.length() - 1);
                if (previous != last) {
                    break;
                }
                inserted = previous + inserted.substring(0, inserted.length() - 1);
                start0--;
            }
        }

        return new NormalizedIndel(chrom, start0, deleted, inserted);
    }

    private boolean validateReference(NormalizedIndel indel, String chromosome) {
        if (indel.deleted.isEmpty()) {
            return indel.start0 >= 0 && indel.start0 <= chromosome.length();
        }

        int end = indel.start0 + indel.deleted.length();
        return indel.start0 >= 0
                && end <= chromosome.length()
                && chromosome.substring(indel.start0, end).equals(indel.deleted);
    }

    private String classifyId83(NormalizedIndel indel, String chromosome) {
        if (indel.isDeletion()) {
            return classifyDeletion(indel, chromosome);
        }

        if (indel.isInsertion()) {
            return classifyInsertion(indel, chromosome);
        }

        return null;
    }

    private String classifyDeletion(NormalizedIndel indel, String chromosome) {
        String deleted = indel.deleted;
        int eventLength = deleted.length();

        if (eventLength == 1) {
            char base = pyrimidineClass(deleted.charAt(0));
            int repeatUnits = homopolymerRepeatUnits(
                    chromosome, indel.start0, deleted.charAt(0), true);
            int repeatClass = cap(repeatUnits - 1, 0, 5);
            return "1:Del:" + base + ":" + repeatClass;
        }

        int repeatUnits = tandemRepeatUnits(
                chromosome, indel.start0, deleted, true);

        if (repeatUnits > 1) {
            return Math.min(eventLength, 5)
                    + ":Del:R:"
                    + cap(repeatUnits - 1, 0, 5);
        }

        int microhomology = deletionJunctionMicrohomology(indel, chromosome);

        if (microhomology > 0) {
            int lengthClass = Math.min(eventLength, 5);
            int maxMhForClass = lengthClass == 2 ? 1
                    : lengthClass == 3 ? 2
                    : lengthClass == 4 ? 3 : 5;
            int mhClass = Math.min(microhomology, maxMhForClass);
            return lengthClass + ":Del:M:" + mhClass;
        }

        return Math.min(eventLength, 5) + ":Del:R:0";
    }

    private String classifyInsertion(NormalizedIndel indel, String chromosome) {
        String inserted = indel.inserted;
        int eventLength = inserted.length();

        if (eventLength == 1) {
            char base = pyrimidineClass(inserted.charAt(0));
            int repeatUnits = insertionRepeatUnits(
                    chromosome, indel.start0, inserted);
            int repeatClass = cap(repeatUnits, 0, 5);
            return "1:Ins:" + base + ":" + repeatClass;
        }

        int repeatUnits = insertionRepeatUnits(
                chromosome, indel.start0, inserted);

        return Math.min(eventLength, 5)
                + ":Ins:R:"
                + cap(repeatUnits, 0, 5);
    }

    private int homopolymerRepeatUnits(
            String chromosome, int start0, char base, boolean deletion) {

        int left = start0 - 1;
        int right = deletion ? start0 + 1 : start0;
        int count = 1;

        while (left >= 0 && chromosome.charAt(left) == base) {
            count++;
            left--;
        }

        while (right < chromosome.length() && chromosome.charAt(right) == base) {
            count++;
            right++;
        }

        return count;
    }

    private int tandemRepeatUnits(
            String chromosome, int start0, String unit, boolean deletion) {

        int length = unit.length();
        int copies = 1;

        int p = start0 - length;
        while (p >= 0 && chromosome.substring(p, p + length).equals(unit)) {
            copies++;
            p -= length;
        }

        p = deletion ? start0 + length : start0;
        while (p + length <= chromosome.length()
                && chromosome.substring(p, p + length).equals(unit)) {
            copies++;
            p += length;
        }

        return copies;
    }

    private int insertionRepeatUnits(String chromosome, int start0, String unit) {
        int length = unit.length();
        int copies = 0;

        int p = start0 - length;
        while (p >= 0 && chromosome.substring(p, p + length).equals(unit)) {
            copies++;
            p -= length;
        }

        p = start0;
        while (p + length <= chromosome.length()
                && chromosome.substring(p, p + length).equals(unit)) {
            copies++;
            p += length;
        }

        return copies;
    }

    private int deletionJunctionMicrohomology(
            NormalizedIndel indel, String chromosome) {

        int start = indel.start0;
        int deletionLength = indel.deleted.length();
        int end = start + deletionLength;
        int maxMicrohomology = Math.min(deletionLength - 1, 5);

        int forwardMh = 0;
        for (int i = 0; i < maxMicrohomology; i++) {
            int rightIndex = end + i;
            if (rightIndex >= chromosome.length()) {
                break;
            }
            if (indel.deleted.charAt(i) != chromosome.charAt(rightIndex)) {
                break;
            }
            forwardMh++;
        }

        int reverseMh = 0;
        for (int i = 0; i < maxMicrohomology; i++) {
            int deletedIndex = deletionLength - 1 - i;
            int leftIndex = start - 1 - i;
            if (leftIndex < 0) {
                break;
            }
            if (indel.deleted.charAt(deletedIndex) != chromosome.charAt(leftIndex)) {
                break;
            }
            reverseMh++;
        }

        return Math.max(forwardMh, reverseMh);
    }

    private char pyrimidineClass(char base) {
        switch (base) {
            case 'C':
            case 'G':
                return 'C';
            case 'A':
            case 'T':
                return 'T';
            default:
                throw new IllegalArgumentException("Unexpected DNA base: " + base);
        }
    }

    private boolean isSequenceAllele(String allele) {
        if (allele == null || allele.isEmpty() || allele.equals(".")) {
            return false;
        }

        for (int i = 0; i < allele.length(); i++) {
            char c = allele.charAt(i);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                return false;
            }
        }

        return true;
    }

    // -------------------------------------------------------------------------
    // Canonical ID83 channels
    // -------------------------------------------------------------------------

    private LinkedHashMap<String, Integer> initializeMatrix() {
        LinkedHashMap<String, Integer> matrix = new LinkedHashMap<>();

        for (String channel : canonicalId83Channels()) {
            matrix.put(channel, 0);
        }

        if (matrix.size() != 83) {
            throw new IllegalStateException(
                    "Internal ID83 channel definition contains " + matrix.size() + " channels");
        }

        return matrix;
    }

    private List<String> canonicalId83Channels() {
        List<String> channels = new ArrayList<>(83);

        for (String op : Arrays.asList("Del", "Ins")) {
            for (String base : Arrays.asList("C", "T")) {
                for (int repeat = 0; repeat <= 5; repeat++) {
                    channels.add("1:" + op + ":" + base + ":" + repeat);
                }
            }
        }

        for (int length = 2; length <= 5; length++) {
            for (String op : Arrays.asList("Del", "Ins")) {
                for (int repeat = 0; repeat <= 5; repeat++) {
                    channels.add(length + ":" + op + ":R:" + repeat);
                }
            }
        }

        channels.add("2:Del:M:1");
        channels.add("3:Del:M:1");
        channels.add("3:Del:M:2");
        channels.add("4:Del:M:1");
        channels.add("4:Del:M:2");
        channels.add("4:Del:M:3");

        for (int mh = 1; mh <= 5; mh++) {
            channels.add("5:Del:M:" + mh);
        }

        return channels;
    }

    private void validateCosmicChannels(List<String> categories) {
        Set<String> expected = new LinkedHashSet<>(canonicalId83Channels());
        Set<String> actual = new LinkedHashSet<>(categories);

        if (categories.size() != 83 || !expected.equals(actual)) {
            Set<String> missing = new LinkedHashSet<>(expected);
            missing.removeAll(actual);

            Set<String> unexpected = new LinkedHashSet<>(actual);
            unexpected.removeAll(expected);

            throw new IllegalArgumentException(
                    "COSMIC ID83 categories do not match canonical 83 channels. " +
                            "Rows=" + categories.size() +
                            ", missing=" + missing +
                            ", unexpected=" + unexpected);
        }
    }

    // -------------------------------------------------------------------------
    // SigProfilerAssignment-compatible sparse assignment
    // -------------------------------------------------------------------------

    private FitResult sparseRefit(double[][] signatures, double[] sample) {
        LinkedHashSet<Integer> all = new LinkedHashSet<>();
        for (int i = 0; i < signatures[0].length; i++) {
            all.add(i);
        }

        // single_sample.fit_signatures(...): NNLS -> normalize to mutation count
        // -> roundConserveSum.
        double[] initial = fitSubsetConserveMutations(signatures, sample, all);

        // Initial backward elimination uses initial_remove_penalty.
        initial = removeAllSingleSignatures(
                signatures, sample, initial, INITIAL_REMOVE_PENALTY);

        LinkedHashSet<Integer> background = positiveIndices(initial);

        double originalDistance = Double.POSITIVE_INFINITY;
        double[] finalActivities = initial;

        while (true) {
            double layerBestDistance = Double.POSITIVE_INFINITY;
            double[] layerBestActivities = null;

            for (int candidate = 0; candidate < signatures[0].length; candidate++) {
                if (background.contains(candidate)) {
                    continue;
                }

                AddResult add = addSingleSignatureSigProfiler(
                        signatures,
                        sample,
                        background,
                        candidate,
                        NNLS_ADD_PENALTY);

                double[] removed = removeAllSingleSignatures(
                        signatures,
                        sample,
                        add.exposures,
                        NNLS_REMOVE_PENALTY);

                double removeDistance = normalizedL2Distance(
                        sample, multiply(signatures, removed));

                boolean sameCompositionAccordingToSigProfiler =
                        sigProfilerCompositionCheck(add.exposures, removed);

                double[] activities = sameCompositionAccordingToSigProfiler
                        ? add.exposures
                        : removed;

                double distance = sameCompositionAccordingToSigProfiler
                        ? add.distance
                        : removeDistance;

                if (distance < layerBestDistance) {
                    layerBestDistance = distance;
                    layerBestActivities = activities;
                }
            }

            if (layerBestActivities != null
                    && layerBestDistance < originalDistance) {
                originalDistance = layerBestDistance;
                finalActivities = layerBestActivities;
                background = positiveIndices(finalActivities);
            } else {
                break;
            }
        }

        finalActivities = roundConserveSum(
                finalActivities, Math.round(sum(sample)));

        return new FitResult(
                finalActivities,
                positiveIndices(finalActivities));
    }

    private AddResult addSingleSignatureSigProfiler(
            double[][] signatures,
            double[] sample,
            Set<Integer> present,
            int candidate,
            double cutoff) {

        LinkedHashSet<Integer> current = new LinkedHashSet<>(present);

        double originalDistance = Double.POSITIVE_INFINITY;
        double[] originalExposure = new double[signatures[0].length];

        if (!current.isEmpty()) {
            SubsetFit originalFit =
                    fitSubsetForAdd(signatures, sample, current);
            originalExposure = originalFit.exposures;
            originalDistance = originalFit.rawDistance;
        }

        LinkedHashSet<Integer> trial = new LinkedHashSet<>(current);
        trial.add(candidate);

        SubsetFit trialFit = fitSubsetForAdd(signatures, sample, trial);

        if (originalDistance - trialFit.rawDistance > cutoff) {
            return new AddResult(
                    trialFit.exposures,
                    trialFit.rawDistance);
        }

        return new AddResult(originalExposure, originalDistance);
    }

    private SubsetFit fitSubsetForAdd(
            double[][] fullMatrix,
            double[] sample,
            Set<Integer> active) {

        int signatureCount = fullMatrix[0].length;
        double[] fullExposure = new double[signatureCount];

        if (active.isEmpty()) {
            return new SubsetFit(
                    fullExposure,
                    Double.POSITIVE_INFINITY);
        }

        List<Integer> indices = new ArrayList<>(active);
        Collections.sort(indices);

        double[][] subset =
                new double[fullMatrix.length][indices.size()];

        for (int row = 0; row < fullMatrix.length; row++) {
            for (int col = 0; col < indices.size(); col++) {
                subset[row][col] =
                        fullMatrix[row][indices.get(col)];
            }
        }

        double[] rawWeights = solveNnls(subset, sample);
        double[] rawReconstruction = multiply(subset, rawWeights);
        double rawDistance =
                normalizedL2Distance(sample, rawReconstruction);

        double[] normalized = normalizeExposureToMutationCount(
                rawWeights, sum(sample));

        // np.round(...) + correction of the largest coefficient.
        double[] rounded = numpyRoundAndConserve(
                normalized, Math.round(sum(sample)));

        for (int i = 0; i < indices.size(); i++) {
            fullExposure[indices.get(i)] = rounded[i];
        }

        return new SubsetFit(fullExposure, rawDistance);
    }

    private double[] removeAllSingleSignatures(
            double[][] signatures,
            double[] sample,
            double[] exposures,
            double cutoff) {

        double[] oldExposures =
                roundConserveSum(exposures, Math.round(sum(sample)));

        LinkedHashSet<Integer> current = positiveIndices(oldExposures);

        if (current.size() <= 1) {
            return oldExposures;
        }

        double[] baseFit =
                fitSubsetRaw(signatures, sample, current);

        double originalDistance = normalizedL2Distance(
                sample, multiply(signatures, baseFit));

        double[] successful = null;

        while (current.size() > 1) {
            double bestDifference = Double.POSITIVE_INFINITY;
            double bestDistance = Double.POSITIVE_INFINITY;
            double[] bestExposure = null;

            for (Integer candidate : new ArrayList<>(current)) {
                LinkedHashSet<Integer> trial =
                        new LinkedHashSet<>(current);
                trial.remove(candidate);

                if (trial.isEmpty()) {
                    continue;
                }

                double[] trialRaw =
                        fitSubsetRaw(signatures, sample, trial);

                double trialDistance = normalizedL2Distance(
                        sample, multiply(signatures, trialRaw));

                double difference =
                        trialDistance - originalDistance;

                if (difference < 0.0) {
                    difference = cutoff + 1e-100;
                }

                if (difference < bestDifference) {
                    bestDifference = difference;
                    bestDistance = trialDistance;
                    bestExposure =
                            normalizeExposureToMutationCount(
                                    trialRaw, sum(sample));
                }
            }

            if (bestExposure == null || bestDifference > cutoff) {
                break;
            }

            successful = bestExposure;
            current = positiveIndices(successful);
            originalDistance = bestDistance;

            if (current.size() <= 1) {
                break;
            }
        }

        if (successful == null) {
            return oldExposures;
        }

        return successful;
    }

    private double[] fitSubsetRaw(
            double[][] fullMatrix,
            double[] sample,
            Set<Integer> active) {

        int signatureCount = fullMatrix[0].length;
        double[] fullExposure = new double[signatureCount];

        if (active.isEmpty()) {
            return fullExposure;
        }

        List<Integer> indices = new ArrayList<>(active);
        Collections.sort(indices);

        double[][] subset =
                new double[fullMatrix.length][indices.size()];

        for (int row = 0; row < fullMatrix.length; row++) {
            for (int col = 0; col < indices.size(); col++) {
                subset[row][col] =
                        fullMatrix[row][indices.get(col)];
            }
        }

        double[] subsetExposure = solveNnls(subset, sample);

        for (int i = 0; i < indices.size(); i++) {
            fullExposure[indices.get(i)] = subsetExposure[i];
        }

        return fullExposure;
    }

    private double[] fitSubsetConserveMutations(
            double[][] fullMatrix,
            double[] sample,
            Set<Integer> active) {

        double[] raw = fitSubsetRaw(fullMatrix, sample, active);
        double[] normalized =
                normalizeExposureToMutationCount(raw, sum(sample));

        return roundConserveSum(
                normalized, Math.round(sum(sample)));
    }

    private double[] normalizeExposureToMutationCount(
            double[] exposure,
            double mutationCount) {

        double total = sum(exposure);
        double[] result = new double[exposure.length];

        if (total <= EPS) {
            return result;
        }

        for (int i = 0; i < exposure.length; i++) {
            result[i] =
                    exposure[i] / total * mutationCount;
        }

        return result;
    }

    private double[] numpyRoundAndConserve(
            double[] values,
            long targetTotal) {

        double[] result = new double[values.length];

        for (int i = 0; i < values.length; i++) {
            result[i] = Math.rint(Math.max(0.0, values[i]));
        }

        long difference =
                targetTotal - Math.round(sum(result));

        if (difference != 0L && result.length > 0) {
            int maxIndex = indexOfMaximum(result);
            result[maxIndex] += difference;
        }

        return result;
    }

    private boolean sigProfilerCompositionCheck(
            double[] a,
            double[] b) {

        List<Integer> ai = positiveIndexList(a);
        List<Integer> bi = positiveIndexList(b);

        return ai.size() == bi.size()
                && numpyAllOfIndices(ai) == numpyAllOfIndices(bi);
    }

    private List<Integer> positiveIndexList(double[] exposures) {
        List<Integer> result = new ArrayList<>();

        for (int i = 0; i < exposures.length; i++) {
            if (exposures[i] > MIN_ACTIVITY) {
                result.add(i);
            }
        }

        return result;
    }

    private boolean numpyAllOfIndices(List<Integer> indices) {
        // np.all(empty integer array) is True.
        for (Integer index : indices) {
            if (index == 0) {
                return false;
            }
        }

        return true;
    }

    private double[] roundConserveSum(
            double[] values,
            long targetTotal) {

        double[] result = new double[values.length];
        List<Integer> order = new ArrayList<>();

        long ceilTotal = 0L;

        for (int i = 0; i < values.length; i++) {
            result[i] = Math.ceil(Math.max(0.0, values[i]));
            ceilTotal += Math.round(result[i]);
            order.add(i);
        }

        Collections.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                double residualA = values[a] - result[a];
                double residualB = values[b] - result[b];
                return Double.compare(residualA, residualB);
            }
        });

        long remove = ceilTotal - targetTotal;

        for (int i = 0; i < remove && i < order.size(); i++) {
            int index = order.get(i);
            if (result[index] > 0.0) {
                result[index] -= 1.0;
            }
        }

        if (Math.round(sum(result)) != targetTotal) {
            throw new IllegalStateException(
                    "Mutation count not conserved during roundConserveSum");
        }

        return result;
    }

    private int indexOfMaximum(double[] values) {
        int index = 0;

        for (int i = 1; i < values.length; i++) {
            if (values[i] > values[index]) {
                index = i;
            }
        }

        return index;
    }

    private double normalizedL2Distance(
            double[] observed,
            double[] reconstructed) {

        return l2(observed, reconstructed)
                / Math.max(vectorNorm(observed), EPS);
    }

    private LinkedHashSet<Integer> positiveIndices(
            double[] exposures) {

        LinkedHashSet<Integer> result = new LinkedHashSet<>();

        for (int i = 0; i < exposures.length; i++) {
            if (exposures[i] > MIN_ACTIVITY) {
                result.add(i);
            }
        }

        return result;
    }

    // -------------------------------------------------------------------------
    // NNLS
    // -------------------------------------------------------------------------

    private double[] solveNnls(double[][] a, double[] b) {
        int m = a.length;
        int n = a[0].length;

        boolean[] skipCol = new boolean[n];

        for (int j = 0; j < n; j++) {
            boolean allZero = true;

            for (int i = 0; i < m; i++) {
                if (Math.abs(a[i][j]) > EPS) {
                    allZero = false;
                    break;
                }
            }

            skipCol[j] = allZero;
        }

        int nActive = 0;
        for (boolean skip : skipCol) {
            if (!skip) {
                nActive++;
            }
        }

        if (nActive == 0) {
            return new double[n];
        }

        double[][] reduced = new double[m][nActive];
        int[] columnMap = new int[nActive];

        int reducedColumn = 0;

        for (int j = 0; j < n; j++) {
            if (!skipCol[j]) {
                for (int i = 0; i < m; i++) {
                    reduced[i][reducedColumn] = a[i][j];
                }

                columnMap[reducedColumn] = j;
                reducedColumn++;
            }
        }

        double[] reducedSolution =
                nnlsLawsonHanson(reduced, b);

        double[] solution = new double[n];

        for (int i = 0; i < nActive; i++) {
            solution[columnMap[i]] =
                    Math.max(0.0, reducedSolution[i]);
        }

        return solution;
    }

    private double[] nnlsLawsonHanson(
            double[][] a,
            double[] b) {

        int m = a.length;
        int n = a[0].length;

        double[] x = new double[n];
        boolean[] passive = new boolean[n];
        double[] w = new double[n];
        double[] residual = new double[m];

        computeResidualsAndW(a, b, x, residual, w);

        int iteration = 0;

        while (true) {
            int entering = -1;
            double maximumW = NNLS_TOL;

            for (int j = 0; j < n; j++) {
                if (!passive[j] && w[j] > maximumW) {
                    maximumW = w[j];
                    entering = j;
                }
            }

            if (entering == -1) {
                break;
            }

            passive[entering] = true;

            while (true) {
                List<Integer> passiveList = new ArrayList<>();

                for (int j = 0; j < n; j++) {
                    if (passive[j]) {
                        passiveList.add(j);
                    }
                }

                if (passiveList.isEmpty()) {
                    break;
                }

                double[] z =
                        solvePassiveLeastSquares(a, b, passiveList);

                boolean allPositive = true;

                for (double value : z) {
                    if (value <= NNLS_TOL) {
                        allPositive = false;
                        break;
                    }
                }

                if (allPositive) {
                    Arrays.fill(x, 0.0);

                    for (int k = 0; k < passiveList.size(); k++) {
                        x[passiveList.get(k)] = z[k];
                    }

                    break;
                }

                double alpha = Double.POSITIVE_INFINITY;

                for (int k = 0; k < passiveList.size(); k++) {
                    int column = passiveList.get(k);

                    if (z[k] <= NNLS_TOL) {
                        double denominator = x[column] - z[k];

                        if (denominator > NNLS_TOL) {
                            alpha = Math.min(
                                    alpha,
                                    x[column] / denominator);
                        }
                    }
                }

                if (!Double.isFinite(alpha)) {
                    for (int k = 0; k < passiveList.size(); k++) {
                        if (z[k] <= NNLS_TOL) {
                            int column = passiveList.get(k);
                            passive[column] = false;
                            x[column] = 0.0;
                        }
                    }

                    continue;
                }

                for (int k = 0; k < passiveList.size(); k++) {
                    int column = passiveList.get(k);
                    x[column] +=
                            alpha * (z[k] - x[column]);
                }

                for (int column : passiveList) {
                    if (x[column] <= NNLS_TOL) {
                        passive[column] = false;
                        x[column] = 0.0;
                    }
                }
            }

            computeResidualsAndW(a, b, x, residual, w);

            iteration++;

            if (iteration > NNLS_MAX_ITER) {
                throw new IllegalStateException(
                        "NNLS failed to converge within "
                                + NNLS_MAX_ITER + " iterations");
            }
        }

        for (int j = 0; j < n; j++) {
            if (x[j] < MIN_ACTIVITY) {
                x[j] = 0.0;
            }
        }

        return x;
    }

    private double[] solvePassiveLeastSquares(
            double[][] a,
            double[] b,
            List<Integer> passiveList) {

        int passiveCount = passiveList.size();

        if (passiveCount == 0) {
            return new double[0];
        }

        RealMatrix passiveMatrix =
                new Array2DRowRealMatrix(a.length, passiveCount);

        for (int i = 0; i < a.length; i++) {
            for (int k = 0; k < passiveCount; k++) {
                passiveMatrix.setEntry(
                        i,
                        k,
                        a[i][passiveList.get(k)]);
            }
        }

        RealVector target = new ArrayRealVector(b, false);

        DecompositionSolver solver =
                new QRDecomposition(passiveMatrix, 1e-15)
                        .getSolver();

        if (!solver.isNonSingular()) {
            solver = new SingularValueDecomposition(passiveMatrix)
                    .getSolver();
        }

        return solver.solve(target).toArray();
    }

    private void computeResidualsAndW(
            double[][] a,
            double[] b,
            double[] x,
            double[] residual,
            double[] w) {

        Arrays.fill(residual, 0.0);

        for (int i = 0; i < a.length; i++) {
            double fitted = 0.0;

            for (int j = 0; j < a[0].length; j++) {
                fitted += a[i][j] * x[j];
            }

            residual[i] = b[i] - fitted;
        }

        Arrays.fill(w, 0.0);

        for (int j = 0; j < a[0].length; j++) {
            double value = 0.0;

            for (int i = 0; i < a.length; i++) {
                value += a[i][j] * residual[i];
            }

            w[j] = value;
        }
    }

    // -------------------------------------------------------------------------
    // FASTA / COSMIC
    // -------------------------------------------------------------------------

    private Genome loadGenome(Path fastaPath) throws Exception {
        LinkedHashMap<String, String> sequences = new LinkedHashMap<>();

        try (BufferedReader reader = Files.newBufferedReader(fastaPath)) {
            String line;
            String chromosome = null;
            StringBuilder sequence = new StringBuilder();

            while ((line = reader.readLine()) != null) {
                line = line.trim();

                if (line.isEmpty()) {
                    continue;
                }

                if (line.charAt(0) == '>') {
                    if (chromosome != null) {
                        sequences.put(chromosome, sequence.toString().toUpperCase(Locale.ROOT));
                    }

                    chromosome = line.substring(1).split("\\s+")[0];
                    sequence = new StringBuilder();
                } else {
                    sequence.append(line);
                }
            }

            if (chromosome != null) {
                sequences.put(chromosome, sequence.toString().toUpperCase(Locale.ROOT));
            }
        }

        if (sequences.isEmpty()) {
            throw new IllegalArgumentException("No FASTA sequences loaded from " + fastaPath);
        }

        return new Genome(sequences);
    }

    private CosmicData loadCosmicSignatures(Path cosmicPath) throws Exception {
        try (BufferedReader reader = Files.newBufferedReader(cosmicPath)) {
            String header = reader.readLine();

            if (header == null) {
                throw new IllegalArgumentException("Empty COSMIC signature file");
            }

            String[] headerColumns = header.split("\t", -1);
            List<String> signatures = new ArrayList<>();

            for (int i = 1; i < headerColumns.length; i++) {
                signatures.add(headerColumns[i].trim());
            }

            List<String> categories = new ArrayList<>();
            List<double[]> rows = new ArrayList<>();

            String line;

            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }

                String[] columns = line.split("\t", -1);

                if (columns.length != signatures.size() + 1) {
                    throw new IllegalArgumentException(
                            "Invalid COSMIC row column count: " + line);
                }

                categories.add(columns[0].trim());

                double[] row = new double[signatures.size()];

                for (int i = 0; i < signatures.size(); i++) {
                    row[i] = Double.parseDouble(columns[i + 1].trim());
                }

                rows.add(row);
            }

            double[][] matrix = new double[rows.size()][];

            for (int i = 0; i < rows.size(); i++) {
                matrix[i] = rows.get(i);
            }

            return new CosmicData(categories, signatures, matrix);
        }
    }

    private double[] alignToCosmic(
            Map<String, Integer> matrix, List<String> cosmicCategories) {

        double[] result = new double[cosmicCategories.size()];

        for (int i = 0; i < cosmicCategories.size(); i++) {
            Integer value = matrix.get(cosmicCategories.get(i));

            if (value == null) {
                throw new IllegalArgumentException(
                        "Missing ID83 channel: " + cosmicCategories.get(i));
            }

            result[i] = value.doubleValue();
        }

        return result;
    }

    // -------------------------------------------------------------------------
    // Metrics
    // -------------------------------------------------------------------------

    private double[] multiply(double[][] a, double[] x) {
        double[] result = new double[a.length];

        for (int row = 0; row < a.length; row++) {
            double value = 0.0;

            for (int col = 0; col < x.length; col++) {
                value += a[row][col] * x[col];
            }

            result[row] = value;
        }

        return result;
    }

    private double cosineSimilarity(double[] a, double[] b) {
        double dot = 0.0;
        double aa = 0.0;
        double bb = 0.0;

        for (int i = 0; i < a.length; i++) {
            dot += a[i] * b[i];
            aa += a[i] * a[i];
            bb += b[i] * b[i];
        }

        if (aa <= EPS || bb <= EPS) {
            return 0.0;
        }

        return dot / (Math.sqrt(aa) * Math.sqrt(bb));
    }

    private double pearson(double[] a, double[] b) {
        double meanA = sum(a) / a.length;
        double meanB = sum(b) / b.length;
        double numerator = 0.0;
        double da2 = 0.0;
        double db2 = 0.0;

        for (int i = 0; i < a.length; i++) {
            double da = a[i] - meanA;
            double db = b[i] - meanB;
            numerator += da * db;
            da2 += da * da;
            db2 += db * db;
        }

        if (da2 <= EPS || db2 <= EPS) {
            return 0.0;
        }

        return numerator / Math.sqrt(da2 * db2);
    }

    private double l1(double[] a, double[] b) {
        double result = 0.0;

        for (int i = 0; i < a.length; i++) {
            result += Math.abs(a[i] - b[i]);
        }

        return result;
    }

    private double l2(double[] a, double[] b) {
        double result = 0.0;

        for (int i = 0; i < a.length; i++) {
            double d = a[i] - b[i];
            result += d * d;
        }

        return Math.sqrt(result);
    }

    private double vectorNorm(double[] a) {
        double value = 0.0;

        for (double v : a) {
            value += v * v;
        }

        return Math.sqrt(value);
    }

    private double[] normalize(double[] values) {
        double total = sum(values);
        double[] result = new double[values.length];

        if (total <= EPS) {
            return result;
        }

        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] / total;
        }

        return result;
    }

    private double klDivergence(double[] p, double[] q) {
        double value = 0.0;

        for (int i = 0; i < p.length; i++) {
            if (p[i] <= 0.0) {
                continue;
            }

            double pi = Math.max(p[i], EPS);
            double qi = Math.max(q[i], EPS);
            value += pi * Math.log(pi / qi);
        }

        return value;
    }

    private double sum(double[] values) {
        double result = 0.0;

        for (double value : values) {
            result += value;
        }

        return result;
    }

    // -------------------------------------------------------------------------
    // Response grouping
    // -------------------------------------------------------------------------

    private Map<String, Map<String, Double>> buildGroupedPercentage(
            Map<String, Integer> matrix, long totalMutations) {

        List<String> groupOrder = Arrays.asList(
                "1bp Deletion C", "1bp Deletion T",
                "1bp Insertion C", "1bp Insertion T",
                ">1bp Deletion 2", ">1bp Deletion 3",
                ">1bp Deletion 4", ">1bp Deletion 5+",
                ">1bp Insertion 2", ">1bp Insertion 3",
                ">1bp Insertion 4", ">1bp Insertion 5+",
                "MH 2", "MH 3", "MH 4", "MH 5+"
        );

        Map<String, Map<String, Double>> grouped =
                new LinkedHashMap<>();

        for (String group : groupOrder) {
            grouped.put(group, new LinkedHashMap<>());
        }

        for (Map.Entry<String, Integer> entry : matrix.entrySet()) {
            String[] parts = entry.getKey().split(":");
            String group = null;
            String subKey = parts[3];

            if ("1".equals(parts[0])) {
                group = "Del".equals(parts[1])
                        ? "1bp Deletion " + parts[2]
                        : "1bp Insertion " + parts[2];
            } else if ("R".equals(parts[2])) {
                String length = "5".equals(parts[0]) ? "5+" : parts[0];
                group = "Del".equals(parts[1])
                        ? ">1bp Deletion " + length
                        : ">1bp Insertion " + length;
            } else if ("M".equals(parts[2])) {
                String length = "5".equals(parts[0]) ? "5+" : parts[0];
                group = "MH " + length;
            }

            if (group != null && grouped.containsKey(group)) {
                double pct = 100.0 * entry.getValue() / Math.max(totalMutations, 1L);
                grouped.get(group).put(subKey, round4(pct));
            }
        }

        return grouped;
    }

    private String contributionLevel(double percentage) {
        if (percentage >= 20.0) {
            return "Dominant";
        }

        if (percentage >= 5.0) {
            return "Intermediate";
        }

        return "Minor";
    }

    private int cap(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }

    private double round4(double value) {
        return Math.round(value * 10000.0) / 10000.0;
    }

    private void increment(MatrixBuildResult result, SkipReason reason) {
        result.skipReasons.put(reason, result.skipReasons.get(reason) + 1L);
    }

    // -------------------------------------------------------------------------
    // Internal models
    // -------------------------------------------------------------------------

    private enum SkipReason {
        MALFORMED_RECORD,
        CHROMOSOME_NOT_FOUND,
        SYMBOLIC_OR_INVALID_ALLELE,
        MULTIALLELIC,
        NOT_INDEL,
        COMPLEX_VARIANT,
        REF_MISMATCH,
        UNCLASSIFIED
    }

    private static final class NormalizedIndel {
        private final String chromosome;
        private final int start0;
        private final String deleted;
        private final String inserted;

        private NormalizedIndel(
                String chromosome, int start0, String deleted, String inserted) {
            this.chromosome = chromosome;
            this.start0 = start0;
            this.deleted = deleted;
            this.inserted = inserted;
        }

        private boolean isDeletion() {
            return !deleted.isEmpty() && inserted.isEmpty();
        }

        private boolean isInsertion() {
            return deleted.isEmpty() && !inserted.isEmpty();
        }
    }

    private static final class MatrixBuildResult {
        private final LinkedHashMap<String, Integer> matrix;
        private final EnumMap<SkipReason, Long> skipReasons =
                new EnumMap<>(SkipReason.class);
        private long totalRecords;
        private long totalAlleles;

        private MatrixBuildResult(LinkedHashMap<String, Integer> matrix) {
            this.matrix = matrix;

            for (SkipReason reason : SkipReason.values()) {
                skipReasons.put(reason, 0L);
            }
        }
    }

    private static final class Genome {
        private final Map<String, String> sequences;

        private Genome(Map<String, String> sequences) {
            this.sequences = sequences;
        }

        private String resolveChromosome(String input) {
            if (sequences.containsKey(input)) {
                return input;
            }

            if (input.startsWith("chr")) {
                String withoutChr = input.substring(3);

                if (sequences.containsKey(withoutChr)) {
                    return withoutChr;
                }
            } else {
                String withChr = "chr" + input;

                if (sequences.containsKey(withChr)) {
                    return withChr;
                }
            }

            return null;
        }

        private String sequence(String chromosome) {
            return sequences.get(chromosome);
        }
    }

    private static final class CosmicData {
        private final List<String> categories;
        private final List<String> signatureNames;
        private final double[][] matrix;

        private CosmicData(
                List<String> categories,
                List<String> signatureNames,
                double[][] matrix) {
            this.categories = categories;
            this.signatureNames = signatureNames;
            this.matrix = matrix;
        }
    }

    private static final class ExposureResult {
        private final String signature;
        private final double activity;
        private final double percentage;

        private ExposureResult(
                String signature, double activity, double percentage) {
            this.signature = signature;
            this.activity = activity;
            this.percentage = percentage;
        }
    }

    private static final class SubsetFit {
        private final double[] exposures;
        private final double rawDistance;

        private SubsetFit(double[] exposures, double rawDistance) {
            this.exposures = exposures;
            this.rawDistance = rawDistance;
        }
    }

    private static final class AddResult {
        private final double[] exposures;
        private final double distance;

        private AddResult(double[] exposures, double distance) {
            this.exposures = exposures;
            this.distance = distance;
        }
    }

    private static final class FitResult {
        private final double[] exposures;
        private final Set<Integer> activeSignatures;

        private FitResult(double[] exposures, Set<Integer> activeSignatures) {
            this.exposures = exposures;
            this.activeSignatures = activeSignatures;
        }
    }
}