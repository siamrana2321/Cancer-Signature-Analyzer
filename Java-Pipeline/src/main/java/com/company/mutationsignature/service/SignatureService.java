package com.company.mutationsignature.service;

import com.company.mutationsignature.model.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;




import java.io.*;
import java.util.*;

@Service
public class SignatureService {
    public SignatureService() throws Exception {

        loadEtiologyTable();

        System.out.println(
                "SBS Etiology Loaded Successfully"
        );
    }
    private final Map<String, String> etiologyMap =
            new HashMap<String, String>();

    private static final String GENOME_PATH =
            "D:\\Mohammad_Siam_Ahmed_Rana\\Mutation_Signature\\Resources\\genome\\GRCh38.fa";
    private static final String COSMIC_PATH =
            "D:\\Mohammad_Siam_Ahmed_Rana\\Mutation_Signature\\Resources\\cosmic\\COSMIC_v3.6_SBS_GRCh38.txt";
    private static final String SBS_Etiology_PATH =
            "D:\\Mohammad_Siam_Ahmed_Rana\\Mutation_Signature\\Resources\\cosmic\\SBS_Etiology.csv";

    private static final Map<Character, Character> COMPLEMENT =
            new LinkedHashMap<Character, Character>();

    static {

        COMPLEMENT.put('A', 'T');
        COMPLEMENT.put('T', 'A');
        COMPLEMENT.put('C', 'G');
        COMPLEMENT.put('G', 'C');

    }

    public SignatureResponse analyze(String filePath) throws Exception {

        List<Mutation> mutations = parseVcf(filePath);

        if (mutations.isEmpty()) {

            throw new RuntimeException(
                    "No SNVs found in the uploaded VCF file."
            );
        }

        Map<String, Integer> sbs96 = buildSbs96(mutations);

        SignatureResponse response = new SignatureResponse();

//        double[] tumorVector = buildTumorVector(sbs96);
        double[] tumorVector =
                buildNormalizedTumorVector(sbs96);
        CosmicMatrix cosmic = loadCosmicMatrix();

        NnlsResult fit =
                fitNnls(
                        tumorVector,
                        cosmic
                );

        List<SignatureResult> results =
                buildNnlsResults(
                        fit,
                        cosmic
                );

        response.setReconstructionCosine(
                fit.getCosine()
        );

        response.setPearson(
                fit.getPearson()
        );

//        response.setSignaturesDetected(
//                results.size()
//        );

        response.setCosmicVersion(
                "COSMIC v3.6 SBS GRCh38"
        );

        response.setRmse(
                fit.getRmse()
        );

        List<BiologicalAnnotation> annotations = buildAnnotations(results);

        response.setBiologicalAnnotations(annotations);

        response.setTotalMutations(mutations.size());

        response.setSbs96Spectrum(sbs96);

        response.setSbs96Grouped(groupSbs96(sbs96));

        response.setCosmicContributions(results);

        response.setClinicalSummary(buildClinicalSummary(results));

        response.setSbs96Percentage(buildSbs96Percentages(sbs96));

        int totalSbs = sbs96.values()
                .stream()
                .mapToInt(Integer::intValue)
                .sum();

        System.out.println("SBS96 Total = " + totalSbs);
        System.out.println(
                "Cosine = "
                        + fit.getCosine()
        );

        System.out.println(
                "Pearson = "
                        + fit.getPearson()
        );

        System.out.println(
                "RMSE = "
                        + fit.getRmse()
        );
//        CosmicMatrix cosmic = loadCosmicMatrix();

        return response;
    }

    private Map<String, Map<String, Double>> buildSbs96Percentages(
            Map<String, Integer> sbs96
    ) {

        double total =
                sbs96.values()
                        .stream()
                        .mapToInt(Integer::intValue)
                        .sum();

        if(total == 0){
            return new LinkedHashMap<>();
        }

        Map<String, Map<String, Double>> grouped =
                new LinkedHashMap<String, Map<String, Double>>();

        for (
                Map.Entry<String, Integer> entry
                : sbs96.entrySet()
        ) {

            String key =
                    entry.getKey();

            String mutation =
                    key.substring(
                            key.indexOf('[') + 1,
                            key.indexOf(']')
                    );

            char left =
                    key.charAt(0);

            char right =
                    key.charAt(
                            key.length() - 1
                    );

            char ref =
                    mutation.charAt(0);

            String context =
                    ""
                            + left
                            + ref
                            + right;

            double percentage =
                    (entry.getValue() / total)
                            * 100.0;

            if (!grouped.containsKey(mutation)) {

                grouped.put(
                        mutation,
                        new LinkedHashMap<String, Double>()
                );
            }

            grouped.get(mutation)
                    .put(
                            context,
                            percentage
                    );
        }

        return grouped;
    }

    private Map<String, Map<String, Integer>> groupSbs96(
            Map<String, Integer> sbs96
    ) {

        Map<String, Map<String, Integer>> grouped =
                new LinkedHashMap<String, Map<String, Integer>>();

        for (
                Map.Entry<String, Integer> entry
                : sbs96.entrySet()
        ) {

            String key =
                    entry.getKey();

            String mutation =
                    key.substring(
                            key.indexOf('[') + 1,
                            key.indexOf(']')
                    );

            char left =
                    key.charAt(0);

            char right =
                    key.charAt(
                            key.length() - 1
                    );

            char ref =
                    mutation.charAt(0);

            String context =
                    ""
                            + left
                            + ref
                            + right;

            if (
                    !grouped.containsKey(
                            mutation
                    )
            ) {

                grouped.put(
                        mutation,
                        new LinkedHashMap<String, Integer>()
                );
            }

            grouped
                    .get(mutation)
                    .put(
                            context,
                            entry.getValue()
                    );
        }

        return grouped;
    }

    private double calculatePearson(
            double[] a,
            double[] b
    ) {

        return new PearsonsCorrelation()
                .correlation(a, b);
    }

    private double[] buildNormalizedTumorVector(
            Map<String, Integer> sbs96
    ) {

        double[] vector = new double[96];

        double total = 0.0;

        for (Integer value : sbs96.values()) {
            total += value;
        }

        if (total == 0.0) {

            throw new RuntimeException(
                    "Empty SBS96 spectrum."
            );
        }

        int index = 0;

        for (Integer value : sbs96.values()) {

            vector[index++] =
                    value / total;
        }

        return vector;
    }

    private void loadEtiologyTable() throws Exception {

        File file =
                new File(
                        SBS_Etiology_PATH
                );

        if (!file.exists()) {

            throw new RuntimeException(
                    "SBS etiology file not found."
            );
        }

        BufferedReader br =
                new BufferedReader(
                        new FileReader(
                                SBS_Etiology_PATH
                        )
                );

        String line;

        br.readLine();

        while ((line = br.readLine()) != null) {

            String[] parts =
                    line.split(",", 2);

            if (parts.length < 2) {
                continue;
            }

            etiologyMap.put(
                    parts[0].trim(),
                    parts[1].trim()
            );
        }

        br.close();

        System.out.println(
                "Loaded SBS Etiology Entries = "
                        + etiologyMap.size()
        );
    }

    private NnlsResult fitNnls(
            double[] tumor,
            CosmicMatrix cosmic
    ) {

        int signatures =
                cosmic.getSignatureNames().size();

        double[][] A =
                cosmic.getMatrix();

        double[] weights =
                new double[signatures];

        for (int iteration = 0; iteration < 500; iteration++) {

            boolean changed = false;

            for (int j = 0; j < signatures; j++) {

                double numerator = 0.0;
                double denominator = 0.0;

                for (int row = 0; row < 96; row++) {

                    double prediction = 0.0;

                    for (int k = 0; k < signatures; k++) {

                        prediction +=
                                A[row][k]
                                        *
                                        weights[k];
                    }

                    numerator +=
                            A[row][j]
                                    *
                                    (
                                            tumor[row]
                                                    -
                                                    prediction
                                                    +
                                                    A[row][j]
                                                            *
                                                            weights[j]
                                    );

                    denominator +=
                            A[row][j]
                                    *
                                    A[row][j];
                }

                if (denominator == 0.0) {
                    continue;
                }

                double updated =
                        numerator
                                /
                                denominator;

                updated =
                        Math.max(
                                0.0,
                                updated
                        );

                if (
                        Math.abs(
                                updated
                                        -
                                        weights[j]
                        )
                                > 1e-10
                ) {
                    changed = true;
                }

                weights[j] = updated;
            }

            if (!changed) {
                break;
            }
        }

        boolean allZero = true;

        for (double w : weights) {

            if (w > 0.0) {

                allZero = false;
                break;
            }
        }

        if (allZero) {

            throw new RuntimeException(
                    "Unable to fit COSMIC signatures."
            );
        }

        double[] reconstructed =
                new double[96];

        for (int row = 0; row < 96; row++) {

            for (int col = 0;
                 col < signatures;
                 col++) {

                reconstructed[row] +=
                        A[row][col]
                                *
                                weights[col];
            }
        }

        NnlsResult result =
                new NnlsResult();

        result.setWeights(weights);

        result.setReconstructed(reconstructed);

        result.setCosine(
                cosineSimilarity(
                        tumor,
                        reconstructed
                )
        );

        result.setPearson(
                calculatePearson(
                        tumor,
                        reconstructed
                )
        );

        result.setRmse(
                calculateRmse(
                        normalize(tumor),
                        normalize(reconstructed)
                )
        );

        return result;
    }

    private double[] normalize(
            double[] vector
    ) {

        double sum = 0.0;

        for(double v : vector){
            sum += v;
        }

        if(sum == 0){
            throw new RuntimeException(
                    "Cannot normalize zero vector."
            );
        }

        double[] normalized =
                new double[vector.length];

        for(int i=0;i<vector.length;i++){

            normalized[i] =
                    vector[i] / sum;
        }

        return normalized;
    }

    private List<SignatureResult> buildNnlsResults(
            NnlsResult fit,
            CosmicMatrix cosmic
    ) {

        List<SignatureResult> results =
                new ArrayList<SignatureResult>();

        double total = 0.0;

        for (double w : fit.getWeights()) {
            total += w;
        }

        if(total == 0){
            return Collections.emptyList();
        }

        for (int i = 0;
             i < fit.getWeights().length;
             i++) {

            double contribution =
                    fit.getWeights()[i];

            if (contribution <= 0.0) {
                continue;
            }

            double percentage =
                    contribution
                            /
                            total
                            *
                            100.0;

            if (percentage < 1.0) {
                continue;
            }

            SignatureResult result =
                    new SignatureResult();

            result.setSignature(
                    cosmic
                            .getSignatureNames()
                            .get(i)
            );

            result.setSimilarity(
                    percentage
            );

            result.setEtiology(
                    getSignatureMeaning(
                            cosmic
                                    .getSignatureNames()
                                    .get(i)
                    )
            );

            results.add(result);
        }

        results.sort(
                (a, b) ->
                        Double.compare(
                                b.getSimilarity(),
                                a.getSimilarity()
                        )
        );

        if (results.size() > 3) {

            return new ArrayList<>(
                    results.subList(0, 3)
            );
        }

        return results;
    }

    private String getSignatureMeaning(
            String signature
    ) {

        String meaning =
                etiologyMap.get(signature);

        if (meaning == null) {
            return "Unknown";
        }

        return meaning;
    }

    private List<BiologicalAnnotation> buildAnnotations(
            List<SignatureResult> signatures
    ) {

        List<BiologicalAnnotation> annotations =
                new ArrayList<BiologicalAnnotation>();

        for (SignatureResult result : signatures) {

            BiologicalAnnotation annotation =
                    new BiologicalAnnotation();

            annotation.setSignature(
                    result.getSignature()
            );

            annotation.setEtiology(
                    getSignatureMeaning(
                            result.getSignature()
                    )
            );

            if(result.getSimilarity() >= 20)
            {
                annotation.setConfidence("High");
            }
            else if(result.getSimilarity() >= 5)
            {
                annotation.setConfidence("Moderate");
            }
            else
            {
                annotation.setConfidence("Low");
            }

            annotations.add(
                    annotation
            );
        }

        return annotations;
    }

    private String buildClinicalSummary(
            List<SignatureResult> signatures
    ) {

        StringBuilder summary =
                new StringBuilder();

        for (SignatureResult result : signatures) {

            if (result.getSimilarity() < 5.0) {
                continue;
            }

            summary.append(
                    result.getSignature()
            );

            summary.append(" (");

            summary.append(
                    String.format(
                            "%.1f",
                            result.getSimilarity()
                    )
            );

            summary.append("%) : ");

            summary.append(
                    result.getEtiology()
            );

            summary.append(". ");
        }

        return summary.toString();
    }

    private List<SignatureResult> buildCosmicResults(
            double[] contributions,
            List<String> signatureNames
    ) {

        List<SignatureResult> results =
                new ArrayList<SignatureResult>();

        double total = 0.0;

        for (double value : contributions) {
            total += value;
        }

        for (int i = 0; i < contributions.length; i++) {

            if (contributions[i] <= 0.0) {
                continue;
            }

            double percentage =
                    (contributions[i] / total) * 100.0;

            if (percentage < 1.0) {
                continue;
            }

            SignatureResult result =
                    new SignatureResult();

            result.setSignature(
                    signatureNames.get(i)
            );

            result.setSimilarity(
                    percentage
            );

            results.add(result);
        }

        results.sort(
                (a, b) ->
                        Double.compare(
                                b.getSimilarity(),
                                a.getSimilarity()
                        )
        );

        return results;
    }



    private double cosineSimilarity(
            double[] a,
            double[] b
    ) {

        double dot = 0.0;
        double normA = 0.0;
        double normB = 0.0;

        for (int i = 0; i < a.length; i++) {

            dot += a[i] * b[i];

            normA += a[i] * a[i];

            normB += b[i] * b[i];
        }

        if (normA == 0 || normB == 0) {
            return 0.0;
        }

        return dot /
                (
                        Math.sqrt(normA)
                                *
                                Math.sqrt(normB)
                );
    }


    private double[] buildTumorVector(
            Map<String, Integer> sbs96
    ) {

        double[] vector =
                new double[96];

        int index = 0;

        for (Integer value : sbs96.values()) {

            vector[index++] = value;
        }

        return vector;
    }

    private double calculateRmse(
            double[] original,
            double[] reconstructed
    ) {

        double sum = 0.0;

        for (int i = 0; i < original.length; i++) {

            double diff =
                    original[i]
                            - reconstructed[i];

            sum += diff * diff;
        }

        return Math.sqrt(
                sum / original.length
        );
    }

    private CosmicMatrix loadCosmicMatrix()
            throws Exception {
        List<String> contexts =
                new ArrayList<String>();

        File file =
                new File(COSMIC_PATH);

        if (!file.exists()) {

            throw new RuntimeException(
                    "COSMIC matrix file not found."
            );
        }

        BufferedReader reader =
                new BufferedReader(
                        new FileReader(COSMIC_PATH)
                );

        String header =
                reader.readLine();

        String[] columns =
                header.split("\t");

        List<String> signatureNames =
                new ArrayList<String>();

        for (int i = 1; i < columns.length; i++) {

            signatureNames.add(columns[i]);
        }

        int signatureCount =
                signatureNames.size();

        double[][] matrix =
                new double[96][signatureCount];

        String line;

        int row = 0;

        while (
                (line = reader.readLine()) != null
        ) {
            if(row >= 96){
                break;
            }

            String[] tokens =
                    line.split("\t");
            contexts.add(tokens[0]);

            for (
                    int col = 1;
                    col < tokens.length;
                    col++
            ) {

                matrix[row][col - 1] =
                        Double.parseDouble(
                                tokens[col]
                        );
            }

            row++;
        }

        reader.close();

        System.out.println(
                "COSMIC rows = " + row
        );

        System.out.println(
                "Signatures = "
                        + signatureCount
        );

        return new CosmicMatrix(
                matrix,
                signatureNames,
                contexts
        );
    }

    private List<Mutation> parseVcf(
            String filePath
    ) throws Exception {

        List<Mutation> mutations =
                new ArrayList<Mutation>();

        BufferedReader reader =
                new BufferedReader(
                        new FileReader(
                                filePath
                        )
                );

        String line;

        while ((line = reader.readLine()) != null) {

            if (line.startsWith("#")) {
                continue;
            }

            String[] cols = line.split("\t");

            if (cols.length < 5) {
                continue;
            }

            String chromosome =
                    cols[0];

            long position =
                    Long.parseLong(cols[1]);

            String ref =
                    cols[3];

            String alt =
                    cols[4];

            if (ref.length() == 1 &&
                    alt.length() == 1) {

                mutations.add(
                        new Mutation(
                                chromosome,
                                position,
                                ref,
                                alt
                        )
                );
            }
        }

        reader.close();

        return mutations;
    }

    private Map<String, Integer> buildSbs96(
            List<Mutation> mutations
    ) throws Exception {

        Map<String, Integer> sbs96 = initializeSbs96();

        IndexedFastaSequenceFile genome =
                new IndexedFastaSequenceFile(
                        new File(GENOME_PATH)
                );

        for (Mutation mutation : mutations) {

            String chromosome =
                    mutation.getChromosome();

            long position =
                    mutation.getPosition();

//            ReferenceSequence referenceSequence =
//                    genome.getSubsequenceAt(
//                            chromosome,
//                            position - 1,
//                            position + 1
//                    );

            String trinucleotide;
            try {

                ReferenceSequence referenceSequence =
                        genome.getSubsequenceAt(
                                chromosome,
                                position - 1,
                                position + 1
                        );

                trinucleotide = new String(referenceSequence.getBases()).toUpperCase();

                if (trinucleotide.length() != 3) {
                    continue;
                }

            }
            catch (Exception e) {

                continue;
            }

            char left = trinucleotide.charAt(0);

            char ref =
                    mutation.getRef()
                            .toUpperCase()
                            .charAt(0);

            char alt =
                    mutation.getAlt()
                            .toUpperCase()
                            .charAt(0);

            char right = trinucleotide.charAt(2);

            String context =
                    normalizeContext(
                            left,
                            ref,
                            alt,
                            right
                    );

            if (sbs96.containsKey(context)) {

                sbs96.put(
                        context,
                        sbs96.get(context) + 1
                );

            }
        }

        genome.close();

        return sbs96;
    }

    private Map<String, Integer> initializeSbs96() {

        Map<String, Integer> sbs96 =
                new LinkedHashMap<String, Integer>();

        char[] bases =
                {'A', 'C', 'G', 'T'};

        String[] substitutions = {

                "C>A",
                "C>G",
                "C>T",
                "T>A",
                "T>C",
                "T>G"
        };

        for (char left : bases) {

            for (String substitution : substitutions) {

                for (char right : bases) {

                    String context =
                            left +
                                    "[" +
                                    substitution +
                                    "]" +
                                    right;

                    sbs96.put(context, 0);
                }
            }
        }

        return sbs96;
    }

    private String normalizeContext(
            char left,
            char ref,
            char alt,
            char right
    ) {

        if (ref == 'C' || ref == 'T') {

            return left +
                    "[" +
                    ref +
                    ">" +
                    alt +
                    "]" +
                    right;
        }

        String reverseComplement =
                reverseComplement(
                        "" +
                                left +
                                ref +
                                right
                );

        char normalizedLeft =
                reverseComplement.charAt(0);

        char normalizedRight =
                reverseComplement.charAt(2);

        char normalizedRef =
                COMPLEMENT.get(ref);

        char normalizedAlt =
                COMPLEMENT.get(alt);

        return normalizedLeft +
                "[" +
                normalizedRef +
                ">" +
                normalizedAlt +
                "]" +
                normalizedRight;
    }

    private String reverseComplement(
            String sequence
    ) {

        StringBuilder builder =
                new StringBuilder();

        for (
                int i = sequence.length() - 1;
                i >= 0;
                i--
        ) {

            builder.append(
                    COMPLEMENT.get(
                            sequence.charAt(i)
                    )
            );
        }

        return builder.toString();
    }
}