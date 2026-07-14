package com.company.mutationsignature.service;

import com.company.mutationsignature.model.SBS96Response;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.springframework.stereotype.Service;

import java.io.*;
import java.util.*;

@Service
public class SBS96SignatureService {

    private static final String GENOME_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\genome\\GRCh38.fa";
    private static final String COSMIC_PATH =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\cosmic\\COSMIC_v3.4_SBS_GRCh38_exome.txt";

    // SigProfilerAssignment defaults used by spa_analyze/decompose_fit.
    private static final double ADD_PENALTY = 0.05;
    private static final double REMOVE_PENALTY = 0.01;
    private static final double EPS = 1e-12;

    private static final Map<Character, Character> COMPLEMENT =
            new HashMap<Character, Character>();

    static {
        COMPLEMENT.put('A', 'T');
        COMPLEMENT.put('T', 'A');
        COMPLEMENT.put('C', 'G');
        COMPLEMENT.put('G', 'C');
    }

    private static final Map<String, String> ETIOLOGY_MAP =
            new HashMap<String, String>();

    static {
        ETIOLOGY_MAP.put("SBS1", "Spontaneous deamination of 5-methylcytosine (clock-like signature)");
        ETIOLOGY_MAP.put("SBS2", "Activity of APOBEC family of cytidine deaminases");
        ETIOLOGY_MAP.put("SBS3", "Defective homologous recombination DNA damage repair");
        ETIOLOGY_MAP.put("SBS4", "Tobacco smoking");
        ETIOLOGY_MAP.put("SBS5", "Unknown (clock-like signature)");
        ETIOLOGY_MAP.put("SBS6", "Defective DNA mismatch repair");

        ETIOLOGY_MAP.put("SBS7a", "Ultraviolet light exposure");
        ETIOLOGY_MAP.put("SBS7b", "Ultraviolet light exposure");
        ETIOLOGY_MAP.put("SBS7c", "Ultraviolet light exposure");
        ETIOLOGY_MAP.put("SBS7d", "Ultraviolet light exposure");

        ETIOLOGY_MAP.put("SBS8", "Unknown");
        ETIOLOGY_MAP.put("SBS9", "Polymerase eta somatic hypermutation");

        ETIOLOGY_MAP.put("SBS10a", "Polymerase epsilon exonuclease domain mutations");
        ETIOLOGY_MAP.put("SBS10b", "Polymerase epsilon exonuclease domain mutations");
        ETIOLOGY_MAP.put("SBS10c", "Polymerase epsilon mutations");
        ETIOLOGY_MAP.put("SBS10d", "Polymerase delta mutations");

        ETIOLOGY_MAP.put("SBS11", "Temozolomide treatment");
        ETIOLOGY_MAP.put("SBS12", "Unknown");
        ETIOLOGY_MAP.put("SBS13", "Activity of APOBEC family of cytidine deaminases");
        ETIOLOGY_MAP.put("SBS14", "Defective DNA mismatch repair and polymerase epsilon mutations");
        ETIOLOGY_MAP.put("SBS15", "Defective DNA mismatch repair");
        ETIOLOGY_MAP.put("SBS16", "Unknown");

        ETIOLOGY_MAP.put("SBS17a", "Unknown");
        ETIOLOGY_MAP.put("SBS17b", "Unknown");

        ETIOLOGY_MAP.put("SBS18", "Reactive oxygen species");
        ETIOLOGY_MAP.put("SBS19", "Unknown");
        ETIOLOGY_MAP.put("SBS20", "Defective DNA mismatch repair and polymerase delta mutations");
        ETIOLOGY_MAP.put("SBS21", "Defective DNA mismatch repair");

        ETIOLOGY_MAP.put("SBS22a", "Aristolochic acid exposure");
        ETIOLOGY_MAP.put("SBS22b", "Aristolochic acid exposure");

        ETIOLOGY_MAP.put("SBS23", "Unknown");
        ETIOLOGY_MAP.put("SBS24", "Aflatoxin exposure");
        ETIOLOGY_MAP.put("SBS25", "Chemotherapy treatment");
        ETIOLOGY_MAP.put("SBS26", "Defective DNA mismatch repair");
        ETIOLOGY_MAP.put("SBS27", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS28", "Unknown");
        ETIOLOGY_MAP.put("SBS29", "Tobacco chewing");
        ETIOLOGY_MAP.put("SBS30", "Defective base excision repair due to NTHL1 mutations");
        ETIOLOGY_MAP.put("SBS31", "Platinum chemotherapy treatment");
        ETIOLOGY_MAP.put("SBS32", "Azathioprine treatment");
        ETIOLOGY_MAP.put("SBS33", "Unknown");
        ETIOLOGY_MAP.put("SBS34", "Unknown");
        ETIOLOGY_MAP.put("SBS35", "Platinum chemotherapy treatment");
        ETIOLOGY_MAP.put("SBS36", "Defective base excision repair due to MUTYH mutations");
        ETIOLOGY_MAP.put("SBS37", "Unknown");
        ETIOLOGY_MAP.put("SBS38", "Ultraviolet light exposure");
        ETIOLOGY_MAP.put("SBS39", "Unknown");

        ETIOLOGY_MAP.put("SBS40a", "Unknown");
        ETIOLOGY_MAP.put("SBS40b", "Unknown");
        ETIOLOGY_MAP.put("SBS40c", "Unknown");

        ETIOLOGY_MAP.put("SBS41", "Unknown");
        ETIOLOGY_MAP.put("SBS42", "Haloalkane exposure");
        ETIOLOGY_MAP.put("SBS43", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS44", "Defective DNA mismatch repair");
        ETIOLOGY_MAP.put("SBS45", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS46", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS47", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS48", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS49", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS50", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS51", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS52", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS53", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS54", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS55", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS56", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS57", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS58", "Unknown");
        ETIOLOGY_MAP.put("SBS59", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS60", "Possible sequencing artefact");

        ETIOLOGY_MAP.put("SBS84", "Activation-induced cytidine deaminase activity");
        ETIOLOGY_MAP.put("SBS85", "Activation-induced cytidine deaminase activity");
        ETIOLOGY_MAP.put("SBS86", "Unknown");
        ETIOLOGY_MAP.put("SBS87", "Thiopurine chemotherapy treatment");
        ETIOLOGY_MAP.put("SBS88", "Colibactin exposure");
        ETIOLOGY_MAP.put("SBS89", "Unknown");
        ETIOLOGY_MAP.put("SBS90", "Polymerase epsilon exonuclease domain mutations");
        ETIOLOGY_MAP.put("SBS91", "Unknown");
        ETIOLOGY_MAP.put("SBS92", "Tobacco smoking");
        ETIOLOGY_MAP.put("SBS93", "Unknown");
        ETIOLOGY_MAP.put("SBS94", "Unknown");
        ETIOLOGY_MAP.put("SBS95", "Possible sequencing artefact");
        ETIOLOGY_MAP.put("SBS96", "Unknown");
        ETIOLOGY_MAP.put("SBS97", "Unknown");
        ETIOLOGY_MAP.put("SBS98", "Unknown");
        ETIOLOGY_MAP.put("SBS99", "Unknown");
        ETIOLOGY_MAP.put("SBS100", "Unknown");
    }

    public SBS96SignatureService() {
    }

    public SBS96Response analyze(String filePath) throws Exception {

        CosmicMatrix cosmic = loadCosmicMatrix();

        // IMPORTANT: SigProfilerMatrixGenerator VCF conversion behaviour:
        // 1-bp REF/ALT -> one SBS.
        // 2-bp REF/ALT -> split into two SBS only when BOTH bases change.
        // everything else -> INDEL/other stream, not SBS96.
        List<Mutation> mutations = parseVcfLikeSigProfiler(filePath);

        if (mutations.isEmpty()) {
            throw new RuntimeException("No SBS-compatible mutations found in VCF.");
        }

        Map<String, Integer> sbs96 = buildSbs96(mutations, cosmic.getContexts());
        double[] sample = buildTumorVector(sbs96, cosmic.getContexts());

        int totalMutations = (int) Math.round(sum(sample));

        AssignmentResult fit = addRemoveSignatures(
                cosmic.getMatrix(),
                sample,
                cosmic.getSignatureNames()
        );

        List<SBS96Response.SignatureContribution> signatureContributions =
                buildSignatureContributions(fit.getActivities(), cosmic.getSignatureNames());

        SBS96Response response = new SBS96Response();
        response.setCosmicVersion("COSMIC v3.4 SBS GRCh38 WES (exome-normalized)");
        response.setTotalMutations(totalMutations);
        response.setSbs96Spectrum(sbs96);
        response.setSbs96Grouped(groupSbs96(sbs96));
        response.setSbs96Percentage(buildSbs96Percentages(sbs96));
        response.setSignatureContributions(signatureContributions);
        response.setClinicalSummary(buildClinicalSummary(signatureContributions));

        response.setReconstructionCosine(fit.getCosine());
        response.setPearson(fit.getPearson());
        response.setRmse(fit.getRmse());
        response.setL1Norm(fit.getL1Norm());
        response.setL1NormPercent(fit.getL1NormPercent());
        response.setL2Norm(fit.getL2Norm());
        response.setL2NormPercent(fit.getL2NormPercent());
        response.setKlDivergence(fit.getKlDivergence());

        System.out.println("SBS96 Total = " + totalMutations);
        System.out.println("Selected signatures = " + fit.getSelectedSignatureNames());
        System.out.println("Cosine = " + fit.getCosine());
        System.out.println("L1 = " + fit.getL1Norm() + " (" + fit.getL1NormPercent() + "%)");
        System.out.println("L2 = " + fit.getL2Norm() + " (" + fit.getL2NormPercent() + "%)");
        System.out.println("KL = " + fit.getKlDivergence());
        System.out.println("Pearson = " + fit.getPearson());

        return response;
    }

    // -------------------------------------------------------------------------
    // SigProfilerMatrixGenerator-compatible VCF conversion
    // -------------------------------------------------------------------------

    private List<Mutation> parseVcfLikeSigProfiler(String filePath) throws Exception {

        List<Mutation> mutations = new ArrayList<Mutation>();
        String previousVariantLine = null;

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String rawLine;

            while ((rawLine = reader.readLine()) != null) {

                if (rawLine.length() == 0 || rawLine.charAt(0) == '#') {
                    continue;
                }

                String[] cols = rawLine.trim().split("\\s+");
                if (cols.length < 5) {
                    continue;
                }

                String chromosome = normalizeChromosome(cols[0]);
                long position;

                try {
                    position = Long.parseLong(cols[1]);
                } catch (NumberFormatException e) {
                    continue;
                }

                String ref = cols[3].toUpperCase(Locale.ROOT);
                String alt = cols[4].toUpperCase(Locale.ROOT);

                // SigProfiler compares the parsed full VCF row with the previous row.
                String comparableLine = rawLine.trim();

                if (comparableLine.equals(previousVariantLine)) {
                    previousVariantLine = comparableLine;
                    continue;
                }

                if (ref.length() == 1 && alt.length() == 1) {

                    if (!isCanonicalBase(ref.charAt(0))
                            || !isCanonicalBase(alt.charAt(0))
                            || ref.equals(alt)) {
                        previousVariantLine = comparableLine;
                        continue;
                    }

                    mutations.add(new Mutation(
                            chromosome, position, ref.charAt(0), alt.charAt(0)
                    ));

                } else if (ref.length() == 2
                        && alt.length() == 2
                        && ref.indexOf('-') < 0
                        && alt.indexOf('-') < 0) {

                    char ref1 = ref.charAt(0);
                    char ref2 = ref.charAt(1);
                    char alt1 = alt.charAt(0);
                    char alt2 = alt.charAt(1);

                    // This intentionally mirrors convertVCF(): if either base is
                    // invalid or unchanged, the whole 2-bp record is skipped.
                    if (!isCanonicalBase(ref1)
                            || !isCanonicalBase(ref2)
                            || !isCanonicalBase(alt1)
                            || !isCanonicalBase(alt2)
                            || ref1 == alt1
                            || ref2 == alt2) {
                        previousVariantLine = comparableLine;
                        continue;
                    }

                    mutations.add(new Mutation(chromosome, position, ref1, alt1));
                    mutations.add(new Mutation(chromosome, position + 1L, ref2, alt2));
                }

                previousVariantLine = comparableLine;
            }
        }

        return mutations;
    }

    private String normalizeChromosome(String chromosome) {
        // SigProfiler's VCF converter removes the first 3 characters when
        // chromosome text length > 2 (e.g. chr1 -> 1, chrX -> X).
        String normalized = chromosome;
        if (normalized.length() > 2) {
            normalized = normalized.substring(3);
        }
        if ("M".equalsIgnoreCase(normalized) || "mt".equalsIgnoreCase(normalized)) {
            normalized = "MT";
        }
        return normalized;
    }

    private boolean isCanonicalBase(char base) {
        return base == 'A' || base == 'C' || base == 'G' || base == 'T';
    }

    // -------------------------------------------------------------------------
    // SBS96 matrix generation
    // -------------------------------------------------------------------------

    private Map<String, Integer> buildSbs96(
            List<Mutation> mutations,
            List<String> cosmicContexts) throws Exception {

        Map<String, Integer> spectrum = new LinkedHashMap<String, Integer>();

        // Use COSMIC row order as the authoritative SBS96 order.
        for (String context : cosmicContexts) {
            spectrum.put(context, 0);
        }

        try (IndexedFastaSequenceFile genome =
                     new IndexedFastaSequenceFile(new File(GENOME_PATH))) {

            for (Mutation mutation : mutations) {

                if (mutation.position <= 1L) {
                    continue;
                }

                ReferenceSequence sequence;
                try {
                    sequence = getReferenceSequence(genome, mutation);
                } catch (Exception e) {
                    continue;
                }

                String tri = new String(sequence.getBases()).toUpperCase(Locale.ROOT);

                if (tri.length() != 3
                        || !isCanonicalBase(tri.charAt(0))
                        || !isCanonicalBase(tri.charAt(1))
                        || !isCanonicalBase(tri.charAt(2))) {
                    continue;
                }

                // Reference validation is essential. Do not classify a VCF REF
                // against a different FASTA reference base.
                if (tri.charAt(1) != mutation.ref) {
                    continue;
                }

                String context = normalizeContext(
                        tri.charAt(0),
                        mutation.ref,
                        mutation.alt,
                        tri.charAt(2)
                );

                Integer old = spectrum.get(context);
                if (old != null) {
                    spectrum.put(context, old + 1);
                }
            }
        }

        return spectrum;
    }

    private ReferenceSequence getReferenceSequence(
            IndexedFastaSequenceFile genome,
            Mutation mutation) {

        try {
            return genome.getSubsequenceAt(
                    mutation.chromosome,
                    mutation.position - 1L,
                    mutation.position + 1L
            );
        } catch (Exception first) {
            // Allows FASTA dictionaries using chr1 while SigProfiler-style
            // internal chromosome names are 1.
            String alternate = mutation.chromosome.startsWith("chr")
                    ? mutation.chromosome.substring(3)
                    : "chr" + mutation.chromosome;

            return genome.getSubsequenceAt(
                    alternate,
                    mutation.position - 1L,
                    mutation.position + 1L
            );
        }
    }

    private String normalizeContext(char left, char ref, char alt, char right) {

        if (ref == 'C' || ref == 'T') {
            return "" + left + "[" + ref + ">" + alt + "]" + right;
        }

        char newLeft = complement(right);
        char newRef = complement(ref);
        char newAlt = complement(alt);
        char newRight = complement(left);

        return "" + newLeft + "[" + newRef + ">" + newAlt + "]" + newRight;
    }

    private char complement(char base) {
        Character value = COMPLEMENT.get(base);
        if (value == null) {
            throw new IllegalArgumentException("Unsupported DNA base: " + base);
        }
        return value;
    }

    private double[] buildTumorVector(
            Map<String, Integer> spectrum,
            List<String> cosmicContexts) {

        double[] sample = new double[cosmicContexts.size()];

        for (int i = 0; i < cosmicContexts.size(); i++) {
            Integer count = spectrum.get(cosmicContexts.get(i));
            sample[i] = count == null ? 0.0 : count.doubleValue();
        }

        return sample;
    }

    // -------------------------------------------------------------------------
    // SigProfilerAssignment-style forward stagewise add/remove fitting
    // -------------------------------------------------------------------------

    private AssignmentResult addRemoveSignatures(
            double[][] W,
            double[] sample,
            List<String> signatureNames) {

        /*
         * SigProfilerAssignment v1.1.4 COSMIC-fit control flow:
         *
         * fit_signatures(all COSMIC signatures)
         *     -> remove_all_single_signatures(cutoff = initial_remove_penalty = 0.05)
         *     -> union(non-zero signatures, SBS1/SBS5)
         *     -> add_remove_signatures(add = 0.05, remove = 0.01)
         *     -> roundConserveSum
         *
         * This is the flow used by make_final_solution/process_sample.
         */

        List<Integer> backgroundSignatures = new ArrayList<Integer>();

        int sbs1 = signatureNames.indexOf("SBS1");
        int sbs5 = signatureNames.indexOf("SBS5");

        if (sbs1 >= 0) {
            backgroundSignatures.add(sbs1);
        }

        if (sbs5 >= 0) {
            backgroundSignatures.add(sbs5);
        }

        // ---------------- Initial fit: ss.fit_signatures(processAvg, genome) ----------------
        double[] initialExposure = fitSignaturesExact(W, sample);

        // ---------------- Initial remove: cutoff = initial_remove_penalty ----------------
        RemoveResult initialRemove = removeAllSingleSignaturesExact(
                W,
                initialExposure,
                sample,
                Collections.<Integer>emptyList(),
                0.05
        );

        double[] exposure = initialRemove.exposure;

        // init_add_sig_idx = union(nonzero(exposure), background_sigs)
        List<Integer> backgroundSigs =
                nonZeroIndicesExact(exposure);

        backgroundSigs =
                unionIndices(
                        backgroundSigs,
                        backgroundSignatures
                );

        // ---------------- ss.add_remove_signatures ----------------
        AddRemoveResult sparse = addRemoveExact(
                W,
                sample,
                backgroundSigs,
                backgroundSignatures,
                signatureNames,
                ADD_PENALTY,
                REMOVE_PENALTY
        );

        double[] finalActivities = roundConserveSumExact(sparse.exposure);

        if (Math.round(sum(finalActivities)) != Math.round(sum(sample))) {
            throw new IllegalStateException(
                    "Mutation count not conserved. Expected "
                            + Math.round(sum(sample))
                            + " but assigned "
                            + Math.round(sum(finalActivities))
            );
        }

        double[] reconstructed = multiply(W, finalActivities);

        AssignmentResult result = new AssignmentResult();
        result.activities = finalActivities;
        result.reconstructed = reconstructed;
        result.cosine = cosineSimilarity(sample, reconstructed);
        result.pearson = pearson(sample, reconstructed);
        result.rmse = rmse(normalize(sample), normalize(reconstructed));
        result.l1Norm = l1Norm(sample, reconstructed);
        result.l1NormPercent =
                result.l1Norm / l1Norm(sample, new double[sample.length]) * 100.0;
        result.l2Norm = l2Norm(sample, reconstructed);
        result.l2NormPercent =
                result.l2Norm / l2Norm(sample, new double[sample.length]) * 100.0;
        result.klDivergence = klDivergence(sample, reconstructed);

        result.selectedSignatureNames = new ArrayList<String>();
        for (int i = 0; i < finalActivities.length; i++) {
            if (finalActivities[i] != 0.0) {
                result.selectedSignatureNames.add(signatureNames.get(i));
            }
        }

        return result;
    }

    /**
     * Exact control flow of single_sample.fit_signatures(..., metric="l2").
     * Distance is calculated from raw NNLS weights. Returned activities are
     * normalized to the sample total and rounded with roundConserveSum.
     */
    private double[] fitSignaturesExact(double[][] W, double[] genome) {

        double[] weights = nnls(W, genome);
        double weightSum = sum(weights);

        if (weightSum == 0.0) {
            throw new IllegalStateException("NNLS returned zero total weight.");
        }

        double[] solution = new double[weights.length];
        double mutationTotal = sum(genome);

        for (int i = 0; i < weights.length; i++) {
            solution[i] = weights[i] / weightSum * mutationTotal;
        }

        double[] exposure = roundConserveSumExact(solution);

        if (Math.round(sum(exposure)) != Math.round(mutationTotal)) {
            throw new IllegalStateException("Mutation count not conserved in fit_signatures.");
        }

        return exposure;
    }

    /**
     * Direct Java translation of SigProfilerAssignment v1.1.4
     * single_sample.add_remove_signatures.
     */
    private AddRemoveResult addRemoveExact(
            double[][] W,
            double[] sample,
            List<Integer> initialBackground,
            List<Integer> permanentSigs,
            List<String> allSigIds,
            double addPenalty,
            double removePenalty) {

        List<Integer> backgroundSigs =
                new ArrayList<Integer>(initialBackground);
        List<Integer> alwaysBackground =
                new ArrayList<Integer>(permanentSigs);

        List<Integer> candidateSigs = new ArrayList<Integer>();
        for (int i = 0; i < W[0].length; i++) {
            candidateSigs.add(i);
        }

        double originalDistance = Double.POSITIVE_INFINITY;
        double[] finalActivities = null;
        List<Integer> finalBackground = null;

        while (true) {

            double layerOriginalDistance = Double.POSITIVE_INFINITY;
            double[] activities = null;
            List<Integer> selectedSignatures = null;

            List<Integer> sigsToBeAdded =
                    differenceIndices(candidateSigs, backgroundSigs);

            for (Integer candidate : sigsToBeAdded) {

                List<Integer> loopSig =
                        Collections.singletonList(candidate);

                backgroundSigs =
                        unionIndices(backgroundSigs, alwaysBackground);

                backgroundSigs =
                        addConnectedSigsExact(backgroundSigs, allSigIds);

                AddResult add = addSignaturesExact(
                        W,
                        sample,
                        new ArrayList<Integer>(backgroundSigs),
                        loopSig,
                        addPenalty
                );

                RemoveResult remove = removeAllSingleSignaturesExact(
                        W,
                        add.exposure,
                        sample,
                        alwaysBackground,
                        removePenalty
                );

                /*
                 * Preserve the upstream condition literally:
                 *
                 * np.nonzero(add_exposures)[0].all()
                 *     == np.nonzero(remove_exposures)[0].all()
                 * and shapes equal
                 *
                 * np.ndarray.all() returns a boolean, not array equality.
                 */
                int[] addNonZero = nonZeroIndexArrayExact(add.exposure);
                int[] removeNonZero = nonZeroIndexArrayExact(remove.exposure);

                boolean sameAccordingToUpstream =
                        numpyAll(addNonZero) == numpyAll(removeNonZero)
                                && addNonZero.length == removeNonZero.length;

                double distance;
                double[] exposure;

                if (sameAccordingToUpstream) {
                    distance = add.distance;
                    exposure = add.exposure;
                } else {
                    distance = remove.distance;
                    exposure = remove.exposure;
                }

                if (distance < layerOriginalDistance) {
                    selectedSignatures = nonZeroIndicesExact(exposure);
                    layerOriginalDistance = distance;
                    activities = exposure;
                }
            }

            if (activities == null) {
                throw new IllegalStateException(
                        "SigProfiler add/remove layer produced no candidate solution."
                );
            }

            if (layerOriginalDistance < originalDistance) {
                originalDistance = layerOriginalDistance;
                backgroundSigs =
                        new ArrayList<Integer>(selectedSignatures);
                finalBackground =
                        new ArrayList<Integer>(selectedSignatures);
                finalActivities =
                        Arrays.copyOf(activities, activities.length);
            } else {
                break;
            }
        }

        if (finalActivities == null) {
            throw new IllegalStateException("No final SigProfiler assignment solution.");
        }

        AddRemoveResult result = new AddRemoveResult();
        result.background = finalBackground;
        result.exposure = roundConserveSumExact(finalActivities);
        result.distance = originalDistance;
        return result;
    }

    /**
     * Direct Java translation of single_sample.add_signatures using NNLS/L2.
     */
    private AddResult addSignaturesExact(
            double[][] W,
            double[] genome,
            List<Integer> presentSignatures,
            List<Integer> toBeAdded,
            double cutoff) {

        List<Integer> notToBeAdded = new ArrayList<Integer>();

        for (int i = 0; i < W[0].length; i++) {
            if (!toBeAdded.contains(i)) {
                notToBeAdded.add(i);
            }
        }

        List<Integer> listed = presentSignatures;
        List<Integer> nonListed = new ArrayList<Integer>();

        for (int i = 0; i < W[0].length; i++) {
            if (!listed.contains(i) && !notToBeAdded.contains(i)) {
                nonListed.add(i);
            }
        }

        double originalSimilarity = Double.POSITIVE_INFINITY;
        double[] finalExposure = null;

        if (!listed.isEmpty()) {
            Collections.sort(listed);

            NnlsFit fit = fitSubsetRaw(W, genome, listed);
            originalSimilarity = relativeL2(genome, fit.rawReconstruction);

            // Upstream add_signatures uses np.round + maximum correction.
            finalExposure = expandExposure(
                    numpyRoundAndCorrectMaximum(
                            fit.normalizedExposure,
                            Math.round(sum(genome))
                    ),
                    listed,
                    W[0].length
            );
        }

        while (true) {

            double bestDifference = Double.NEGATIVE_INFINITY;
            double bestSimilarity = Double.POSITIVE_INFINITY;
            double[] loopExposure = null;
            int bestSignature = -1;

            for (Integer sig : new ArrayList<Integer>(nonListed)) {

                List<Integer> loopListed =
                        new ArrayList<Integer>(listed);

                if (!listed.isEmpty()) {
                    loopListed.add(sig);
                    Collections.sort(loopListed);
                } else {
                    loopListed.add(sig);
                }

                NnlsFit fit = fitSubsetRaw(W, genome, loopListed);
                double newSimilarity =
                        relativeL2(genome, fit.rawReconstruction);

                double[] newExposure = expandExposure(
                        numpyRoundAndCorrectMaximum(
                                fit.normalizedExposure,
                                Math.round(sum(genome))
                        ),
                        loopListed,
                        W[0].length
                );

                double difference =
                        originalSimilarity - newSimilarity;

                if (difference > bestDifference) {
                    bestDifference = difference;
                    bestSimilarity = newSimilarity;
                    loopExposure = newExposure;
                    bestSignature = sig;
                }
            }

            if (bestSignature >= 0
                    && originalSimilarity - bestSimilarity > cutoff) {

                originalSimilarity = bestSimilarity;
                listed.add(bestSignature);
                nonListed.remove(Integer.valueOf(bestSignature));
                Collections.sort(listed);
                finalExposure = loopExposure;

                if (!nonListed.isEmpty()) {
                    continue;
                }
            }

            break;
        }

        if (finalExposure == null) {
            finalExposure = new double[W[0].length];
        }

        AddResult result = new AddResult();
        result.exposure = finalExposure;
        result.distance = originalSimilarity;
        return result;
    }

    /**
     * Direct Java translation of single_sample.remove_all_single_signatures.
     */
    private RemoveResult removeAllSingleSignaturesExact(
            double[][] W,
            double[] H,
            double[] genomes,
            List<Integer> backgroundSigs,
            double cutoff) {

        List<Integer> backgroundSig =
                new ArrayList<Integer>(backgroundSigs);

        List<Integer> baseHIndex = nonZeroIndicesExact(H);

        if (baseHIndex.isEmpty()) {
            RemoveResult result = new RemoveResult();
            result.exposure = new double[H.length];
            result.distance = Double.POSITIVE_INFINITY;
            return result;
        }

        NnlsFit baseFit = fitSubsetRaw(W, genomes, baseHIndex);
        double originalSimilarity =
                relativeL2(genomes, baseFit.rawReconstruction);

        double[] oldExposures = roundConserveSumExact(H);

        if (nonZeroIndicesExact(oldExposures).size() <= 1) {
            RemoveResult result = new RemoveResult();
            result.exposure = oldExposures;
            result.distance = originalSimilarity;
            return result;
        }

        double[] successExposure = null;
        double successDistance = 0.0;

        while (true) {

            double[] activeExposure =
                    successExposure == null
                            ? oldExposures
                            : successExposure;

            List<Integer> initialZerosIdx =
                    zeroIndicesExact(activeExposure);

            List<Integer> selectableIdx =
                    positiveIndicesExact(activeExposure);

            if (selectableIdx.size() <= 1) {
                break;
            }

            backgroundSig = getChangedBackgroundSigIdxExact(
                    oldExposures,
                    backgroundSig
            );

            double bestDifference = Double.POSITIVE_INFINITY;
            double[] bestExposure = null;
            double bestSimilarity = 1.0;

            double[][] Winit =
                    subsetColumns(W, selectableIdx);

            int l = Winit[0].length;

            for (int localIndex = 0; localIndex < l; localIndex++) {

                if (backgroundSig.contains(localIndex)) {
                    continue;
                }

                List<Integer> loopSelection = new ArrayList<Integer>();
                for (int i = 0; i < l; i++) {
                    if (i != localIndex) {
                        loopSelection.add(i);
                    }
                }

                if (loopSelection.isEmpty()) {
                    continue;
                }

                double[][] W1 =
                        subsetColumns(Winit, loopSelection);

                double[] weights = nnls(W1, genomes);
                double[] newSample = multiply(W1, weights);

                double weightSum = sum(weights);
                if (weightSum == 0.0) {
                    continue;
                }

                double[] normalized = new double[weights.length];
                for (int i = 0; i < weights.length; i++) {
                    normalized[i] =
                            weights[i] / weightSum * sum(genomes);
                }

                List<Double> expanded = new ArrayList<Double>();
                for (double value : normalized) {
                    expanded.add(value);
                }

                // newExposure.insert(i, 0)
                expanded.add(localIndex, 0.0);

                // insert original zero positions in sorted order
                Collections.sort(initialZerosIdx);
                for (Integer zero : initialZerosIdx) {
                    expanded.add(zero, 0.0);
                }

                double[] newExposure = new double[expanded.size()];
                for (int i = 0; i < expanded.size(); i++) {
                    newExposure[i] = expanded.get(i);
                }

                double newSimilarity =
                        relativeL2(genomes, newSample);

                double difference =
                        newSimilarity - originalSimilarity;

                // Preserve upstream behavior exactly.
                if (difference < 0.0) {
                    difference = cutoff + 1e-100;
                }

                if (difference < bestDifference) {
                    bestDifference = difference;
                    bestExposure = newExposure;
                    bestSimilarity = newSimilarity;
                }
            }

            if (bestExposure == null || bestDifference > cutoff) {
                break;
            }

            successExposure = bestExposure;
            successDistance = bestSimilarity;

            if (nonZeroIndicesExact(successExposure).size() == 1) {
                break;
            }

            backgroundSig = getChangedBackgroundSigIdxExact(
                    successExposure,
                    backgroundSig
            );

            originalSimilarity = bestSimilarity;
        }

        RemoveResult result = new RemoveResult();

        if (successExposure == null) {
            result.exposure = oldExposures;
            result.distance = originalSimilarity;
        } else {
            result.exposure = successExposure;
            result.distance = successDistance;
        }

        return result;
    }

    /**
     * Preserves get_changed_background_sig_idx semantics, including value-based
     * lookup through list.index as used by decompose_subroutines.get_indeces.
     */
    private List<Integer> getChangedBackgroundSigIdxExact(
            double[] exposures,
            List<Integer> backgroundSigs) {

        List<Double> backgroundSigValues =
                new ArrayList<Double>();

        for (Integer index : backgroundSigs) {

            if (index >= 0 && index < exposures.length) {
                backgroundSigValues.add(exposures[index]);
            }
        }

        List<Double> tempExposures =
                new ArrayList<Double>();

        for (double exposure : exposures) {

            if (exposure != 0.0) {
                tempExposures.add(exposure);
            }
        }

        /*
         * SigProfiler removes background signatures whose
         * exposure became zero.
         */
        List<Double> activeBackgroundValues =
                new ArrayList<Double>();

        for (Double value : backgroundSigValues) {

            if (value != 0.0) {
                activeBackgroundValues.add(value);
            }
        }

        List<Integer> remappedBackground =
                new ArrayList<Integer>();

        for (Double backgroundValue : activeBackgroundValues) {

            int localIndex =
                    tempExposures.indexOf(backgroundValue);

            if (localIndex >= 0) {
                remappedBackground.add(localIndex);
            }
        }

        return remappedBackground;
    }

    private List<Integer> addConnectedSigsExact(
            List<Integer> backgroundSigs,
            List<String> allSigIds) {

        String[][] connected = new String[][]{
                {"SBS2", "SBS13"},
                {"SBS7a", "SBS7b", "SBS7c", "SBS7d"},
                {"SBS10a", "SBS10b"},
                {"SBS17a", "SBS17b"}
        };

        Set<String> backgroundNames = new HashSet<String>();

        for (Integer index : backgroundSigs) {
            if (index >= 0 && index < allSigIds.size()) {
                backgroundNames.add(allSigIds.get(index));
            }
        }

        List<String> connectNames = new ArrayList<String>();

        for (String[] group : connected) {
            boolean intersects = false;

            for (String name : group) {
                if (backgroundNames.contains(name)) {
                    intersects = true;
                    break;
                }
            }

            if (intersects) {
                connectNames.addAll(Arrays.asList(group));
            }
        }

        backgroundNames.addAll(connectNames);

        List<Integer> result = new ArrayList<Integer>();

        for (String name : backgroundNames) {
            int index = allSigIds.indexOf(name);
            if (index >= 0) {
                result.add(index);
            }
        }

        Collections.sort(result);
        return result;
    }

    private NnlsFit fitSubsetRaw(
            double[][] W,
            double[] genome,
            List<Integer> selected) {

        double[][] subset = subsetColumns(W, selected);
        double[] weights = nnls(subset, genome);
        double[] rawReconstruction = multiply(subset, weights);

        double weightSum = sum(weights);
        double[] normalized = new double[weights.length];

        if (weightSum != 0.0) {
            for (int i = 0; i < weights.length; i++) {
                normalized[i] =
                        weights[i] / weightSum * sum(genome);
            }
        }

        NnlsFit result = new NnlsFit();
        result.rawReconstruction = rawReconstruction;
        result.normalizedExposure = normalized;
        return result;
    }

    private double[][] subsetColumns(
            double[][] matrix,
            List<Integer> columns) {

        double[][] selected =
                new double[matrix.length][columns.size()];

        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < columns.size(); col++) {
                selected[row][col] =
                        matrix[row][columns.get(col)];
            }
        }

        return selected;
    }

    private double[] expandExposure(
            double[] compact,
            List<Integer> selected,
            int fullLength) {

        double[] result = new double[fullLength];

        for (int i = 0; i < selected.size(); i++) {
            result[selected.get(i)] = compact[i];
        }

        return result;
    }

    private double[] numpyRoundAndCorrectMaximum(
            double[] values,
            long maxMutation) {

        double[] rounded = new double[values.length];

        if (values.length == 0) {
            return rounded;
        }

        int maxIndex = 0;

        for (int i = 0; i < values.length; i++) {
            rounded[i] = Math.rint(values[i]);

            if (values[i] > values[maxIndex]) {
                maxIndex = i;
            }
        }

        if (Math.round(sum(rounded)) != maxMutation) {
            rounded[maxIndex] =
                    Math.rint(rounded[maxIndex])
                            + maxMutation
                            - sum(rounded);
        }

        return rounded;
    }

    /**
     * Exact roundConserveSum:
     * total = np.round(np.sum(x))
     * x_out = np.ceil(x)
     * order = np.argsort(x - x_out)
     * x_out[order[:sum(x_out)-total]] -= 1
     */
    private double[] roundConserveSumExact(double[] values) {

        long total = Math.round(Math.rint(sum(values)));

        final double[] residual = new double[values.length];
        double[] output = new double[values.length];
        List<Integer> order = new ArrayList<Integer>();

        long ceilTotal = 0L;

        for (int i = 0; i < values.length; i++) {
            output[i] = Math.ceil(values[i]);
            residual[i] = values[i] - output[i];
            ceilTotal += (long) output[i];
            order.add(i);
        }

        Collections.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                int c = Double.compare(residual[a], residual[b]);
                if (c != 0) return c;
                return Integer.compare(a, b);
            }
        });

        int removeCount =
                (int) (ceilTotal - total + 1e-10);

        for (int i = 0; i < removeCount; i++) {
            output[order.get(i)] -= 1.0;
        }

        return output;
    }

    private List<Integer> nonZeroIndicesExact(double[] values) {
        List<Integer> result = new ArrayList<Integer>();

        for (int i = 0; i < values.length; i++) {
            if (values[i] != 0.0) {
                result.add(i);
            }
        }

        return result;
    }

    private List<Integer> positiveIndicesExact(double[] values) {
        List<Integer> result = new ArrayList<Integer>();

        for (int i = 0; i < values.length; i++) {
            if (values[i] > 0.0) {
                result.add(i);
            }
        }

        return result;
    }

    private List<Integer> zeroIndicesExact(double[] values) {
        List<Integer> result = new ArrayList<Integer>();

        for (int i = 0; i < values.length; i++) {
            if (values[i] == 0.0) {
                result.add(i);
            }
        }

        return result;
    }

    private int[] nonZeroIndexArrayExact(double[] values) {
        List<Integer> indices = nonZeroIndicesExact(values);
        int[] result = new int[indices.size()];

        for (int i = 0; i < indices.size(); i++) {
            result[i] = indices.get(i);
        }

        return result;
    }

    private boolean numpyAll(int[] values) {
        for (int value : values) {
            if (value == 0) {
                return false;
            }
        }

        return true;
    }

    private List<Integer> unionIndices(
            List<Integer> first,
            List<Integer> second) {

        Set<Integer> union = new HashSet<Integer>(first);
        union.addAll(second);

        List<Integer> result = new ArrayList<Integer>(union);
        Collections.sort(result);

        return result;
    }

    private List<Integer> differenceIndices(
            List<Integer> first,
            List<Integer> second) {

        Set<Integer> difference = new HashSet<Integer>(first);
        difference.removeAll(second);

        List<Integer> result = new ArrayList<Integer>(difference);
        Collections.sort(result);

        return result;
    }

    private double relativeL2(
            double[] original,
            double[] reconstructed) {

        return l2Norm(original, reconstructed)
                / l2Norm(original, new double[original.length]);
    }

    // Lawson-Hanson style active-set NNLS, matching scipy.optimize.nnls semantics.
// -------------------------------------------------------------------------
// Lawson-Hanson active-set NNLS.
// Designed to follow scipy.optimize.nnls active-set behaviour more closely.
// -------------------------------------------------------------------------

    private double[] nnls(double[][] a, double[] b) {

        int n = a[0].length;

        double[] x = new double[n];
        boolean[] passive = new boolean[n];

        /*
         * scipy.optimize.nnls default maxiter is approximately 3 * n.
         *
         * Do not use 30*n*n here. Continuing the active-set search far beyond
         * SciPy's iteration behaviour can produce a different active signature
         * set when COSMIC signatures are highly correlated.
         */
        int maxIterations = Math.max(1, 3 * n);
        int iterations = 0;

        while (iterations < maxIterations) {

            double[] residual = subtract(b, multiply(a, x));
            double[] w = transposeMultiply(a, residual);

            int entering = -1;
            double maxW = EPS;

            /*
             * np.argmax-style behaviour:
             * keep the first index when values are equal.
             */
            for (int j = 0; j < n; j++) {

                if (!passive[j] && w[j] > maxW) {
                    maxW = w[j];
                    entering = j;
                }
            }

            if (entering < 0) {
                break;
            }

            passive[entering] = true;

            while (true) {

                int[] passiveSet = passiveIndices(passive);

                if (passiveSet.length == 0) {
                    break;
                }

                double[][] passiveMatrix =
                        selectColumns(a, passiveSet);

                double[] passiveSolution =
                        leastSquares(passiveMatrix, b);

                double[] z = new double[n];

                for (int i = 0; i < passiveSet.length; i++) {
                    z[passiveSet[i]] = passiveSolution[i];
                }

                boolean allPositive = true;

                for (int index : passiveSet) {

                    if (z[index] <= 0.0) {
                        allPositive = false;
                        break;
                    }
                }

                if (allPositive) {
                    x = z;
                    break;
                }

                double alpha = Double.POSITIVE_INFINITY;

                for (int index : passiveSet) {

                    if (z[index] <= 0.0) {

                        double denominator =
                                x[index] - z[index];

                        if (denominator > 0.0) {

                            double candidateAlpha =
                                    x[index] / denominator;

                            if (candidateAlpha < alpha) {
                                alpha = candidateAlpha;
                            }
                        }
                    }
                }

                if (Double.isInfinite(alpha)) {
                    alpha = 0.0;
                }

                for (int j = 0; j < n; j++) {

                    x[j] =
                            x[j]
                                    + alpha
                                    * (z[j] - x[j]);

                    if (Math.abs(x[j]) <= EPS) {
                        x[j] = 0.0;
                    }
                }

                for (int index : passiveSet) {

                    if (x[index] <= EPS) {
                        passive[index] = false;
                        x[index] = 0.0;
                    }
                }

                iterations++;

                if (iterations >= maxIterations) {
                    break;
                }
            }

            iterations++;

            if (iterations >= maxIterations) {
                break;
            }
        }

        for (int i = 0; i < x.length; i++) {

            if (x[i] < 0.0 && x[i] > -1e-9) {
                x[i] = 0.0;
            }
        }

        return x;
    }


    private double[] leastSquares(double[][] a, double[] b) {

        RealMatrix matrix =
                new Array2DRowRealMatrix(a, false);

        RealVector vector =
                new ArrayRealVector(b, false);

        /*
         * scipy.optimize.nnls internally solves active-set least-squares
         * systems using QR-oriented numerical routines.
         *
         * QR is a closer numerical choice here than SVD.
         */
        DecompositionSolver qr =
                new QRDecomposition(
                        matrix,
                        1e-12
                ).getSolver();

        if (qr.isNonSingular()) {
            return qr.solve(vector).toArray();
        }

        /*
         * COSMIC signatures can be nearly collinear.
         * Use SVD only as a rank-deficient fallback.
         */
        DecompositionSolver svd =
                new SingularValueDecomposition(matrix)
                        .getSolver();

        return svd.solve(vector).toArray();
    }

    private int[] passiveIndices(boolean[] passive) {
        int count = 0;

        for (boolean value : passive) {
            if (value) count++;
        }

        int[] indices = new int[count];
        int position = 0;

        for (int i = 0; i < passive.length; i++) {
            if (passive[i]) {
                indices[position++] = i;
            }
        }

        return indices;
    }

    private double[][] selectColumns(
            double[][] matrix,
            int[] columns) {

        double[][] selected =
                new double[matrix.length][columns.length];

        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < columns.length; col++) {
                selected[row][col] =
                        matrix[row][columns[col]];
            }
        }

        return selected;
    }

    // -------------------------------------------------------------------------
    // COSMIC matrix
    // -------------------------------------------------------------------------

    private CosmicMatrix loadCosmicMatrix() throws Exception {

        File file = new File(COSMIC_PATH);

        if (!file.exists()) {
            throw new RuntimeException("COSMIC SBS matrix not found: " + COSMIC_PATH);
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {

            String header = reader.readLine();

            if (header == null) {
                throw new RuntimeException("COSMIC SBS matrix is empty.");
            }

            String[] headerTokens = header.split("\\t", -1);
            List<String> signatureNames = new ArrayList<String>();

            for (int i = 1; i < headerTokens.length; i++) {
                signatureNames.add(headerTokens[i].trim());
            }

            List<String> contexts = new ArrayList<String>();
            List<double[]> rows = new ArrayList<double[]>();

            String line;

            while ((line = reader.readLine()) != null) {

                if (line.trim().length() == 0) {
                    continue;
                }

                String[] tokens = line.split("\\t", -1);

                if (tokens.length != signatureNames.size() + 1) {
                    throw new RuntimeException(
                            "Invalid COSMIC row column count: " + line
                    );
                }

                contexts.add(tokens[0].trim());

                double[] row = new double[signatureNames.size()];

                for (int i = 1; i < tokens.length; i++) {
                    row[i - 1] = Double.parseDouble(tokens[i]);
                }

                rows.add(row);
            }

            if (rows.size() != 96) {
                throw new RuntimeException(
                        "Expected 96 COSMIC SBS contexts, found " + rows.size()
                );
            }

            double[][] matrix = new double[rows.size()][];

            for (int i = 0; i < rows.size(); i++) {
                matrix[i] = rows.get(i);
            }

            return new CosmicMatrix(matrix, signatureNames, contexts);
        }
    }

    // -------------------------------------------------------------------------
    // Response construction
    // -------------------------------------------------------------------------

    private List<SBS96Response.SignatureContribution> buildSignatureContributions(
            double[] activities,
            List<String> signatureNames) {

        List<SBS96Response.SignatureContribution> results =
                new ArrayList<SBS96Response.SignatureContribution>();

        double total = sum(activities);

        for (int i = 0; i < activities.length; i++) {

            if (activities[i] <= 0.0) {
                continue;
            }

            double contributionPercentage =
                    activities[i] / Math.max(total, EPS) * 100.0;

            SBS96Response.SignatureContribution result =
                    new SBS96Response.SignatureContribution();

            result.setSignature(signatureNames.get(i));
            result.setContributionPercentage(contributionPercentage);
            result.setEtiology(getSignatureMeaning(signatureNames.get(i)));
            result.setContributionLevel(
                    contributionLevel(contributionPercentage)
            );

            results.add(result);
        }

        Collections.sort(
                results,
                new Comparator<SBS96Response.SignatureContribution>() {
                    @Override
                    public int compare(
                            SBS96Response.SignatureContribution a,
                            SBS96Response.SignatureContribution b) {
                        return Double.compare(
                                b.getContributionPercentage(),
                                a.getContributionPercentage()
                        );
                    }
                }
        );

        return results;
    }

    private String contributionLevel(double contributionPercentage) {

        if (contributionPercentage >= 20.0) {
            return "Dominant";
        }

        if (contributionPercentage >= 5.0) {
            return "Intermediate";
        }

        return "Minor";
    }

    private String buildClinicalSummary(
            List<SBS96Response.SignatureContribution> signatures) {

        StringBuilder summary = new StringBuilder();

        for (SBS96Response.SignatureContribution result : signatures) {
            summary.append(result.getSignature())
                    .append(" (")
                    .append(String.format(
                            Locale.US,
                            "%.1f",
                            result.getContributionPercentage()
                    ))
                    .append("%) : ")
                    .append(result.getEtiology())
                    .append(". ");
        }

        return summary.toString();
    }

    private Map<String, Map<String, Integer>> groupSbs96(
            Map<String, Integer> sbs96) {

        Map<String, Map<String, Integer>> grouped =
                new LinkedHashMap<String, Map<String, Integer>>();

        for (Map.Entry<String, Integer> entry : sbs96.entrySet()) {

            String key = entry.getKey();
            String mutation = key.substring(
                    key.indexOf('[') + 1,
                    key.indexOf(']')
            );

            char left = key.charAt(0);
            char ref = mutation.charAt(0);
            char right = key.charAt(key.length() - 1);

            String context = "" + left + ref + right;

            if (!grouped.containsKey(mutation)) {
                grouped.put(
                        mutation,
                        new LinkedHashMap<String, Integer>()
                );
            }

            grouped.get(mutation).put(context, entry.getValue());
        }

        return grouped;
    }

    private Map<String, Map<String, Double>> buildSbs96Percentages(
            Map<String, Integer> sbs96) {

        Map<String, Map<String, Double>> grouped =
                new LinkedHashMap<String, Map<String, Double>>();

        double total = 0.0;

        for (Integer value : sbs96.values()) {
            total += value;
        }

        if (total == 0.0) {
            return grouped;
        }

        for (Map.Entry<String, Integer> entry : sbs96.entrySet()) {

            String key = entry.getKey();
            String mutation = key.substring(
                    key.indexOf('[') + 1,
                    key.indexOf(']')
            );

            char left = key.charAt(0);
            char ref = mutation.charAt(0);
            char right = key.charAt(key.length() - 1);

            String context = "" + left + ref + right;

            if (!grouped.containsKey(mutation)) {
                grouped.put(
                        mutation,
                        new LinkedHashMap<String, Double>()
                );
            }

            grouped.get(mutation).put(
                    context,
                    entry.getValue() / total * 100.0
            );
        }

        return grouped;
    }


    private String getSignatureMeaning(String signature) {
        return ETIOLOGY_MAP.getOrDefault(signature, "Unknown");
    }

    // -------------------------------------------------------------------------
    // Metrics
    // -------------------------------------------------------------------------

    private double cosineSimilarity(double[] a, double[] b) {
        double dot = 0.0;
        double aa = 0.0;
        double bb = 0.0;

        for (int i = 0; i < a.length; i++) {
            dot += a[i] * b[i];
            aa += a[i] * a[i];
            bb += b[i] * b[i];
        }

        if (aa == 0.0 || bb == 0.0) {
            return 0.0;
        }

        return dot / (Math.sqrt(aa) * Math.sqrt(bb));
    }

    private double pearson(double[] a, double[] b) {
        if (a.length < 2) {
            return 0.0;
        }

        double value = new PearsonsCorrelation().correlation(a, b);

        return Double.isNaN(value) ? 0.0 : value;
    }

    private double rmse(double[] a, double[] b) {
        double value = 0.0;

        for (int i = 0; i < a.length; i++) {
            double difference = a[i] - b[i];
            value += difference * difference;
        }

        return Math.sqrt(value / a.length);
    }

    private double l1Norm(double[] a, double[] b) {
        double value = 0.0;

        for (int i = 0; i < a.length; i++) {
            value += Math.abs(a[i] - b[i]);
        }

        return value;
    }

    private double l2Norm(double[] a, double[] b) {
        double value = 0.0;

        for (int i = 0; i < a.length; i++) {
            double difference = a[i] - b[i];
            value += difference * difference;
        }

        return Math.sqrt(value);
    }

    private double klDivergence(double[] original, double[] reconstructed) {

        double originalTotal = sum(original);
        double reconstructedTotal = sum(reconstructed);

        if (originalTotal <= 0.0 || reconstructedTotal <= 0.0) {
            return 0.0;
        }

        double kl = 0.0;

        for (int i = 0; i < original.length; i++) {

            if (original[i] <= 0.0) {
                continue;
            }

            double p = original[i] / originalTotal;
            double q = reconstructed[i] / reconstructedTotal;

            if (q <= 0.0) {
                return Double.POSITIVE_INFINITY;
            }

            kl += p * Math.log(p / q);
        }

        return kl;
    }

    private double[] normalize(double[] vector) {
        double total = sum(vector);
        double[] normalized = new double[vector.length];

        if (total == 0.0) {
            return normalized;
        }

        for (int i = 0; i < vector.length; i++) {
            normalized[i] = vector[i] / total;
        }

        return normalized;
    }

    // -------------------------------------------------------------------------
    // Matrix/vector helpers
    // -------------------------------------------------------------------------

    private double[] multiply(double[][] matrix, double[] vector) {
        double[] result = new double[matrix.length];

        for (int row = 0; row < matrix.length; row++) {
            double value = 0.0;

            for (int col = 0; col < vector.length; col++) {
                value += matrix[row][col] * vector[col];
            }

            result[row] = value;
        }

        return result;
    }

    private double[] transposeMultiply(double[][] matrix, double[] vector) {
        double[] result = new double[matrix[0].length];

        for (int col = 0; col < matrix[0].length; col++) {
            double value = 0.0;

            for (int row = 0; row < matrix.length; row++) {
                value += matrix[row][col] * vector[row];
            }

            result[col] = value;
        }

        return result;
    }

    private double[] subtract(double[] a, double[] b) {
        double[] result = new double[a.length];

        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    private double sum(double[] values) {
        double total = 0.0;

        for (double value : values) {
            total += value;
        }

        return total;
    }

    // -------------------------------------------------------------------------
    // Inner classes
    // -------------------------------------------------------------------------

    private static class Mutation {
        private final String chromosome;
        private final long position;
        private final char ref;
        private final char alt;

        private Mutation(
                String chromosome,
                long position,
                char ref,
                char alt) {
            this.chromosome = chromosome;
            this.position = position;
            this.ref = ref;
            this.alt = alt;
        }
    }

    private static class CosmicMatrix {
        private final double[][] matrix;
        private final List<String> signatureNames;
        private final List<String> contexts;

        private CosmicMatrix(
                double[][] matrix,
                List<String> signatureNames,
                List<String> contexts) {
            this.matrix = matrix;
            this.signatureNames = signatureNames;
            this.contexts = contexts;
        }

        private double[][] getMatrix() {
            return matrix;
        }

        private List<String> getSignatureNames() {
            return signatureNames;
        }

        private List<String> getContexts() {
            return contexts;
        }
    }

    private static class NnlsFit {
        private double[] rawReconstruction;
        private double[] normalizedExposure;
    }

    private static class AddResult {
        private double[] exposure;
        private double distance;
    }

    private static class RemoveResult {
        private double[] exposure;
        private double distance;
    }

    private static class AddRemoveResult {
        private List<Integer> background;
        private double[] exposure;
        private double distance;
    }

    private static class AssignmentResult {
        private double[] activities;
        private double[] reconstructed;
        private double cosine;
        private double pearson;
        private double rmse;
        private double l1Norm;
        private double l1NormPercent;
        private double l2Norm;
        private double l2NormPercent;
        private double klDivergence;
        private List<String> selectedSignatureNames;

        private double[] getActivities() {
            return activities;
        }

        private double getCosine() {
            return cosine;
        }

        private double getPearson() {
            return pearson;
        }

        private double getRmse() {
            return rmse;
        }

        private double getL1Norm() {
            return l1Norm;
        }

        private double getL1NormPercent() {
            return l1NormPercent;
        }

        private double getL2Norm() {
            return l2Norm;
        }

        private double getL2NormPercent() {
            return l2NormPercent;
        }

        private double getKlDivergence() {
            return klDivergence;
        }

        private List<String> getSelectedSignatureNames() {
            return selectedSignatureNames;
        }
    }
}
