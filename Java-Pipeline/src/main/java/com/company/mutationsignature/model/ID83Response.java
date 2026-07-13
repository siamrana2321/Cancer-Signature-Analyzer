package com.company.mutationsignature.model;

import java.util.List;
import java.util.Map;

public class ID83Response {

    private String cosmicVersion;
    private long totalMutations;

    private Map<String, Integer> id83Spectrum;
    private Map<String, Double> id83Percentage;
    private Map<String, Map<String, Double>> id83Grouped;

    private List<SignatureContribution> signatureContributions;

    private String clinicalSummary;

    private double reconstructionCosine;
    private double pearson;
    private double l1Norm;
    private double l1NormPercent;
    private double l2Norm;
    private double l2NormPercent;
    private double klDivergence;

    // -------------------------------------------------------------------------
    // Signature contribution
    // -------------------------------------------------------------------------

    public static class SignatureContribution {

        public String signature;
        public double contributionPercentage;
        public String etiology;
        public String contributionLevel;

        public SignatureContribution(
                String signature,
                double contributionPercentage,
                String etiology,
                String contributionLevel) {

            this.signature = signature;
            this.contributionPercentage = contributionPercentage;
            this.etiology = etiology;
            this.contributionLevel = contributionLevel;
        }
    }

    // -------------------------------------------------------------------------
    // Getters / Setters
    // -------------------------------------------------------------------------

    public String getCosmicVersion() {
        return cosmicVersion;
    }

    public void setCosmicVersion(String cosmicVersion) {
        this.cosmicVersion = cosmicVersion;
    }

    public long getTotalMutations() {
        return totalMutations;
    }

    public void setTotalMutations(long totalMutations) {
        this.totalMutations = totalMutations;
    }

    public Map<String, Integer> getId83Spectrum() {
        return id83Spectrum;
    }

    public void setId83Spectrum(Map<String, Integer> id83Spectrum) {
        this.id83Spectrum = id83Spectrum;
    }

    public Map<String, Double> getId83Percentage() {
        return id83Percentage;
    }

    public void setId83Percentage(Map<String, Double> id83Percentage) {
        this.id83Percentage = id83Percentage;
    }

    public Map<String, Map<String, Double>> getId83Grouped() {
        return id83Grouped;
    }

    public void setId83Grouped(
            Map<String, Map<String, Double>> id83Grouped) {
        this.id83Grouped = id83Grouped;
    }

    public List<SignatureContribution> getSignatureContributions() {
        return signatureContributions;
    }

    public void setSignatureContributions(
            List<SignatureContribution> signatureContributions) {
        this.signatureContributions = signatureContributions;
    }

    public String getClinicalSummary() {
        return clinicalSummary;
    }

    public void setClinicalSummary(String clinicalSummary) {
        this.clinicalSummary = clinicalSummary;
    }

    public double getReconstructionCosine() {
        return reconstructionCosine;
    }

    public void setReconstructionCosine(double reconstructionCosine) {
        this.reconstructionCosine = reconstructionCosine;
    }

    public double getPearson() {
        return pearson;
    }

    public void setPearson(double pearson) {
        this.pearson = pearson;
    }

    public double getL1Norm() {
        return l1Norm;
    }

    public void setL1Norm(double l1Norm) {
        this.l1Norm = l1Norm;
    }

    public double getL1NormPercent() {
        return l1NormPercent;
    }

    public void setL1NormPercent(double l1NormPercent) {
        this.l1NormPercent = l1NormPercent;
    }

    public double getL2Norm() {
        return l2Norm;
    }

    public void setL2Norm(double l2Norm) {
        this.l2Norm = l2Norm;
    }

    public double getL2NormPercent() {
        return l2NormPercent;
    }

    public void setL2NormPercent(double l2NormPercent) {
        this.l2NormPercent = l2NormPercent;
    }

    public double getKlDivergence() {
        return klDivergence;
    }

    public void setKlDivergence(double klDivergence) {
        this.klDivergence = klDivergence;
    }
}