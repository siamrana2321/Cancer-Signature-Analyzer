package com.company.mutationsignature.model;

import java.util.List;
import java.util.Map;

public class SBS96Response {

    private String cosmicVersion;
    private Integer totalMutations;
    private Map<String, Integer> sbs96Spectrum;
    private Map<String, Map<String, Integer>> sbs96Grouped;
    private Map<String, Map<String, Double>> sbs96Percentage;
    private List<SignatureContribution> signatureContributions;
    private String clinicalSummary;

    private Double reconstructionCosine;
    private Double rmse;
    private Double pearson;

    // SigProfilerAssignment reconstruction metrics
    private Double l1Norm;
    private Double l1NormPercent;
    private Double l2Norm;
    private Double l2NormPercent;
    private Double klDivergence;

    public String getCosmicVersion() { return cosmicVersion; }
    public void setCosmicVersion(String cosmicVersion) { this.cosmicVersion = cosmicVersion; }

    public Integer getTotalMutations() { return totalMutations; }
    public void setTotalMutations(Integer totalMutations) { this.totalMutations = totalMutations; }

    public Map<String, Integer> getSbs96Spectrum() { return sbs96Spectrum; }
    public void setSbs96Spectrum(Map<String, Integer> sbs96Spectrum) { this.sbs96Spectrum = sbs96Spectrum; }

    public Map<String, Map<String, Integer>> getSbs96Grouped() { return sbs96Grouped; }
    public void setSbs96Grouped(Map<String, Map<String, Integer>> sbs96Grouped) { this.sbs96Grouped = sbs96Grouped; }

    public Map<String, Map<String, Double>> getSbs96Percentage() { return sbs96Percentage; }
    public void setSbs96Percentage(Map<String, Map<String, Double>> sbs96Percentage) { this.sbs96Percentage = sbs96Percentage; }

    public List<SignatureContribution> getSignatureContributions() { return signatureContributions; }
    public void setSignatureContributions(List<SignatureContribution> signatureContributions) {
        this.signatureContributions = signatureContributions;
    }

    public String getClinicalSummary() { return clinicalSummary; }
    public void setClinicalSummary(String clinicalSummary) { this.clinicalSummary = clinicalSummary; }

    public Double getReconstructionCosine() { return reconstructionCosine; }
    public void setReconstructionCosine(Double reconstructionCosine) { this.reconstructionCosine = reconstructionCosine; }

    public Double getRmse() { return rmse; }
    public void setRmse(Double rmse) { this.rmse = rmse; }

    public Double getPearson() { return pearson; }
    public void setPearson(Double pearson) { this.pearson = pearson; }

    public Double getL1Norm() { return l1Norm; }
    public void setL1Norm(Double l1Norm) { this.l1Norm = l1Norm; }

    public Double getL1NormPercent() { return l1NormPercent; }
    public void setL1NormPercent(Double l1NormPercent) { this.l1NormPercent = l1NormPercent; }

    public Double getL2Norm() { return l2Norm; }
    public void setL2Norm(Double l2Norm) { this.l2Norm = l2Norm; }

    public Double getL2NormPercent() { return l2NormPercent; }
    public void setL2NormPercent(Double l2NormPercent) { this.l2NormPercent = l2NormPercent; }

    public Double getKlDivergence() { return klDivergence; }
    public void setKlDivergence(Double klDivergence) { this.klDivergence = klDivergence; }

    public static class SignatureContribution {
        private String signature;
        private double contributionPercentage;
        private String etiology;
        private String contributionLevel;

        public String getSignature() { return signature; }
        public void setSignature(String signature) { this.signature = signature; }

        public double getContributionPercentage() { return contributionPercentage; }
        public void setContributionPercentage(double contributionPercentage) {
            this.contributionPercentage = contributionPercentage;
        }

        public String getEtiology() { return etiology; }
        public void setEtiology(String etiology) { this.etiology = etiology; }

        public String getContributionLevel() { return contributionLevel; }
        public void setContributionLevel(String contributionLevel) {
            this.contributionLevel = contributionLevel;
        }
    }
}
