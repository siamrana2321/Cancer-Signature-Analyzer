package com.company.mutationsignature.model;

import java.util.List;
import java.util.Map;

public class SBS96Response {

    private String cosmicVersion;
    private Integer totalMutations;
    private Map<String, Integer> sbs96Spectrum;
    private Map<String, Map<String, Integer>> sbs96Grouped;
    private Map<String, Map<String, Double>> sbs96Percentage;
    private List<SignatureResult> cosmicContributions;
    private List<BiologicalAnnotation> biologicalAnnotations;
    private String clinicalSummary;
    private Double reconstructionCosine;
    private Double rmse;
    private Double pearson;

    // Getters and Setters
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

    public List<SignatureResult> getCosmicContributions() { return cosmicContributions; }
    public void setCosmicContributions(List<SignatureResult> cosmicContributions) { this.cosmicContributions = cosmicContributions; }

    public List<BiologicalAnnotation> getBiologicalAnnotations() { return biologicalAnnotations; }
    public void setBiologicalAnnotations(List<BiologicalAnnotation> biologicalAnnotations) { this.biologicalAnnotations = biologicalAnnotations; }

    public String getClinicalSummary() { return clinicalSummary; }
    public void setClinicalSummary(String clinicalSummary) { this.clinicalSummary = clinicalSummary; }

    public Double getReconstructionCosine() { return reconstructionCosine; }
    public void setReconstructionCosine(Double reconstructionCosine) { this.reconstructionCosine = reconstructionCosine; }

    public Double getRmse() { return rmse; }
    public void setRmse(Double rmse) { this.rmse = rmse; }

    public Double getPearson() { return pearson; }
    public void setPearson(Double pearson) { this.pearson = pearson; }

    // ---------- Nested DTOs ----------
    public static class SignatureResult {
        private String signature;
        private double similarity;
        private String etiology;

        public String getSignature() { return signature; }
        public void setSignature(String signature) { this.signature = signature; }

        public double getSimilarity() { return similarity; }
        public void setSimilarity(double similarity) { this.similarity = similarity; }

        public String getEtiology() { return etiology; }
        public void setEtiology(String etiology) { this.etiology = etiology; }
    }

    public static class BiologicalAnnotation {
        private String signature;
        private String etiology;
        private String confidence;

        public String getSignature() { return signature; }
        public void setSignature(String signature) { this.signature = signature; }

        public String getEtiology() { return etiology; }
        public void setEtiology(String etiology) { this.etiology = etiology; }

        public String getConfidence() { return confidence; }
        public void setConfidence(String confidence) { this.confidence = confidence; }
    }
}