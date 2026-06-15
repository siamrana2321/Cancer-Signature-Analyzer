package com.company.mutationsignature.model;

import java.util.List;
import java.util.Map;

public class SignatureResponse {

    private String cosmicVersion;
    private Integer totalMutations;
    private Map<String, Integer> sbs96Spectrum;
    private Map<String, Map<String, Integer>> sbs96Grouped;
    private Map<String,Map<String,Double>> sbs96Percentage;
    private List<SignatureResult> cosmicContributions;
    private List<BiologicalAnnotation> biologicalAnnotations;
    //    private Integer signaturesDetected;
    private String clinicalSummary;
    private Double reconstructionCosine;
    private Double rmse;
    private Double pearson;


    public String getCosmicVersion() {
        return cosmicVersion;
    }
    public void setCosmicVersion(String cosmicVersion) {
        this.cosmicVersion = cosmicVersion;
    }

    public Integer getTotalMutations() {
        return totalMutations;
    }

    public void setTotalMutations(Integer totalMutations) {
        this.totalMutations = totalMutations;
    }

    public Map<String, Integer> getSbs96Spectrum() {
        return sbs96Spectrum;
    }

    public void setSbs96Spectrum(Map<String, Integer> sbs96Spectrum) {
        this.sbs96Spectrum = sbs96Spectrum;
    }

    public Map<String, Map<String, Double>> getSbs96Percentage() {
        return sbs96Percentage;
    }

    public void setSbs96Percentage(Map<String, Map<String, Double>> sbs96Percentage) {
        this.sbs96Percentage = sbs96Percentage;
    }

    public Map<String, Map<String, Integer>> getSbs96Grouped() {
        return sbs96Grouped;
    }

    public void setSbs96Grouped(Map<String, Map<String, Integer>> sbs96Grouped) {
        this.sbs96Grouped = sbs96Grouped;
    }

    public List<SignatureResult> getCosmicContributions() {
        return cosmicContributions;
    }
    public void setCosmicContributions(List<SignatureResult> cosmicContributions) {
        this.cosmicContributions = cosmicContributions;
    }

    public List<BiologicalAnnotation> getBiologicalAnnotations() {
        return biologicalAnnotations;
    }

    public void setBiologicalAnnotations(List<BiologicalAnnotation> biologicalAnnotations) {
        this.biologicalAnnotations = biologicalAnnotations;
    }

    public String getClinicalSummary() {
        return clinicalSummary;
    }

    public void setClinicalSummary(String clinicalSummary) {
        this.clinicalSummary = clinicalSummary;
    }

//    public Integer getSignaturesDetected() {
//        return signaturesDetected;
//    }
//
//    public void setSignaturesDetected(Integer signaturesDetected) {
//        this.signaturesDetected = signaturesDetected;
//    }

    public Double getReconstructionCosine() {
        return reconstructionCosine;
    }

    public void setReconstructionCosine(Double reconstructionCosine) {
        this.reconstructionCosine = reconstructionCosine;
    }

    public Double getRmse() {
        return rmse;
    }

    public void setRmse(Double rmse) {
        this.rmse = rmse;
    }

    public Double getPearson() {
        return pearson;
    }

    public void setPearson(Double pearson) {
        this.pearson = pearson;
    }

}