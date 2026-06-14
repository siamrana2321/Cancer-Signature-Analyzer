package com.company.mutationsignature.model;

public class NnlsResult {

    private double[] weights;

    private double[] reconstructed;

    private double cosine;

    private double rmse;

    private double pearson;

    public double[] getWeights() {
        return weights;
    }

    public void setWeights(double[] weights) {
        this.weights = weights;
    }

    public double[] getReconstructed() {
        return reconstructed;
    }

    public void setReconstructed(double[] reconstructed) {
        this.reconstructed = reconstructed;
    }

    public double getCosine() {
        return cosine;
    }

    public void setCosine(double cosine) {
        this.cosine = cosine;
    }

    public double getRmse() {
        return rmse;
    }

    public void setRmse(double rmse) {
        this.rmse = rmse;
    }

    public double getPearson() {
        return pearson;
    }

    public void setPearson(double pearson) {
        this.pearson = pearson;
    }
}