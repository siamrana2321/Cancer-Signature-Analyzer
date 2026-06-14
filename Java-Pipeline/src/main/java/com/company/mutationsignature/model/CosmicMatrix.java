package com.company.mutationsignature.model;

import java.util.List;

public class CosmicMatrix {

    private double[][] matrix;

    private List<String> signatureNames;

    private List<String> contexts;

    public CosmicMatrix(
            double[][] matrix,
            List<String> signatureNames,
            List<String> contexts
    ) {
        this.matrix = matrix;
        this.signatureNames = signatureNames;
        this.contexts = contexts;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public List<String> getSignatureNames() {
        return signatureNames;
    }

    public List<String> getContexts() {
        return contexts;
    }
}