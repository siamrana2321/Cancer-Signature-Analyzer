package com.company.mutationsignature.controller;

import com.company.mutationsignature.model.ID83Response;
import com.company.mutationsignature.service.MutationSignatureService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.nio.file.Paths;

@RestController
@RequestMapping("/api/pipeline")
public class PipelineController {

    @Autowired
    private MutationSignatureService service;

    @Value("${mutation.signature.vcf.path}")
    private String defaultVcfPath;

    @Value("${mutation.signature.genome.path}")
    private String defaultGenomePath;

    @Value("${mutation.signature.cosmic.path}")
    private String defaultCosmicPath;

    @GetMapping("/ID83run")
    public ResponseEntity<ID83Response> ID83RunDefault() {
        try {
            ID83Response response = service.runPipeline(
                    Paths.get(defaultVcfPath),
                    Paths.get(defaultGenomePath),
                    Paths.get(defaultCosmicPath)
            );
            return ResponseEntity.ok(response);
        } catch (Exception e) {
            return ResponseEntity.internalServerError().build();
        }
    }
}