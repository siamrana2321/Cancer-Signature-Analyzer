package com.company.mutationsignature.controller;

import com.company.mutationsignature.model.SBS96Response;
import com.company.mutationsignature.service.SBS96SignatureService;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.Collections;
import com.company.mutationsignature.model.ID83Response;
import com.company.mutationsignature.service.ID83SignatureService;
import org.springframework.beans.factory.annotation.Autowired;
import java.nio.file.Paths;

@RestController
@RequestMapping("/api/MutationSignature")
public class MutationSignatureController {
    private static final Logger log = LoggerFactory.getLogger(MutationSignatureController.class);

    @Autowired
    private final SBS96SignatureService SBS96service;
    @Autowired
    private ID83SignatureService ID83service;

    private static final String SAMPLE_DIRECTORY =
            "D:\\horizondb\\resources\\refDB\\MutationSignatureDB\\sample\\";

    public MutationSignatureController(SBS96SignatureService SBS96service, ID83SignatureService ID83service) {
        this.SBS96service = SBS96service;
        this.ID83service = ID83service;
    }

    @GetMapping("/SBS96analyze")
    public ResponseEntity<?> analyze(@RequestParam("vcfFileName") String fileName) {
        final String filePath = SAMPLE_DIRECTORY + fileName;
        try {
            SBS96Response response = SBS96service.analyze(filePath);
            return ResponseEntity.ok(response);
        } catch (Exception e) {
            log.error("Analysis failed", e);
            return ResponseEntity
                    .status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body(Collections.singletonMap("error", e.getMessage()));
        }
    }

    @GetMapping("/ID83analyze")
    public ResponseEntity<ID83Response> ID83RunDefault(@RequestParam("vcfFileName") String fileName) {
        final String filePath = SAMPLE_DIRECTORY + fileName;
        try {
            ID83Response response = ID83service.runPipeline(
                    Paths.get(filePath)
            );
            return ResponseEntity.ok(response);
        } catch (Exception e) {
            return ResponseEntity.internalServerError().build();
        }
    }

}
