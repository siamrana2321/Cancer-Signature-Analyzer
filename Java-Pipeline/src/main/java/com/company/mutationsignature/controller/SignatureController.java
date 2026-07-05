package com.company.mutationsignature.controller;

import com.company.mutationsignature.model.SignatureResponse;
import com.company.mutationsignature.service.SignatureService;

import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Collections;

@RestController
@RequestMapping("/api/signature")
public class SignatureController {
    private static final Logger log = LoggerFactory.getLogger(SignatureController.class);

    private final SignatureService service;
    private static final String SAMPLE_DIRECTORY =
            "D:\\Mohammad_Siam_Ahmed_Rana\\Mutation_Signature\\Resources\\sample\\";

    public SignatureController(
            SignatureService service
    ) {
        this.service = service;
    }

//    @PostMapping("/analyze")
//    public SignatureResponse analyze(
//            @RequestParam("file")
//            MultipartFile file
//    ) throws Exception {
//
//        return service.analyze(file);
//
//    }

    @PostMapping("/analyze")
    public ResponseEntity<?> analyze(
            @RequestParam("vcfFileName") String fileName
    ) {
        final String filePath = SAMPLE_DIRECTORY + fileName;
        try {

            SignatureResponse response =
                    service.analyze(filePath);

            return ResponseEntity.ok(response);

        }
        catch (Exception e) {

            return ResponseEntity
                    .status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body(
                            Collections.singletonMap(
                                    "error",
                                    e.getMessage()
                            )
                    );
        }
    }
}