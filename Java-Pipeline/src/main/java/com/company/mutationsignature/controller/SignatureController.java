package com.company.mutationsignature.controller;

import com.company.mutationsignature.model.SignatureResponse;
import com.company.mutationsignature.service.SignatureService;

import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.util.Collections;

@RestController
@RequestMapping("/api/signature")
public class SignatureController {

    private final SignatureService service;

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
            @RequestParam("file") MultipartFile file
    ) {

        try {

            SignatureResponse response = service.analyze(file);

            return ResponseEntity.ok(response);

        } catch (Exception e) {

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