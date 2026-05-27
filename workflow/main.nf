#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TEST {
    output:
    stdout

    script:
    """
    echo "Bryophyte PCD pipeline initialized successfully"
    """
}

workflow {
    TEST()
}