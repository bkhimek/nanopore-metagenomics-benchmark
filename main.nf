#!/usr/bin/env nextflow
// Nanopore metagenomics benchmark: core path (Kraken2 -> score). See README, section "Nextflow".
nextflow.enable.dsl = 2
include { CORE } from './workflows/core'
workflow { CORE() }
