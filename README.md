# YB Sensor Analysis Pipeline

Scripts for analyzing piRNA biogenesis sensor libraries in the context of **A Naïve RNA Sampling Core Enables Adaptive piRNA Specificity Against Transposable Elements**.

## Overview

This repository contains the pipeline for processing, quantifying and summarizing sensor library readouts. It is designed to work downstream of the internal **AnnotationPipeline** for small RNA and/or sensor datasets, but also supports a standalone mode using input FASTA libraries.

Typical use cases include:

- Quantifying sensor read counts across libraries  
- Separating sense and antisense signals  
- Comparing sensor activity across experimental conditions  
- Generating summary tables for downstream plotting and statistics  

## Repository Structure

```text
├── script-files/        # Core analysis scripts
├── utility-files/       # Helper scripts, config and utility resources
├── submit.sh            # Main submission script to run the analysis
└── README.md            # This file
```

## Pipeline Components

### 1. Input Modes

The pipeline supports two mutually exclusive input modes:

1. **AnnotationPipeline mode (`-i`)**
   - `-i` is a hyperlink (or path-derived string) pointing to an existing AnnotationPipeline run.
   - Intended for internal runs processed via the established AP framework.

2. **FASTA directory mode (`-I`)**
   - `-I` is a directory containing sensor library FASTA files.
   - Library names are inferred from the FASTA filenames (file stem without extension).
   - Useful for running the sensor analysis on pre-defined libraries without a full AP run.

Only one of `-i` or `-I` may be used per analysis run.

### 2. Quantification and Strand Assignment

Core scripts in `script-files/` handle:

- Mapping or counting reads against sensor libraries  
- Separating counts by **sense** and **antisense** orientation  
- Aggregating per-library statistics  

Exact tools and parameters are configured within the scripts and can be adapted if needed.

### 3. Output and Summary

The pipeline produces:

- Per-library count tables (sense and antisense)  
- Combined summary tables for all libraries  
- Intermediate files in a temporary analysis directory  
- Final outputs in a dedicated `sensor-analysis/<analysis-name>/` folder  

These outputs are intended to be directly usable for downstream plotting and statistics.

## Requirements

- **Apptainer** or **Singularity** (as used on the Brennecke lab infrastructure)

The link to the containers can be found on https://github.com/BrenneckeLab/Handler_2026-YB

If running outside the original cluster setup, you will likely need to adapt hard-coded paths in the scripts (e.g. storage locations, AP base paths, Singularity image directories).

## Usage

The main entry point is:

```bash
bash submit.sh -h

  OPTIONS:
      -i, --input         hyperlink to the AnnotationPipeline run to analyze
                          please copy the full link from the top page of the analysis
                          (exclusive with -I)

      -I, --input-dir     directory containing fasta files to analyze
                          library names will be extracted from filenames
                          (exclusive with -i)

      -N --name           name of the analysis
      -a, --all           set flag for local processing of all libraries
                          by default only the first one will be analyzed if option -C is set
      -f, --force         delete all old data to start a fresh analysis
                          also sample to plasmid reference file will be deleted
      -h, --help          Show this message
      -C, --computing     set flag for local processing (use only if multiple cores available)
      -D, --debug         set debug mode - does not trigger git commit
```

## Related Resources

### Main Publication Repository
https://github.com/BrenneckeLab/Handler_2026-YB


## Citation

Please find the proper citation in https://github.com/BrenneckeLab/Handler_2026-YB


## Contact

For questions or additional information, please contact:
dominik.handler@imba.oeaw.ac.at

## License

MIT License
