# PROJECT_DESCRIPTION

This document is optimized for fast code navigation and safe edits in this repository.

## 1) What this project is

- Name: `EncyclopeDIA` (Java/Maven monorepo-style app).
- Domain: DIA proteomics search/quantification plus related tools.
- Core mode: library search of DIA data (`Encyclopedia`).
- Additional modes:
  - `Walnut`/`Pecanpie`: FASTA-driven DIA search (PECAN-style).
  - `XCorDIA`: DIA search with XCorr-like scoring (variant mode by default).
  - `Thesaurus`: phosphopeptide positional isomer localization.
  - `Scribe`: DDA library search.
  - `CLIConverter`: format conversion utilities.
  - `Batch`: run queued jobs from `.encxml`.

Primary launcher routing lives in:
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/Main.java`

## 2) Fast mental model

Common execution pattern across tools:

1. Parse CLI args (`CommandLineParser`, `InstrumentSpecificSearchParameters`).
2. Build typed `SearchParameters` (`SearchParameterParser` / `PecanParameterParser` / tool-specific).
3. Build `JobData` object (contains IO paths, percolator files, parameters, factories).
4. Read DIA/MZML into internal stripe representation (`StripeFile*`, `DIAFileReader`).
5. Score candidates (`algorithms/*` scoring tasks/factories).
6. Run/consume Percolator outputs.
7. Write reports + optional library artifacts (`.elib`, `.dlib`, `.blib` conversions).

## 3) Entrypoint map (what to run / where to debug)

| Mode | Class | Trigger flag(s) from `Main` | Main output suffix |
|---|---|---|---|
| Global launcher | `Main` | N/A | delegates |
| EncyclopeDIA | `Encyclopedia` | default or `-encyclopedia` | `.encyclopedia.txt` |
| Walnut/PECAN | `Pecanpie` / `Walnut` | `-walnut`, `-pecan` | `.pecan.txt` |
| XCorDIA | `XCorDIA` | `-xcordia` | `.xcordia.txt` |
| Thesaurus | `Thesaurus` | `-thesaurus` | `.thesaurus.txt` |
| Scribe | `Scribe` | `-scribe` | `.scribe.txt` |
| Converter hub | `CLIConverter` | `-convert` | converter-specific |
| Batch XML | `Batch` | `-batch` | job-dependent |

Reference classes:
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/Encyclopedia.java`
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/Pecanpie.java`
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/XCorDIA.java`
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/Thesaurus.java`
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/Scribe.java`
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/CLIConverter.java`
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/Batch.java`

## 4) Package navigation index

Use this as the primary "where should I look?" map.

- `src/main/java/edu/washington/gs/maccoss/encyclopedia/algorithms`
  - scoring, alignment, quantitation, percolator integration, model/prediction code.
  - subpackages to know quickly:
    - `library`, `pecan`, `xcordia`, `phospho`, `scribe`, `percolator`, `alignment`, `quantitation`.
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/filereaders`
  - DIA/MZML parsing, library IO, stripe file processing, converters from external formats.
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/filewriters`
  - report writers, library exporters, downstream file generation.
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/datastructures`
  - core domain models and parameter containers.
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/gui`
  - Swing UI framework, parameter panels, browser tooling.
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/cli`
  - individual converter command implementations.
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/jobs`
  - queued/batch execution model (`WorkerJob`, XML driver support).
- `src/main/java/edu/washington/gs/maccoss/encyclopedia/utils`
  - logging, threading, mass spec utilities, parsing helpers.

## 5) File format + artifact cheatsheet

Key extensions and their roles:

- Raw/processed inputs: `.mzML`, `.DIA`
- Libraries:
  - `.dlib` (`LibraryFile.DLIB`)
  - `.elib` (`LibraryFile.ELIB`)
  - `.blib` (`BlibFile.BLIB`)
- Batch driver: `.encxml`
- Common output reports:
  - EncyclopeDIA: `.encyclopedia.txt`
  - EncyclopeDIA v2: `.encyclopedia2.txt`
  - Walnut/PECAN: `.pecan.txt`
  - XCorDIA: `.xcordia.txt`
  - Thesaurus: `.thesaurus.txt`
  - Scribe: `.scribe.txt`
- Common sidecars: `.log`, `.features.txt`, decoy/protein reports (mode-specific).

## 6) Converter capability map

`CLIConverter` routes to:

- Prosit/Spectronaut/DIA-NN CSV -> library (`-prositCSVToLibrary`)
- BLIB -> library (`-blibToLib`)
- MSP/SPTXT -> library (`-mspToLib`)
- OpenSwath TSV -> library (`-openswathTSVToLibrary`)
- Library -> BLIB (`-libraryToBlib`)
- Merge libraries (`-mergeLibraries`)
- FASTA -> Prosit CSV (`-fastaToPrositCSV`)
- FASTA -> Koina/Prosit library (`-fastaToKoinaLibrary`)
- DIA preprocessing/merge (`-processDIA`)
- Library PTM adjustment (`-adjustLibraryForPTMs`)

## 7) Build/test navigation

Build system:
- Maven (`pom.xml`)
- Java source/target: 17.

Notable Maven behavior:
- Integration tests are skipped by default via `-DskipITs=true`.
- Enable ITs with `-DskipITs=false`.
- Shaded executable jars are configured via build profiles (`buildGlobal`, `buildEncyclopedia`, and additional tool profiles in `pom.xml`).

Useful commands:

```bash
# compile quickly
mvn -DskipITs=true -DskipTests compile

# unit tests
mvn test

# include integration tests
mvn verify -DskipITs=false

# build executable jar profile (example)
mvn package -PbuildEncyclopedia
```

Test footprint:
- `src/test/java` contains ~190+ test files (unit + integration + end-to-end patterns).
- End-to-end entry anchor: `AbstractEndToEndIT`.

## 8) Grep-first navigation shortcuts

Use these patterns first before deep reading:

```bash
# list top-level executables
rg -n "public static void main" src/main/java/edu/washington/gs/maccoss/encyclopedia

# find mode routing decisions
rg -n "containsKey\\(\"-.*\"\\)" src/main/java/edu/washington/gs/maccoss/encyclopedia/Main.java

# find output filename constants
rg -n "OUTPUT_FILE_SUFFIX|DECOY_FILE_SUFFIX|FEATURE_FILE_SUFFIX" src/main/java/edu/washington/gs/maccoss/encyclopedia/algorithms

# locate converter handlers
rg -n "class Convert|main\\(String\\[] args\\)" src/main/java/edu/washington/gs/maccoss/encyclopedia/cli

# locate percolator touchpoints
rg -n "Percolator" src/main/java/edu/washington/gs/maccoss/encyclopedia

# locate GUI parameter panels
rg -n "ParametersPanel|SearchPanel" src/main/java/edu/washington/gs/maccoss/encyclopedia/gui
```

## 9) Change-impact guide

If changing X, verify Y:

- CLI flags / argument parsing:
  - files: `Main.java`, mode class (`Encyclopedia.java`, etc.), parser classes.
  - verify: help text, required-arg checks, GUI compatibility.
- scoring logic:
  - files: `algorithms/*Scoring*`, `*JobData`, `Percolator*`.
  - verify: feature outputs, threshold handling, decoy/target behavior.
- file formats / library IO:
  - files: `filereaders/LibraryFile.java`, `filereaders/BlibFile.java`, converter classes.
  - verify: round-trip reading/writing and version metadata compatibility.
- batch / XML serialization:
  - files: `jobs/XMLDriverFactory.java`, `*JobData.writeToXML/readFromXML`.
  - verify: old XML compatibility + null/path guards.
- GUI panels:
  - files: `gui/framework/*ParametersPanel*`, `SearchPanel*`.
  - verify: default values match CLI defaults.

## 10) Pragmatic pitfalls to remember

- Many classes are static-heavy and orchestrator-heavy. Prefer minimal, targeted edits and verify side effects.
- Similar logic is duplicated across tool entrypoints (`Encyclopedia`, `Scribe`, `XCorDIA`, `Pecanpie`, `Thesaurus`); check sibling modes for parity when editing shared behavior.
- `bin/` mirrors many repository files; treat `src/main/java` as source of truth for code changes.
- `target/` artifacts may exist and be stale; do not trust generated outputs over source constants.

## 11) Personal quick-start sequence for future sessions

When re-entering this repo, do:

1. `rg -n "public static void main" src/main/java/edu/washington/gs/maccoss/encyclopedia`
2. `rg -n "class .*JobData|OUTPUT_FILE_SUFFIX" src/main/java/edu/washington/gs/maccoss/encyclopedia/algorithms`
3. Open `Main.java` then mode-specific entrypoint for current task.
4. Trace into `filereaders` and `algorithms` package matching that mode.
5. Run smallest relevant tests first, then broader suite if behavior changed.

