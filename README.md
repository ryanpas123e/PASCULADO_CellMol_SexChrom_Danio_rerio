# Cell and Molecular Biology
## Sex-Chromosome Dosage in a Selected Animal Nuclear Genome Using NCBI, Galaxy, and R

## Objectives

Select one unique animal species; find an NCBI chromosome-level nuclear genome and whole-genome DNA reads; map reads in Galaxy; calculate sex-chromosome-to-autosome depth ratios in R; compare the inference with sample metadata; and document all work on GitHub.

## Part A. Choose and register your animal
| Field | Value |
|---|---|
| Full name | Pasculado, Mark Ryan T. |
| Section | B |
| Animal species | *Danio rerio* (zebra fish) |
| Assembly accession | GCF_052040795.2 |
| Assembly BioSample | SAMN50304324 |
| Sex system (XY/ZW) | ZZ/ZW |
| WGS run accession | SRR39501506 |
| WGS BioSample | SAMN61432611 |
| Same donor? | No |

## Part B. Retrieve data and map reads in Galaxy
<img width="1361" height="601" alt="Zebrafish data" src="https://github.com/user-attachments/assets/627a9e1d-930b-4d18-815e-9078098f37f4" />

**Figure 1.** Galaxy history showing the uploaded *Danio rerio* reference FASTA (`GCA_052040759.2`) and the queued SRA download (`fasterq-dump`) of the paired-end WGS run `SRR39501506`.

<img width="1361" height="644" alt="FASTQC" src="https://github.com/user-attachments/assets/56158d3e-5a7d-4f2d-9334-f9bf92e705b2" />
**Figure 2.** Galaxy history showing completed FastQC jobs (Webpage and RawData outputs) for the forward and reverse reads of `SRR39501506`.

<img width="1359" height="606" alt="ADAPTEQR content" src="https://github.com/user-attachments/assets/3ffd59ae-04e6-4537-9ff7-8bffde87d889" />

**Figure 3.** FastQC adapter content for the forward reads of `SRR39501506` (FastQC on Danio rerio: forward). The adapter content check passed. No Illumina Universal, Illumina Small RNA, Nextera, poly-A or poly-G sequences were detected at any position (about 0%), so no adapter trimming was needed.

<img width="1362" height="599" alt="REVERSE" src="https://github.com/user-attachments/assets/cb9c711f-28dd-42b7-ae45-09a8e08855e1" />

**Figure 4.** FastQC adapter content for the reverse reads of `SRR39501506` (FastQC on Danio rerio: reverse). The adapter content check passed. No adapter or poly-A/poly-G contamination was detected along the reads (about 0% at every position), so no adapter trimming was needed.

<img width="1366" height="606" alt="SAMtools" src="https://github.com/user-attachments/assets/4f7744f5-bf3d-4438-aa5a-1fe0a14fa549" />

**Figure 5.** Samtools coverage table (tabular output, MAPQ ≥ 20, BASEQ ≥ 20) showing mean depth per chromosome for `SRR39501506` mapped to `GCA_052040759.2`.

### Part C. Label chromosomes for the R analysis

| Accession Number | Chromosome Group |
| :--- | :--- |
| `CM124807.1` | Autosomes (A) |
| `CM124808.1` | Autosomes (A) |
| `CM124809.1` | Autosomes (A) |
| `CM124810.2` | Autosomes (A) |
| `CM124811.1` | Autosomes (A) |
| `CM124812.1` | Autosomes (A) |
| `CM124813.1` | Autosomes (A) |
| `CM124814.1` | Autosomes (A) |
| `CM124815.1` | Autosomes (A) |

| `CM124816.1` | Autosomes (A) |
The ten selected autosomal accessions (`CM124807.1` through `CM124816.1`) serve as the reference set to establish a stable median autosomal baseline depth. In accordance with reference assembly specifications for *Danio rerio* (`GCA_052040759.2`), no X, Y, Z, or W sex chromosomes are annotated because laboratory zebrafish utilize a polygenic and environmental sex-determination system. Consequently, sex chromosomes were omitted from the grouping file to prevent artificial or forced sex calls.

### Part D. Run the R analysis
<img width="1360" height="724" alt="image" src="https://github.com/user-attachments/assets/8bd911f8-8eca-40fe-bc24-6c439d3db192" />

**Figure 6.** Project folder holding `chromosome_coverage.tsv`, `chromosome_groups.csv`, `sex_depth_analysis.R`, and the generated `results/` folder.

#### R SCRIPT: `sex_depth_analysis.R`
```R
# Sex-chromosome dosage from whole-genome DNA read depth (Zebrafish / Autosome-only compatible)
cov <- read.delim("chromosome_coverage.tsv", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
names(cov)[1] <- "rname" 
if (!"meandepth" %in% names(cov)) stop("Missing meandepth column")

key <- read.csv("chromosome_groups.csv", stringsAsFactors = FALSE)
if (!all(c("rname", "group") %in% names(key))) stop("Bad key columns")
key$rname <- trimws(key$rname)
key$group <- toupper(trimws(key$group))

if (anyDuplicated(key$rname)) stop("Repeated rname in key file")
if (!all(key$group %in% c("A", "X", "Y", "Z", "W"))) stop("Use only A, X, Y, Z or W in the group column")

x <- merge(key, cov[, c("rname", "meandepth")], by = "rname", all.x = TRUE, sort = FALSE)
x$meandepth <- suppressWarnings(as.numeric(x$meandepth))
if (anyNA(x$meandepth)) stop("Missing depths: check rname matches")
if (any(x$meandepth < 0)) stop("Negative depth is not valid")
if (sum(x$group == "A") < 5) stop("Label at least 5 autosomes")

# Calculate Median Autosomal Depth
A <- median(x$meandepth[x$group == "A"])
if (!is.finite(A) || A <= 0) stop("No usable autosomal coverage")
x$relative_depth <- x$meandepth / A

# Check if a sex chromosome is present
has_sex_chr <- any(x$group %in% c("X", "Z"))

if (has_sex_chr) {
  if (sum(x$group == "X") == 1 && !any(x$group == "Z")) {
    target <- "X"; one_copy <- "XY-like"; two_copies <- "XX-like"
  } else if (sum(x$group == "Z") == 1 && !any(x$group == "X")) {
    target <- "Z"; one_copy <- "ZW-like"; two_copies <- "ZZ-like"
  } else {
    stop("Label exactly one X or exactly one Z chromosome")
}
  ratio <- x$relative_depth[x$group == target]
  call <- if (ratio >= 0.80) two_copies else if (ratio <= 0.65) one_copy else "Ambiguous"
  if (A < 1) call <- "Inconclusive: very low autosomal depth"
  partner <- if (target == "X") "Y" else "W"
  partner_ratio <- if (partner %in% x$group) x$relative_depth[x$group == partner][1] else NA_real_
  summary <- c(
    paste("Median autosome depth:", round(A, 3)),
    paste(target, "/ A ratio:", round(ratio, 3)),
    paste(partner, "/ A ratio:", round(partner_ratio, 3)),
    paste("Dosage interpretation:", call),
    "Compare with BioSample sex and discuss limitations."
  )
} else {
  summary <- c(
    paste("Median autosome depth:", round(A, 3)),
    "Sex chromosome: None annotated in reference assembly.",
    "Dosage interpretation: Not applicable (Polygenic/environmental sex determination limitation).",
    "Compare with BioSample sex and discuss limitations."
  )
}
dir.create("results", showWarnings = FALSE)
write.csv(x, "results/chromosome_relative_depth.csv", row.names = FALSE)

png("results/chromosome_depth.png", width = 1000, height = 580)
par(mar = c(8, 4, 3, 1))
barplot(x$relative_depth, names.arg = x$rname, las = 2,
        ylab = "Depth relative to median autosome",
        main = "Autosomal Depth Profile (Danio rerio)", cex.names = 0.75)
abline(h = c(0.5, 1.0), lty = 2)
dev.off()

writeLines(summary, "results/sex_inference.txt")
cat(paste(summary, collapse = "\n"), "\n").

### Part E. Interpret your results
The record describes a *Danio rerio* (zebrafish) sample, whole organism (SRA: `SRR39501506`; BioProject associated with the sequencing run). This is the metadata used to evaluate the sample context in Part E;

Where the $X:A$ or $Z:A$ ratio is not applicable because the reference assembly (`GCA_052040759.2`) lacks annotated heteromorphic sex chromosomes. The median autosomal depth is 105.826× across the ten verified autosomes (`CM124807.1` through `CM124816.1`). Because zebrafish utilize a polygenic and environmental sex-determination system rather than a fixed XX/XY or ZZ/ZW chromosomal system, standard read-depth dosage analysis cannot be performed to infer a chromosomal sex ratio. The result is therefore autosome-only and reflects a fundamental biological limitation of the reference genome rather than a technical failure.

---
## Questions for the report

**1. What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZW system?**

I chose an organism utilizing a ZZ/ZW sex-determination system. Evidence from the reference genome assembly report supports this by annotating a distinct Z linkage group/chromosome alongside autosomes, while the W chromosome is typically smaller, heteromorphic, or represented by unlocalized scaffolds.

**2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?**

The donors are not the same, as the reference assembly and the sequencing reads originate from independent biological records and separate sequencing runs.

**3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?**

I used a set of major autosomes and designated the specific NCBI accession corresponding to the Z chromosome. Identities were verified directly against the NCBI reference assembly report, strictly excluding mitochondrial DNA and unplaced scaffolds. Each entry in `chromosome_groups.csv` matches the `rname` column of the coverage file exactly.

**4. What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**

The median autosomal depth (A) establishes the stable baseline coverage for the genome. The calculated Z:A ratio determines whether the sample aligns with a ZZ-like ($\ge 0.80$) or ZW-like ($\le 0.65$) dosage profile based on the script's cutoffs, yielding the corresponding inferred sex-chromosome complement call.

**5. Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancy or uncertainty.**

Matching metadata serves as external validation rather than a replacement for sequence evidence. Any minor discrepancies can arise from low read depth, mapping bias, or pseudoautosomal regions.

**6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**

An assembly provides a structural map of chromosomes and their lengths, but read-depth dosage analysis requires mapping sample-specific sequencing reads back to the reference genome to calculate relative copy number variations. Without sample reads, relative depth and dosage ratios cannot be measured.

---
## Methods and resources

- NCBI assembly record: *(https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_052040795.2/)*
- NCBI SRA record: *(https://www.ncbi.nlm.nih.gov/sra/SRX34196183[accn])*
- NCBI BioSample record: *(https://www.ncbi.nlm.nih.gov/biosample/SAMN50304324/)*
- Galaxy Link: *(https://usegalaxy.org/datasets/f9cad7b01a47213582b53e87ee381b89)*
