# 🧬 bamCoverage

`bamCoverage` is a Java-based tool for analyzing BAM files with a reference genome. It provides high-performance coverage calculation, supporting large datasets and customizable memory usage.

---

## 🚀 Quick Start

### a. Prerequisites

1. Install **Java 21**
2. Install **Maven 3.9.9**

### b. Download and Compile the Code

```bash
git clone https://github.com/anjiyuan/bamCoverage.git
cd bamCoverage
mvn clean package
```
### c. Run bamCoverage
```bash
java -Xms20G -Djava.library.path=target/lib -cp target/bamCoverage-1.0.jar qut.bamcoverage.bamcoverage xxx.genome.fasta xxx.bam
java -Xms20G -Djava.library.path=target/lib -cp target/bamCoverage-1.0.jar qut.bamcoverage.codingRegion_SNP xxx.gff3 xxx.bam
```
Replace xxx.genome.fasta xxx.gff3 and xxx.bam with your actual genome FASTA, annotation file and BAM file paths.
