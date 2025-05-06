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
java -Xms20G -Djava.library.path=target/lib -cp target/bamCoverage-1.0.jar qut.bamcoverage.bamcoverage_GPT xxx.genome.fasta xxx.bam
```
Replace xxx.genome.fasta and xxx.bam with your actual genome FASTA and BAM file paths.
