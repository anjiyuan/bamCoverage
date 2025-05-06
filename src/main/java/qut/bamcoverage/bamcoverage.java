/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package qut.bamcoverage;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

import java.io.*;
import java.util.*;
import org.qut.post_blast.include.genome_str;
import qut.bamcoverage.include.BaseStats;
import qut.bamcoverage.include.LocNucl;
import java.nio.charset.StandardCharsets;
/**
 *
 * @author an
 */
public class bamcoverage {
    

    public static void main(String[] args) throws IOException {
        System.out.println("args.length="+args.length);
        String ref_fn = "C:\\Jiyuan\\sourceCode\\NB_annotation\\NbLx02\\NbLx02.genome.fasta";
        String bam_fn = "C:\\Jiyuan\\sourceCode\\chamilka\\NbLx02_LX01_Lall.bam";
        int minReads = 20;
        String selectChr = "";
        if (args.length >= 2) {
            ref_fn = args[0];
            bam_fn = args[1];
            if (args.length >= 3) {
                ref_fn = args[0];
                bam_fn = args[1];
                minReads = Integer.parseInt(args[2]);
                if(args.length == 4){
                    selectChr = args[3];
                }
            }
        }        
        new bamcoverage().proc(bam_fn, ref_fn, minReads, selectChr);
    }
    void proc(String bam_fn, String ref_fn, int minReads, String selectChr) throws IOException {    
//        bam_fn = "C:\\Jiyuan\\sourceCode\\chamilka\\test\\debug.bam";
//        bam_fn = "C:\\Jiyuan\\sourceCode\\chamilka\\test\\NbLx02.test.10k.1passAligned.sortedByCoord.out.bam";
        
        
        File bamFile = new File(bam_fn);
//        File refFasta = new File(ref_fn);

        SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
//        ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFasta);
        SAMFileHeader header = reader.getFileHeader();
        genome_str ref_genome = new genome_str();
        for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
            String chr = seq.getSequenceName();
            if(!selectChr.isEmpty() && !chr.equals(selectChr)){
                continue;
            }
            int chrLength = seq.getSequenceLength();
            System.out.println("Processing: " + chr);

            int[] coverage = new int[chrLength];
            Map<String, LocNucl> insertions = new HashMap<>();
            Map<String, BaseStats> mismatches = new HashMap<>();
            Map<String, LocNucl> deletions = new HashMap<>();

            Map<String, StringBuilder> chr_seq = ref_genome.proc(ref_fn, chr);
            byte[] refBases = chr_seq.get(chr).toString().getBytes(StandardCharsets.UTF_8);
//            byte[] refBases = chr_seq.get(chr)refFile.getSequence(chr).getBases();

            SAMRecordIterator it = reader.query(chr, 0, 0, false);
            while (it.hasNext()) {
                SAMRecord record = it.next();
                if (record.getReadUnmappedFlag()) continue;

                byte[] readBases = record.getReadBases();
                int refPos = record.getAlignmentStart();
                int readPos = 0;
//if(record.getReadName().contains("A00808:169:HNTKVDSXX:4:1473:29116:10019")){
//    int debug = 0;
//}
                for (CigarElement ce : record.getCigar().getCigarElements()) {
                    CigarOperator op = ce.getOperator();
                    int len = ce.getLength();

                    switch (op) {
                        case M:
                        case EQ:
                        case X:
                            for (int i = 0; i < len; i++) {
                                if (refPos - 1 < coverage.length) {
                                    coverage[refPos - 1]++;
                                }
                                if ((refPos - 1) < refBases.length && readPos < readBases.length) {
                                    char refBase = (char) refBases[refPos - 1];
                                    char readBase = (char)readBases[readPos];
//if(66305470 == refPos - 1){
//    int debug = 0;
//    System.out.println("ref="+refBase+" read="+readBase+" id="+record.getReadName());
//}                              
//                                    if (record.getReadNegativeStrandFlag()) {
//                                        readBase = (byte) complement((char) readBase);
//                                    }
                                    
//                                    char readBase = (char) readBases[readPos];

                                    final String mismatchKey = chr + "_" + (refPos - 1);
                                    final String mismatchChr = chr;
                                    final int mismatchLoc = refPos - 1;
                                    final char mismatchRefBase = refBase;

                                    mismatches.computeIfAbsent(mismatchKey, k -> new BaseStats(mismatchChr, mismatchLoc, mismatchRefBase))
                                            .increment(readBase);
//                                    mismatches.computeIfAbsent(mismatchKey, k -> new BaseStats(mismatchChr, mismatchLoc, mismatchRefBase))
//                                            .increment(readBase);mismatches.get(2369878)
                                }
                            refPos++;
                            readPos++;
                        }
                            break;

                        case I: {
                            StringBuilder insertSeq = new StringBuilder();
                            for (int i = 0; i < len && readPos < readBases.length; i++) {
                                insertSeq.append((char) readBases[readPos++]);
                            }
//if(refPos - 1 == 120279477){
//    int debug = 0;
//}                           
//                            boolean isReverse = record.getReadNegativeStrandFlag();
//                            String insSeq = isReverse
//                                    ? reverseComplement(insertSeq.toString())
//                                    : insertSeq.toString();
                             String insSeq =insertSeq.toString();
                            final int insertLoc = refPos - 1;
                            final String insertKey = chr + "_" + insertLoc + "_" + insSeq;
                            final String insertChr = chr;

                            insertions.compute(insertKey, (k, v) -> {
                                if (v == null) {
                                    return new LocNucl(insertChr, insertLoc, insSeq);
                                }
                                v.increment();
                                return v;
                            });
//                            String ins = record.getReadNegativeStrandFlag()
//                                    ? reverseComplement(insertSeq.toString())
//                                    : insertSeq.toString();
//                            String key = chr + "_" + (refPos - 1) + "_" + ins;
//                            insertions.compute(key, (k, v) -> {
//                                if (v == null) return new LocNucl(chr, refPos - 1, ins);
//                                v.increment(); return v;
//                            });
                            break;
                        }

                        case D: {
                            for (int i = 0; i < len; i++) {
                                if ((refPos - 1) < refBases.length) {
                                    final int delLoc = refPos - 1;
                                    final char delBase = (char) refBases[delLoc];
                                    final String delKey = chr + "_" + delLoc + "_" + delBase;
                                    final String delChr = chr;

                                    deletions.compute(delKey, (k, v) -> {
                                        if (v == null) {
                                            return new LocNucl(delChr, delLoc, String.valueOf(delBase));
                                        }
                                        v.increment();
                                        return v;
                                    });
                                }
                                refPos++;
                            }
                            break;
                        }

                        case S: readPos += len; break;
                        case N: case H: case P: refPos += len; break;
                    }
                }
            }
            it.close();

            writeCoverageToFile(bam_fn+"."+chr + "_coverage.txt", chr, coverage);
            List<LocNucl> sortedInsertions = new ArrayList<>(insertions.values());
            sortedInsertions.sort(Comparator
                    .comparing((LocNucl l) -> l.chr)
                    .thenComparingInt(l -> l.loc));
            writeLocNuclToFile(bam_fn+"."+chr + "_insertions.txt", sortedInsertions, coverage, "chr\tloc\t#reads\tinserted_nucl\tcount", minReads);
            
            List<LocNucl> sortedDeletions = new ArrayList<>(deletions.values());
            sortedDeletions.sort(Comparator
                    .comparing((LocNucl l) -> l.chr)
                    .thenComparingInt(l -> l.loc));
            writeLocNuclToFile(bam_fn+"."+chr + "_deletions.txt", sortedDeletions, coverage, "chr\tloc\t#reads\tdeleted_nucl\tcount", minReads);
            writeMismatchToFile(bam_fn+"."+chr + "_mismatches.txt", mismatches.values(), minReads);
//            break;
        }

        reader.close();
    }

    public void writeCoverageToFile(String filename, String chr, int[] coverage) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            writer.write("chr\tloc\tcoverage\n");
            for (int i = 0; i < coverage.length; i++) {
                writer.write(chr + "\t" + i + "\t" + coverage[i] + "\n");
            }
        }
    }

    public void writeLocNuclToFile(String filename, List<LocNucl> list, int[] coverage, String header, int minReads) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            writer.write(header + "\n");
//            Collections.sort(list, new Comparator<LocNucl>(){
//                @Override
//                public int compare(LocNucl o1, LocNucl o2) {
//                    int ret = o1.chr.compareTo(o2.chr);
//                    if(ret == 0){
//                        ret = o1.loc - o2.loc;
//                    }
//                    return ret;
//                }
//            });
            for (LocNucl ln : list) {
                if ((coverage[ln.loc] > minReads) && (ln.count > coverage[ln.loc] * 0.2)) {
                    writer.write(ln.chr + "\t"
                            + ln.loc + "\t"
                            + coverage[ln.loc] + "\t"
                            + ln.nucl + "\t"
                            + ln.count + "\n");
                }
//                writer.write(ln.toString() + "\n");
            }
        }
    }

    public void writeMismatchToFile(String filename, Collection<BaseStats> stats, int minReads) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            writer.write("chr\tloc\tref:count\tA:count\tC:count\tG:count\tT:count\n");
            List<BaseStats> sorted_stats = new ArrayList<>(stats);
            Collections.sort(sorted_stats, (BaseStats o1, BaseStats o2) -> {
                int ret = o1.chr.compareTo(o2.chr);
                if(ret == 0){
                    ret = o1.loc - o2.loc;
                }
                return ret;
            });
            for (BaseStats b : sorted_stats) {
                if ((b.totalReads() > minReads) && (b.refCount() < b.totalReads() * 0.8)) {
                    writer.write(b.format() + "\n");
                }
            }
        }
    }

    public String reverseComplement(String seq) {
        StringBuilder sb = new StringBuilder();
        for (int i = seq.length() - 1; i >= 0; i--) {
            char c = seq.charAt(i);
            sb.append(complement(c));
        }
        return sb.toString();
    }

    private static char complement(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            default:
                return 'N';
        }
    }
}
