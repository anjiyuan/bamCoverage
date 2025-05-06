/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package qut.bamcoverage;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import qut.libgff3.full_locus;
import qut.libgff3.gene_str;
import qut.libgff3.gff3;

/**
 *
 * @author an
 */
public class codingRegion_SNP {
//   String path = "C:\\Jiyuan\\sourceCode\\NB_annotation\\NbLx02\\SNP_zuba";
    
    public static void main(String[] args) throws IOException {
        String gff3_fn = "C:\\Jiyuan\\sourceCode\\NB_annotation\\NbLx02\\SNP_zuba\\NbLx02.gff3";
        String gene_prefix = "LX";
        String SNP_prefix = "C:\\Jiyuan\\sourceCode\\NB_annotation\\NbLx02\\SNP_zuba\\NbLx02.Zuba.2pass.Aligned.sortedByCoord.out.bam";
//        System.out.println(SNP_prefix);
        gff3 mygff3 = new gff3();
        mygff3.readFullIn(gff3_fn);
        new codingRegion_SNP().proc(mygff3, SNP_prefix, gene_prefix);
    }
    
    void proc(gff3 mygff3, String SNP_prefix, String gene_prefix) throws IOException{
//        for (int chr_no = 1; chr_no <=19; chr_no++){
//            String chrom = gene_prefix+String.format("%02d", chr_no);
//            String mismatch_fn = SNP_prefix+String.format("%02d", chr_no)+"_mismatches.txt";
//            System.out.println(mismatch_fn);
//            one_chrom(mygff3, mismatch_fn);
//        }
//            String mismatch_fn = SNP_prefix+"."+chrom+"_mismatches.txt";
//            System.out.println(mismatch_fn);
            one_chrom(mygff3, SNP_prefix);
        
    }
    
    void one_chrom(gff3 mygff3, String SNP_prefix) throws IOException{
        for(String chrom : mygff3.genes.keySet()){
            String mismatch_fn = SNP_prefix+"."+chrom+"_mismatches.txt";
//        for (int chr_no = 1; chr_no <=19; chr_no++){
//            String chr_str = String.format("%02d", chr_no);
            List<gene_str> genes = mygff3.genes.get(chrom);
            List<rang> CDS = new ArrayList<>();
            for(gene_str gs : genes){
                for (int i = 0; i < gs.mRNA.size(); i++){
                    for(full_locus fl : gs.CDSs(i)){
                        CDS.add(new rang(fl.start, fl.stop));
                    }
                }
            }
            Collections.sort(CDS, new Comparator<rang>(){
                @Override
                public int compare(rang o1, rang o2) {
                    int ret = o1.stop - o2.stop;
                    if(ret == 0){
                        ret = o1.start - o2.start;
                    }
                    return ret;
                }
            });
//            String fn = path+"/"+SNP_prefix+chr_str+"_mismatches.txt";
//            BufferedWriter bw = new BufferedWriter(new FileWriter(fn+".CDS.txt"));
            BufferedReader br = new BufferedReader(new FileReader(mismatch_fn));
            br.readLine();
            String line;
            int homo = 0;
            int heter = 0;
            while ((line = br.readLine()) != null) {
                mismatch_str ms = new mismatch_str(line);
                int po = Collections.binarySearch( CDS, new rang(ms.pos, ms.pos),new Comparator<rang>() {
                    @Override
                    public int compare(rang o1, rang o2) {
                        int ret = o1.stop - o2.stop;
                        if (ret == 0) {
                            ret = o1.start - o2.start;
                        }
                        return ret;
                    }
                });
                if (po < 0){
                    po = -po - 1;
                }
                if ((po < CDS.size())&& (ms.pos > CDS.get(po).start)) {//located in CDS
                    if (ms.max_ACGT() >= ms.total() * 0.8) {
                        homo++;
                    } else {
                        heter++;
                    }
                }
            }
            br.close();
            System.out.println(chrom + ":homo=\t" + homo + "\theter=\t"+heter);
        }
    }
}

class rang{
    public int start;
    public int stop;

    public rang(int start, int stop) {
        this.start = start;
        this.stop = stop;
    }
    
}
class mismatch_str{
    public String chr;
    public int pos;
    public int ref_num;
    public int[] ACGT_num = new int[4];
    public mismatch_str(String line) {
        String[] strarray = line.split("\t");
        chr = strarray[0];
        pos = Integer.parseInt(strarray[1]);
        ref_num = Integer.parseInt(strarray[2].split(":")[1]);
        for(int i = 0; i < ACGT_num.length; i++){
            ACGT_num[i] = Integer.parseInt(strarray[3+i].split(":")[1]);
        }
    }

    public int total() {
        int ret = 0;
        for (int i = 0; i < ACGT_num.length; i++) {
            ret += ACGT_num[i];
        }
        return ret;
    }
    public int max_ACGT(){
        int ret = ACGT_num[0];
        for(int i = 1; i < ACGT_num.length; i++){
            if(ret < ACGT_num[i]){
                ret = ACGT_num[i];
            }
        }
        return ret;
    }
}
