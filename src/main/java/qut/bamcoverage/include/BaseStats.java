/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package qut.bamcoverage.include;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author an
 */
public class BaseStats {
        public String chr;
    public int loc;
    public char refBase;
    public Map<Character, Integer> counts;

    public BaseStats(String chr, int loc, char refBase) {
        this.chr = chr;
        this.loc = loc;
        this.refBase = refBase;
        counts = new HashMap<>();
        for (char b : new char[]{'A', 'C', 'G', 'T'}) {
            counts.put(b, 0);
        }
    }
    public int totalReads(){
        return counts.get('A') + counts.get('C') + counts.get('G') + counts.get('T');
    }
    public int refCount(){
        return counts.getOrDefault(refBase, 0);
    }

    public void increment(char base) {
        counts.computeIfPresent(base, (k, v) -> v + 1);
    }

    public String format() {
        return chr + "\t" + loc + "\t" +
                refBase + ":" + counts.getOrDefault(refBase, 0) + "\t" +
                "A:" + counts.get('A') + "\t" +
                "C:" + counts.get('C') + "\t" +
                "G:" + counts.get('G') + "\t" +
                "T:" + counts.get('T');
    }
}
