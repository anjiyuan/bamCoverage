/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package qut.bamcoverage.include;

/**
 *
 * @author an
 */
public class LocNucl {
        public String chr;
    public int loc;
    public String nucl;
    public int count;

    public LocNucl(String chr, int loc, String nucl) {
        this.chr = chr;
        this.loc = loc;
        this.nucl = nucl;
        this.count = 1;
    }

    public void increment() {
        this.count++;
    }

    @Override
    public String toString() {
        return chr + "\t" + loc + "\t" + nucl + "\t" + count;
    }
}
