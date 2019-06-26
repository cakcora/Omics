package structures;

import java.util.Objects;

public class Row {

    private final Gene gene;
    private final PTMType ptmType;
    private final Position position;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Row row = (Row) o;
        return Objects.equals(gene, row.gene) &&
                Objects.equals(ptmType, row.ptmType) &&
                Objects.equals(position, row.position);
    }

    @Override
    public int hashCode() {

        return Objects.hash(gene, ptmType, position);
    }

    @Override
    public String toString() {
        return "Row{" +
                "gene=" + gene +
                ", ptmType=" + ptmType +
                ", position=" + position +
                '}';
    }

    private final double[] data;

    public Row(Gene gene, PTMType mod, Position pos, String[] strings) {
    this.gene = gene;
    this.ptmType= mod;
    this.position= pos;
    this.data = new double[strings.length];
        for (int i = 0; i < strings.length; i++) {
            if(strings[i].equalsIgnoreCase("NA")){
                data[i] = Double.NaN;
            }
            else data[i]=Double.parseDouble(strings[i]);
        }

    }

    public double[] getData() {
        return data;
    }


    public String getPTMName() {
        return gene+""+ptmType+""+position;
    }
}
