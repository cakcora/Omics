import structures.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class PTMLoader {
    String sep = "\t";
    String regex = ";";
    String reg2 = " ";
    public PTMLoader() {


    }

    public Dataset parse(String input, int startCol, int length) throws IOException {
        String line ="";
        BufferedReader dataFile = new BufferedReader(new FileReader(input));
        Dataset set = new Dataset();
        String header = dataFile.readLine();
        while((line=dataFile.readLine())!=null){
            String arr[] = line.split(this.sep);
            String ptmName = arr[0];

            if(ptmName.contains(regex)){
                //multiple ptms on the same row
                for(String pname: ptmName.split(regex)){
                    Row row = addPtm(pname.trim(),arr,startCol,length);
                    set.addRow(row);
                }
            }
           else {
                Row row = addPtm(ptmName.trim(),arr,startCol,length);
                set.addRow(row);
            }

        }
        dataFile.close();
        return set;
    }

    private Row addPtm(String pname, String[] arr, int startCol, int length) {
        String a []= pname.split(reg2);
        Gene gene = new Gene(a[0]);
        PTMType mod = new PTMType(a[1]);
        Position pos = new Position(a[2]);
        Row row = new Row(gene,mod,pos,Arrays.copyOfRange(arr,startCol,startCol+length));
       return row;
    }


}
