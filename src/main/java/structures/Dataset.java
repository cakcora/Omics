package structures;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Dataset {
    ArrayList<Row> rows = new ArrayList();
    private static final Logger logger = LoggerFactory.getLogger(Dataset.class);
    public void addRow(Row row) {
        rows.add(row);
    }

    public void removeAllNans() {
        List<Integer> remove =new ArrayList<>();
        int index=0;
        for(Row r:rows){
            double[] a= Arrays.stream(r.getData()).distinct().toArray();
            if(a.length==1&&Double.isNaN(a[0])){
               remove.add(index);
            }
            index++;
        }
        //logger.info(remove.size()+"/"+rows.size()+" PTms will be removed");
        Collections.sort(remove, Collections.reverseOrder());
        for(int r=0;r<remove.size();r++){

            int o = remove.get(r);
            rows.remove(o);
        }
        logger.info(rows.size()+" PTms exist in the dataset");
    }

    public ArrayList<Row> getRows() {
        return rows;
    }
}
