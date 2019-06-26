package v1;

import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class DegreeDeviation {
    private static final Logger logger = LoggerFactory.getLogger(DegreeDeviation.class);
    public static void main(String[] args) throws IOException {
        int start[] = new int[]{1,6,11,16,21,26,31,36,41, 46,51,56,61,66,72};
        int len[] = new int[]{5,5,5,5,5,5,5,5,5,5,5,5,5,6,6};
        String dataFile = args[0];
        String line = "";
        BufferedReader br = new BufferedReader(new FileReader(dataFile));
        String header[] = br.readLine().replaceAll("\"", "").split("\t");
        int count = 0;
        int valNA = -1;
        while (br.readLine() != null) {
            count++;
        }
        logger.info(header.length+" columns");
        double[][] data = new double[count][header.length];
        int[][] proteinDegrees = new int[count][start.length];
        Map<Integer,String> proteinNames = new HashMap<Integer,String>();
        br = new BufferedReader(new FileReader(dataFile));
        br.readLine();
        int index = 0;
        int d=0;
        while ((line = br.readLine()) != null) {
            String arr[] = line.replaceAll("\"", "").split("\t");
            String protein = arr[0];
            proteinNames.put(index,protein);
            for (int i = 1; i < header.length; i++) {
                String v = arr[i];
                if (v.equalsIgnoreCase("NA")) {
                    data[index][i-1] = valNA;
                }
                else {
                    double value = Double.parseDouble(v);
                    data[index][i-1] = (value);
                }
            }
            index++;
        }
        logger.info(index+ "proteins");
        double corThreshold = 0.7d;
        {
            for (int layer = 0; layer <15; layer++) {
                int edgeCount = 0;
                UndirectedSparseGraph<Integer,Integer> graph = new UndirectedSparseGraph<Integer,Integer>();
                int sIndex = start[layer];
                int dlen = len[layer];
                for (int j = 0; j < count; j++) {
                    double[] a1 = Arrays.copyOfRange(data[j], sIndex , sIndex + dlen);
                    if(allNA(a1,valNA)) {
                        continue;
                    }
                    if (!graph.containsVertex(j)) graph.addVertex(j);
                    for (int k = 1 + j; k < count; k++) {
                        double[] a2 = Arrays.copyOfRange(data[k], sIndex , sIndex + dlen);
                        if(allNA(a2,valNA))continue;
                        if (!graph.containsVertex(k)) graph.addVertex(k);
                        if(a1.length!=dlen)logger.error(layer+"th layer has "+a1.length+", but we need "+dlen);
                        double corr = new PearsonsCorrelation().correlation(a1, a2);
                        if (!Double.isNaN(corr) && corr > corThreshold) {
                            graph.addEdge(edgeCount++, j, k);
///                       logger.info(layer+"\t"+j+"\t"+k+"\t"+corr);
                        }
                    }
                }
                for (int v =0;v<count;v++) {
                    int neighborCount=-1;
                    if(graph.containsVertex(v))
                        neighborCount= graph.getNeighborCount(v);
                    proteinDegrees[v][layer] = neighborCount;
                }

            }
        }
        for (int proteinID = 0; proteinID < proteinDegrees.length; proteinID++) {
            String proteinName = proteinNames.get(proteinID);
            double meanDegree = 0d;
            int exists=0;
            for (int layer = 0; layer < proteinDegrees[0].length; layer++) {
                if(proteinDegrees[proteinID][layer]!=valNA){
                    exists++;
                    meanDegree+= proteinDegrees[proteinID][layer];
                }
            }
            meanDegree = meanDegree/exists;
            double deviation=0d;
            for (int layer = 0; layer < proteinDegrees[0].length; layer++) {
                if (proteinDegrees[proteinID][layer] != valNA) {
                    deviation += Math.pow(proteinDegrees[proteinID][layer] - meanDegree, 2);
                }
            }
            deviation = Math.pow(deviation/exists, 0.5);
            if(!Double.isNaN(deviation)) {
                String str = proteinName + "\t" + deviation + "\t" + meanDegree+"\t"+exists;
                System.out.println(str);
            }
        }



    }

    private static boolean allNA(double[] a1, int valNA) {
        if((int)a1[0]==valNA&&(int)a1[1]==valNA){
            return true;
        }
        return false;
    }
}
