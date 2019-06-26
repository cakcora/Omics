package v1;

import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;


import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class Initial {
    private static final Logger logger = LoggerFactory.getLogger(Initial.class);
    public static void main(String[] args) throws IOException {
        int start[] = new int[]{1,6,11,16,21,26,31,36,41, 46,51,56,61,66,72};
        int len[] = new int[]{5,5,5,5,5,5,5,5,5,5,5,5,5,6,6};
        String dataFile = args[0];
        String out = args[1];
        FileUtils.cleanDirectory(new File(out));
        logger.info(out);
        int valNA = -1;
        String line = "";
        BufferedReader br = new BufferedReader(new FileReader(dataFile));
        String header[] = br.readLine().replaceAll("\"", "").split("\t");
        int count = 0;
        while (br.readLine() != null) {
            count++;
        }
        double[][] data = new double[count][header.length];
        Map<Integer,String> proteinNames = new HashMap<Integer,String>();
        br = new BufferedReader(new FileReader(dataFile));
        br.readLine();
        int index = 0;
        int d=0;
        while ((line = br.readLine()) != null) {
            String arr[] = line.replaceAll("\"", "").split("\t");
            d=arr.length;
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
        logger.info(index+ "length "+d);
        for (int step = 3; step < 10; step=step+2) {
            double corThreshold = step/10d;
            for (int layer = 0; layer <15; layer++) {
                int edgeCount = 0;
                UndirectedSparseGraph<Integer,Integer> graph = new UndirectedSparseGraph<Integer,Integer>();
                int sIndex = start[layer];
                int dlen = len[layer];
                for (int node1 = 0; node1 < count; node1++) {
                    double[] a1 = Arrays.copyOfRange(data[node1], sIndex , sIndex + dlen);
                    if(a1.length!=dlen)logger.error("Data is not "+dlen+" length");
                    if(allNA(a1,valNA)) {
                        continue;
                    }
                    if (!graph.containsVertex(node1)) graph.addVertex(node1);
                    for (int node2 = 1 + node1; node2 < count; node2++) {
                        double[] a2 = Arrays.copyOfRange(data[node2], sIndex , sIndex + dlen);
                        if(a2.length!=dlen)logger.error("Data is not "+dlen+" length");
                        if(allNA(a2,valNA))continue;
                        if (!graph.containsVertex(node2)) graph.addVertex(node2);
                        double corr = new PearsonsCorrelation().correlation(a1, a2);
//                        logger.info(proteinNames.get(node1)+"\t"+proteinNames.get(node2)+"\t"+corr);
                        if (!Double.isNaN(corr) && corr > corThreshold) {
                            graph.addEdge(edgeCount++, node1, node2);
                        }
                    }
                }
                DescriptiveStatistics ds = new DescriptiveStatistics();
                int zeroDegreeNodeCount = 0;
                for (int v : graph.getVertices()) {
                    int neighborCount = graph.getNeighborCount(v);
                    ds.addValue(neighborCount);
                    if(neighborCount==0)zeroDegreeNodeCount++;
                }
                int edgeCount1 = graph.getEdgeCount();
                logger.info("\t"+step+"\t"+layer+"\t"+ graph.getVertexCount() + "\t" + edgeCount1 +"\t"+edgeCount1/((edgeCount1*edgeCount1-1d)/2d)+ "\t" + ds.getMean() + "\t" + ds.getPercentile(50) + "\t" + ds.getMax()+"\t"+zeroDegreeNodeCount);

                BufferedWriter wr = new BufferedWriter(new FileWriter(out+"layer"+layer+"threshold"+ step +".txt"));
                for(int s:graph.getEdges()){
                    Pair<Integer> endpoints = graph.getEndpoints(s);
                    String source =   proteinNames.get(endpoints.getFirst());
                    String target =   proteinNames.get(endpoints.getSecond());
                    wr.write(source+"\t"+target+"\r\n");
                }
                wr.close();
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
