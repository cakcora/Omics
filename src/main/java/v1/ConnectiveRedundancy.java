package v1;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

public class ConnectiveRedundancy {
    private static final Logger logger = LoggerFactory.getLogger(ConnectiveRedundancy.class);
    public static void main(String[] args) throws IOException {
        int start[] = new int[]{1,6,11,16,21,26,31,36,41, 46,51,56,61,66,72};
        int len[] = new int[]{5,5,5,5,5,5,5,5,5,5,5,5,5,6,6};
        String dataFile = args[0];
        String out = args[1];
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
        Map<Integer,String> proteinNames = new HashMap<Integer,String>();
        br = new BufferedReader(new FileReader(dataFile));
        br.readLine();
        int proteinCount = 0;
        int d=0;
        while ((line = br.readLine()) != null) {
            String arr[] = line.replaceAll("\"", "").split("\t");
            String protein = arr[0];
            proteinNames.put(proteinCount,protein);
            for (int i = 1; i < header.length; i++) {
                String v = arr[i];
                if (v.equalsIgnoreCase("NA")) {
                    data[proteinCount][i-1] = valNA;
                }
                else {
                    double value = Double.parseDouble(v);
                    data[proteinCount][i-1] = (value);
                }
            }
            proteinCount++;
        }
        logger.info(proteinCount+ " ptms");
        double corThreshold = 0.7d;
        HashMap<Integer, Graph> graphs = new HashMap<Integer, Graph>();
        for (int layer = 0; layer <start.length; layer++) {
            int edgeCount = 0;
            UndirectedSparseGraph<Integer,Integer> graph = new UndirectedSparseGraph<Integer,Integer>();
            int sIndex = start[layer]-1;
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
            graphs.put(layer,graph);
        }
        BufferedWriter wr = new BufferedWriter(new FileWriter(out+"redundancy.txt"));
        wr.write("structures.PTM\tlayerCount\tRedundancy\ttotaldegree\n");
        for (int proteinID = 0; proteinID <proteinCount; proteinID++) {
            String proteinName = proteinNames.get(proteinID);
            double totalDegree = 0d;
            int exists =0;
            HashSet<Integer> neighbors = new HashSet<Integer>();
            for(int layer = 0;layer<start.length;layer++){
                Graph graph = graphs.get(layer);
                if(graph.containsVertex(proteinID)) {
                    exists++;
                    Collection graphNeighbors = graph.getNeighbors(proteinID);
                    totalDegree += graphNeighbors.size();
                    neighbors.addAll(graphNeighbors);
                }
            }
            if(totalDegree>0){
                double redundancy = 1d- (neighbors.size()/totalDegree);
                wr.append(proteinName+"\t"+exists+"\t"+redundancy+"\t"+totalDegree+"\r\n");
            }
        }
        wr.close();
        //layer-wise redundancy
        BufferedWriter wr2 = new BufferedWriter(new FileWriter(out+"layerPairs.txt"));
        wr2.write("layer1\tlayer2\tmeanredundancy\tmedianredundancy\tstdredundancy\tptmcount\n");
        for(int l1 = 0;l1<start.length;l1++){
            Graph<Integer,Integer> graph1 = graphs.get(l1);
            for(int l2= l1+1;l2<start.length;l2++){
                DescriptiveStatistics ds = new DescriptiveStatistics();
                Graph graph2 = graphs.get(l2);
                for (int proteinID = 0; proteinID <proteinCount; proteinID++) {
                   if(!graph1.containsVertex(proteinID)||!graph2.containsVertex(proteinID))continue;
                    double degree=graph1.getNeighborCount(proteinID)+graph2.getNeighborCount(proteinID);
                    HashSet neighbors = new HashSet(graph1.getNeighbors(proteinID));
                    neighbors.addAll(graph2.getNeighbors(proteinID));
                    int unique = neighbors.size();
                    if(unique>degree)logger.error("More neighbors than degree");
                    double redundancy = 1d- (unique/degree);
                    if(unique!=0)
                   ds.addValue(redundancy);
                }
                wr2.append(l1+"\t"+l2+"\t"+ds.getMean()+"\t"+ds.getPercentile(50)+"\t"+ds.getStandardDeviation()+"\t"+ds.getN()+"\r\n");
            }

        }
        wr2.close();
    }

    private static boolean allNA(double[] a1, int valNA) {
        if((int)a1[0]==valNA&&(int)a1[1]==valNA){
            return true;
        }
        return false;
    }
}
