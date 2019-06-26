import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

import static java.util.Comparator.comparing;

public class FlattenedConnectiveRedundancy {
    private static final Logger logger = LoggerFactory.getLogger(FlattenedConnectiveRedundancy.class);
    public static void main(String[] args) throws IOException {
//          int start[] = new int[]{1,46,56,66};
//         int len[] = new int[]{45,10,10,12};
        int start[] = new int[]{1,6,11,16,21,26,31,36,41,46,51,56,61,66,72};
        int len[] = new int[]{5,5,5,5,5,5,5,5,5,5,5,5,5,6,6};


        String dataFile = args[0];
        String out = args[1];
        double corThreshold = 0.8d;
        String line = "";
        BufferedReader br = new BufferedReader(new FileReader(dataFile));
        String header[] = br.readLine().replaceAll("\"", "").split("\t");
        int count = 0;
        double valNA = -1d;
        int d=0;
        while (br.readLine() != null) {
            count++;
        }
        int topK = 500000;
        logger.info(header.length+" columns "+topK+"\t"+corThreshold);
        double[][] data = new double[count][header.length];
        Map<Integer,String> proteinNames = new HashMap<Integer,String>();
        br = new BufferedReader(new FileReader(dataFile));
        br.readLine();
        int proteinCount = 0;
        d=0;
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

        HashMap<Integer, Graph> graphs = new HashMap<Integer, Graph>();
        for (int layer = 0; layer <start.length; layer++) {
            int edgeCount = 0;
            UndirectedSparseGraph<Integer,Integer> graph = new UndirectedSparseGraph<Integer,Integer>();
            int sIndex = start[layer]-1;
            int dlen = len[layer];


            for (int j = 0; j < count; j++) {
                double[] a1 = Arrays.copyOfRange(data[j], sIndex , sIndex + dlen);

                if (!graph.containsVertex(j)) graph.addVertex(j);
                double[] valQ = new double [count];
                for (int k = 1 + j; k < count; k++) {
                    double[] a2 = Arrays.copyOfRange(data[k], sIndex , sIndex + dlen);
                    if (!graph.containsVertex(k)) graph.addVertex(k);
                    if(a1.length!=dlen)logger.error(layer+"th layer has "+a1.length+", but we need "+dlen);
                    double sim = new CosineSimilarity().similarity(valNA,a1, a2);

                    if (!Double.isNaN(sim)) {
                        valQ[k]=sim;
                    }
                }

                for(int z:topN(valQ, topK)){
                    if(valQ[z]>corThreshold){
                        graph.addEdge(edgeCount++, j, z);
                    }
                }
            }
            logger.info(layer+"\t"+graph.getVertexCount()+"\t"+graph.getEdgeCount());
            TreeMap<Integer,Integer> degrees = new TreeMap<Integer,Integer>();
            for(int ver:graph.getVertices()){
                int deg = graph.getNeighborCount(ver);
                if(!degrees.containsKey(deg)) degrees.put(deg,1);
                else {
                    degrees.put(deg,1+ degrees.get(deg));
                }
            }
            BufferedWriter wr3 = new BufferedWriter(new FileWriter(layer+"degreeDist.txt"));
            int max = degrees.lastKey();
            for(int de=0;de<=max;de++){

                Integer integer = 0;
                if(degrees.containsKey(de)) {
                    integer = degrees.get(de);
                }
                wr3.append(de+"\t"+ integer +"\r\n");
            }
            wr3.close();
            graphs.put(layer,graph);
        }
        //layer-wise redundancy
        BufferedWriter wr2 = new BufferedWriter(new FileWriter(out+"FlattenedRedundancy.txt"));
        wr2.write("layer1\tlayer2\tproteinCount\tmeanredundancy\tmedianredundancy\tstdredundancy\tptmcount\n");
        for(int l1 = 0;l1<start.length;l1++){

            Graph<Integer,Integer> graph1 = graphs.get(l1);
            for(int l2= l1+1;l2<start.length;l2++){
            int howMany =0;
                DescriptiveStatistics ds = new DescriptiveStatistics();
                Graph graph2 = graphs.get(l2);
                for (int proteinID = 0; proteinID <proteinCount; proteinID++) {
                   if(graph1.getNeighborCount(proteinID)==0||graph2.getNeighborCount(proteinID)==0)continue;
                   else {
                       double degree = graph1.getNeighborCount(proteinID) + graph2.getNeighborCount(proteinID);
                       howMany++;
                       HashSet neighbors = new HashSet(graph1.getNeighbors(proteinID));
                       neighbors.addAll(graph2.getNeighbors(proteinID));
                       int unique = neighbors.size();
                       if (unique > degree) logger.error("More neighbors than degree");
                       double redundancy = 1d - (unique / degree);
                       if (unique != 0)
                           ds.addValue(redundancy);
                   }

                }
                wr2.append(l1+"\t"+l2+"\t"+ds.getMean()+"\t"+howMany+"\t"+ds.getPercentile(50)+"\t"+ds.getStandardDeviation()+"\t"+ds.getN()+"\r\n");
            }

        }
        wr2.close();

    }
    public static int[] topN(final double[] input, final int n) {
        return IntStream.range(0, input.length)
                .boxed()
                .sorted(comparing(i -> -input[i]))
                .mapToInt(i -> i)
                .limit(n)
                .toArray();
    }

    private static boolean allNA(double[] a1, int valNA) {
        if((int)a1[0]==valNA&&(int)a1[1]==valNA){
            return true;
        }
        return false;
    }
}
