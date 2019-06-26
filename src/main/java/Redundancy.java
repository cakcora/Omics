import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import structures.Dataset;
import structures.Row;

import java.io.*;
import java.util.*;


public class Redundancy {
    private static final Logger logger = LoggerFactory.getLogger(Redundancy.class);
    static HashMap<String, Integer> names = new HashMap();
    static HashMap<Integer,String> ids = new HashMap();
    public static void main(String[] args) throws IOException {

        if(false)perExperiment(args);

        String markFile = args[6];
        String outDir = args[1];
        String redFile = args[4];
        String kgDataFile = args[7];
        HashMap<Integer, Graph> graphs = new HashMap<Integer, Graph>();
        graphs.put(0,createCFNGraph(markFile,1,2));
        graphs.put(1,createCFNGraph(kgDataFile,1,2));
        String layerRedFile = args[5];
        BufferedWriter wr = new BufferedWriter(new FileWriter(outDir+redFile));
        multipleRedundancy(graphs, wr, "CFNs", 0, 1);
        wr.close();


    }

    private static void perExperiment(String[] args) throws IOException {
        String markFile = args[0];
        String kgDataFile = args[3];
        String outDir = args[1];
        double simThreshold = Double.parseDouble(args[2]);
        String redFile = args[4];
        String layerRedFile = args[5];

        Dataset[] sets = Loader.getDatasets(markFile, kgDataFile);

        CosineSimilarity simFunction = new CosineSimilarity();
        HashMap<Integer, Graph> graphs = createGraphs(simFunction,simThreshold, sets);
        BufferedWriter wr = new BufferedWriter(new FileWriter(outDir+redFile));


        multipleRedundancy(graphs, wr, "M1to9", 0, 8);
        multipleRedundancy(graphs, wr, "M10to15", 9, 14);
        multipleRedundancy(graphs, wr, "KG1to6", 15, 18);
        wr.close();
        if(false) {
            HashMap<Integer, int[]> redVals = layerwiseRedundancy(outDir, layerRedFile, graphs);
            printRedMatrix(outDir, layerRedFile, redVals);
        }
        // new experiment fro crizotinib
        Dataset set1 = new PTMLoader().parse(kgDataFile, 1, 3);
        Dataset set2 = new PTMLoader().parse(markFile, 62, 2);
        Dataset[] sets2 = {set1, set2};
        HashMap<Integer, Graph> graphs2 = createGraphs(simFunction,simThreshold, sets2);
        HashMap<Integer, int[]> redVals = layerwiseRedundancy(outDir, "2"+layerRedFile,  graphs2);
        printRedMatrix(outDir, "2"+layerRedFile, redVals);
    }

    private static void multipleRedundancy(HashMap<Integer, Graph> graphs, BufferedWriter wr, String setName, int start, int end) throws IOException {
        HashMap<Integer, ArrayList<Integer>> neighbors = new HashMap<>();
        HashMap<Integer, Integer> exists = new HashMap<>();
        for(int layer1 = start; layer1<=end; layer1++){
            Graph graph = graphs.get(layer1);
            for(Object n: graph.getVertices()){
                int node = (int) n;
                if(!neighbors.containsKey(node)){
                    neighbors.put(node,new ArrayList<>());
                    exists.put(node,0);
                }
                neighbors.get(node).addAll(graph.getNeighbors(n));
                exists.put(node,1+exists.get(node));
            }
        }

        for(int n:exists.keySet()){
            if(exists.get(n)>1){

                ArrayList<Integer> ints = neighbors.get(n);
                int size = ints.size();
                int size1 = new HashSet<Integer>(ints).size();
                double red = 100*(1.0- size1 /(double) size);

                wr.append(ids.get(n)+"\t"+ setName +"\t"+exists.get(n)+"\t"+red+"\t"+ size1 +"\r\n");
            }
        }
    }

    private static HashMap<Integer, int[]> layerwiseRedundancy(String out, String layerRedFile, HashMap<Integer, Graph> graphs) throws IOException {
        //layer-wise redundancy
        BufferedWriter wr2 = new BufferedWriter(new FileWriter(out+ layerRedFile));
        wr2.write("layer1\tlayer2\tproteinCount\tmeanredundancy\tmedianredundancy\tstdredundancy\tmaxredundancy\tptmcount\n");
        HashMap<Integer,int[]> redVals = new HashMap<>();
        int numOfLayers = graphs.size();
        int currentLevel=0;
        for(int layer1 = 0; layer1< numOfLayers; layer1++){

            Graph<Integer,Integer> graph1 = graphs.get(layer1);
            for(int layer2 = layer1+1; layer2< numOfLayers; layer2++){

                int howMany =0;
                DescriptiveStatistics ds = new DescriptiveStatistics();
                Graph graph2 = graphs.get(layer2);
                for (int node :graph1.getVertices()) {
                    if(!graph2.containsVertex(node)||graph1.getNeighborCount(node)==0||graph2.getNeighborCount(node)==0)continue;
                    else {
                        double degree = graph1.getNeighborCount(node) + graph2.getNeighborCount(node);
                        howMany++;
                        HashSet neighbors = new HashSet(graph1.getNeighbors(node));
                        neighbors.addAll(graph2.getNeighbors(node));
                        int unique = neighbors.size();
                        if (unique > degree) logger.error("More neighbors than degree");
                        double redundancy = 1d - (unique / degree);
                        if (unique != 0)
                            ds.addValue(redundancy);

                        if(!redVals.containsKey(node)){
                            int [] arr = new int[(int)Math.ceil(numOfLayers*(numOfLayers-1.0)/2.0)];
                            redVals.put(node,arr);
                        }
                        redVals.get(node)[currentLevel]= (int) (100*redundancy);
                    }

                }

                HashSet all = new HashSet();
                all.addAll(graph1.getVertices());
                all.addAll(graph2.getVertices());
                wr2.append(layer1+"\t"+layer2+"\t"+howMany+"\t"+ds.getMean()+"\t"+ds.getPercentile(50)+"\t"+ds.getStandardDeviation()+"\t"+ds.getMax()+"\t"+all.size()+"\r\n");
                currentLevel++;
            }
        }
        wr2.close();
        return redVals;
    }

    private static void printRedMatrix(String out, String file, HashMap<Integer, int[]> redVals) throws IOException {

        BufferedWriter wr3 = new BufferedWriter(new FileWriter(out+ file));
        for(int node:redVals.keySet()){
            String name = ids.get(node);
            wr3.append(name);
            int[] arr = redVals.get(node);
            for(long a: arr){
                wr3.append("\t"+a);
            }
            wr3.append("\r\n");
        }
        wr3.close();
    }

    private static HashMap<Integer, Graph> createGraphs(CosineSimilarity similarity, double thres, Dataset[] sets) {
        HashMap<Integer, Graph> graphs = new HashMap<Integer, Graph>();

        for(int ind =0;ind<sets.length;ind++){
            Dataset d =sets[ind];
            d.removeAllNans();
            UndirectedSparseGraph<Integer,Integer> graph = new UndirectedSparseGraph<Integer,Integer>();

            int edgeCount = 0;
            ArrayList<Row> rows = d.getRows();
            for (int i = 0; i < rows.size(); i++) {
                Row r1= rows.get(i);
                int r1ID = getID(names,r1.getPTMName());

                double []d1 = r1.getData();
                for (int i2 = 0; i2 < rows.size(); i2++){
                    Row r2= rows.get(i2);
                    double []d2 = r2.getData();
                    double sim = similarity.similarity(d1, d2);
                    if(sim>thres){
                        int r2ID = getID(names,r2.getPTMName());
                        graph.addEdge(edgeCount++, r1ID, r2ID);
                    }
                }
            }
            logger.info("Graph details: "+graph.getVertexCount()+" vertices, "+graph.getEdgeCount()+" edges");
            graphs.put(ind,graph);
        }
        return graphs;
    }
    static Graph<Integer,Integer> createCFNGraph(String input, int i, int o) throws IOException {
        HashMap<Integer, Graph> graphs = new HashMap<Integer, Graph>();
        String line ="";
        String sep ="\t";
        BufferedReader dataFile = new BufferedReader(new FileReader(input));
        String header = dataFile.readLine();
        int edgeCount = 0;
        UndirectedSparseGraph<Integer,Integer> graph = new UndirectedSparseGraph<Integer,Integer>();
        while((line=dataFile.readLine())!=null){
            String arr[] = line.split(sep);
            String ptm1 = arr[i];
            String ptm2 = arr[o];
            int r1ID = getID(names,ptm1);
            int r2ID = getID(names,ptm2);
            graph.addEdge(edgeCount++, r1ID, r2ID);

        }
        logger.info("Graph details: "+graph.getVertexCount()+" vertices, "+graph.getEdgeCount()+" edges");
        return graph;

    }

    private static int getID(HashMap<String,Integer> names, String r1) {
        int r1ID;
        if(names.containsKey(r1)){
            r1ID = names.get(r1);
        }
        else {
            names.put(r1, names.size());
            r1ID = names.get(r1);
        }
        ids.put(r1ID,r1);
        return r1ID;
    }
}
