import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import structures.Dataset;
import structures.Row;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class Deviation {
    private static final Logger logger = LoggerFactory.getLogger(Deviation.class);
    static HashMap<String, Integer> names = new HashMap();
    static HashMap<Integer,String> ids = new HashMap();
    public static void main(String[] args) throws IOException {

        if(false)perExperiment(args);
        String markFile = args[5];
        String outDir = args[1];
        String redFile = args[4];
        String kgDataFile = args[6];
        HashMap<Integer, Graph> graphs = new HashMap<Integer, Graph>();
        Redundancy red = new Redundancy();
        graphs.put(0,red.createCFNGraph(markFile,1,2));
        graphs.put(1,red.createCFNGraph(kgDataFile,1,2));
        String layerRedFile = args[5];
        BufferedWriter wr = new BufferedWriter(new FileWriter(outDir+redFile));
        HashMap<Integer,int[]> degVals = new HashMap<>();
        int numOfLayers = graphs.size();
        ids = Redundancy.ids;
        for(int layer1 = 0; layer1< numOfLayers; layer1++){
            Graph<Integer,Integer> graph1 = graphs.get(layer1);
            for (int node :graph1.getVertices()) {
                int degree = graph1.getNeighborCount(node)  ;
                if(!degVals.containsKey(node)){
                    degVals.put(node,new int[numOfLayers]);
                }
                degVals.get(node)[layer1] = degree;
            }
        }
        String out = args[1];
        String degreeFile = args[4];
        printMatrix(out,degreeFile,degVals);
        wr.close();
    }

    private static void perExperiment(String[] args) throws IOException {
        String dataFile = args[0];
        String out = args[1];
        double simThreshold = Double.parseDouble(args[2]);
        String kgDataFile = args[3];
        String degreeFile = args[4];

        Dataset[] sets = Loader.getDatasets(dataFile, kgDataFile);

        CosineSimilarity simFunction = new CosineSimilarity();
        HashMap<Integer, Graph> graphs = createGraphs(simFunction,simThreshold, sets);

        //layer-wise redundancy

        HashMap<Integer,int[]> degVals = new HashMap<>();
        int numOfLayers = sets.length;
        for(int layer1 = 0; layer1< numOfLayers; layer1++){
            Graph<Integer,Integer> graph1 = graphs.get(layer1);
                for (int node :graph1.getVertices()) {
                        int degree = graph1.getNeighborCount(node)  ;
                        if(!degVals.containsKey(node)){
                            degVals.put(node,new int[numOfLayers]);
                        }
                        degVals.get(node)[layer1] = degree;
                }
        }

        printMatrix(out,degreeFile,degVals);
        return;
    }

    private static void printMatrix(String out, String file, HashMap<Integer, int[]> degrees) throws IOException {

        BufferedWriter wr3 = new BufferedWriter(new FileWriter(out+ file));
        for(int node:degrees.keySet()){
            String name = ids.get(node);
            wr3.append(name);
            int[] arr = degrees.get(node);
            for(long a: arr){
                wr3.append(","+a);
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
