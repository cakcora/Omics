import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CosineSimilarity {
    private static final Logger logger = LoggerFactory.getLogger(CosineSimilarity.class);
    public double similarity(double valNA, double[] a1, double[] a2) {
        if(a1.length!=a2.length) throw new RuntimeException("Array lengths should be the same in Cosine Similarity");
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for(int i = 0;i<a1.length;i++) {

            if (a1[i] == valNA) {
                if (a2[i] != valNA) {
                    normB += Math.pow(a2[i], 2);
                }
            } else {
                if (a2[i] == valNA) {
                    normA += Math.pow(a1[i], 2);
                } else {
                    dotProduct += a1[i] * a2[i];
                    normA += Math.pow(a1[i], 2);
                    normB += Math.pow(a2[i], 2);
                }
            }
        }
        double sim =0d;
        if(dotProduct>0d) {
           sim = dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
            return sim;
        }
        return sim;
    }

    public double similarity(double[] a1, double[] a2) {
        if(a1.length!=a2.length) throw new RuntimeException("Array lengths should be the same in Cosine Similarity");
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for(int i = 0;i<a1.length;i++) {

            if (Double.isNaN(a1[i] )) {
                if (!Double.isNaN(a2[i])) {
                    normB += Math.pow(a2[i], 2);
                }
            } else {
                if (Double.isNaN(a2[i])) {
                    normA += Math.pow(a1[i], 2);
                } else {
                    dotProduct += a1[i] * a2[i];
                    normA += Math.pow(a1[i], 2);
                    normB += Math.pow(a2[i], 2);
                }
            }
        }
        double sim =0d;
        if(dotProduct>0d) {
            sim = dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
            if(sim>1.05||sim<0.0) logger.error("Cosine similarity value is invalid: "+sim);
        }
        return sim;
    }
}
