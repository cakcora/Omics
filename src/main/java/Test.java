import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;

import static java.util.Comparator.comparing;

public class Test {
    private static final Logger logger = LoggerFactory.getLogger(Test.class);
    public static void main(String[] args) throws IOException {
        double[] a1 = {-1d, 3, 2};
        double[] a2 = {3, -1d, 2};
        double sim = new CosineSimilarity().similarity(-1d,a1,a2);
        logger.info(sim+" ");

    }


}
