import structures.Dataset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Loader {
    static Dataset[] getDatasets(String dataFile, String kgDataFile) throws IOException {
        Dataset set1 = new PTMLoader().parse(dataFile, 1, 5);
        Dataset set2 = new PTMLoader().parse(dataFile, 6, 5);
        Dataset set3 = new PTMLoader().parse(dataFile, 11, 5);
        Dataset set4 = new PTMLoader().parse(dataFile, 16, 5);
        Dataset set5 = new PTMLoader().parse(dataFile, 21, 5);
        Dataset set6 = new PTMLoader().parse(dataFile, 26, 5);
        Dataset set7 = new PTMLoader().parse(dataFile, 31, 5);
        Dataset set8 = new PTMLoader().parse(dataFile, 36, 5);
        Dataset set9 = new PTMLoader().parse(dataFile, 41, 5);

        Dataset set10 = new PTMLoader().parse(dataFile, 46, 5);
        Dataset set11 = new PTMLoader().parse(dataFile, 51, 5);
        Dataset set12 = new PTMLoader().parse(dataFile, 56, 5);
        Dataset set13 = new PTMLoader().parse(dataFile, 61, 5);
        Dataset set14 = new PTMLoader().parse(dataFile, 66, 6);
        Dataset set15 = new PTMLoader().parse(dataFile, 72, 6);

        Dataset set16 = new PTMLoader().parse(kgDataFile, 1, 3);
        Dataset set17 = new PTMLoader().parse(kgDataFile, 4, 3);
        Dataset set18 = new PTMLoader().parse(kgDataFile, 7, 3);
        Dataset set19 = new PTMLoader().parse(kgDataFile, 10, 3);

        Dataset[] sets = new Dataset[19];
        sets[0] = set1;
        sets[1] = set2;
        sets[2] = set3;
        sets[3] = set4;
        sets[4] = set5;
        sets[5] = set6;
        sets[6] = set7;
        sets[7] = set8;
        sets[8] = set9;
        sets[9] = set10;
        sets[10] = set11;
        sets[11] = set12;
        sets[12] = set13;
        sets[13] = set14;
        sets[14] = set15;
        sets[15] = set16;
        sets[16] = set17;
        sets[17] = set18;
        sets[18] = set19;
        return sets;
    }



}