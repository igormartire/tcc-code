import thesis.ageing.data.build.DataType;
import thesis.ageing.data.build.GOStructureType;
import weka.core.Instances;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * This class is responsible for all data related methods that might be used across the entire project.
 *
 * @author pablonsilva
 * @version 20151203
 */
public class DataUtils {
    /**
     * This method load a weka.Instances dataset given the path where the arff file is saved.
     *
     * @param path for the arff file
     * @return weka.core.Instances file with a loaded dataset.
     */
    public static Instances loadDataSet(String path) {
        Instances data = null;

        try {
            FileReader fr = new FileReader(path);
            BufferedReader br = new BufferedReader(fr);
            data = new Instances(br);

            fr.close();
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        data.setClassIndex(data.numAttributes() - 1);
        return data;
    }

    /**
     * @param data
     * @return
     */
    public static double degreeOfImbalance(Instances data) {
        double[] v = DataUtils.numInstancesPerClass(data);
        if (v[0] > v[1]) {
            return 1 - (v[1] / v[0]);
        } else {
            return 1 - (v[0] / v[1]);
        }
    }

    /**
     * @param data
     * @return
     */
    public static double[] numInstancesPerClass(Instances data) {
        double[] v = new double[2];
        v[0] = 0;
        v[1] = 0;

        for (int i = 0; i < data.numInstances(); i++) {
            if (data.instance(i).classValue() == 1)
                v[0] += 1;
            else
                v[1] += 1;
        }

        return v;
    }

    public static Instances reduce(Instances data, boolean[] subset) throws Exception {
        String indices = "";
        for (int i = 0; i < subset.length; i++) {
            if (subset[i]) {
                indices += (i + 1) + ",";
            }
        }
        indices += (data.classIndex() + 1);

        Instances tData = null;
        Remove rmv = new Remove();
        rmv.setInvertSelection(true);
        rmv.setAttributeIndices(indices);
        rmv.setInputFormat(data);

        tData = Filter.useFilter(data, rmv);
        tData.setClassIndex(tData.numAttributes() - 1);

        return tData;
    }

}
