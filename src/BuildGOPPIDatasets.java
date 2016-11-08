import weka.core.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class BuildGOPPIDatasets {

    //arquivo contendo todos os pares geneid (ou entrezid) e string_id
    private String _path_geneid_stringid;

    //pasta contendo os arquivos com as interacoes entre proteinas de cada organismo:
    // ex de arquivo: CE.protein.links.detailed.v10.txt
    private String _pathPPIData;

    // pasta contendo todas as bases com GO
    private String _pathGoDatasets;

    // pasta de saida
    private String _path_output;

    private HashMap<String, String> string_gene = new HashMap<>();

    public BuildGOPPIDatasets(String go, String gene_string, String ppiData, String output) {
        _path_geneid_stringid = gene_string;
        _pathPPIData = ppiData;
        _pathGoDatasets = go;
        _path_output = output;
    }

    public static void main(String[] args) {
        BuildGOPPIDatasets b = new BuildGOPPIDatasets(
                "/Users/pablonsilva/Desktop/Base-gerada-nova/Threshold/3/",
                "/Users/pablonsilva/Google Drive/Doutorado/Biology of Ageing/Data/PPI/PPIs - STRING/entrez_string.txt",
                "/Users/pablonsilva/Google Drive/Doutorado/Biology of Ageing/Data/PPI/PPIs - STRING/",
                "/Users/pablonsilva/Desktop/nova-base/");
        b.buildData(400);
        b.buildData(500);
        b.buildData(600);
        b.buildData(700);
        b.buildData(800);
        b.buildData(900);
    }

    public void buildData(double score_min) {

        //MAP gene_id (entrezid) with stringID
        try {
            BufferedReader br = new BufferedReader(new FileReader(_path_geneid_stringid));
            String line = null;
            br.readLine();
            do {
                line = br.readLine();
                if (line != null) {
                    String[] splits = line.split(",");
                    String geneId = splits[0];
                    String stringId = splits[1];

                    string_gene.put(stringId, geneId);
                }
            } while (line != null);
            br.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

        //Pega todas as interacoes de cada organismo com os genes da base GO
        HashMap<String, ArrayList<String>> map_geneId_interactions_ce = getInteractions("CE", score_min);
        HashMap<String, ArrayList<String>> map_geneId_interactions_dm = getInteractions("DM", score_min);
        HashMap<String, ArrayList<String>> map_geneId_interactions_mm = getInteractions("MM", score_min);
        HashMap<String, ArrayList<String>> map_geneId_interactions_sc = getInteractions("SC", score_min);

        //GEra as bases GO+PPI para cada organismo e estrutura hierarquica
        String[] var = {"BP", "CC", "MF", "BP+CC", "BP+MF", "CC+MF", "BP+CC+MF"};

        for (int j = 0; j < var.length; j++) {
            //Gera bases
            Instances d_ce = buildArff(DataUtils.loadDataSet(_pathGoDatasets + "CE-" + var[j] + "-threshold-3.arff"), map_geneId_interactions_ce);
            Instances d_dm = buildArff(DataUtils.loadDataSet(_pathGoDatasets + "DM-" + var[j] + "-threshold-3.arff"), map_geneId_interactions_dm);
            Instances d_mm = buildArff(DataUtils.loadDataSet(_pathGoDatasets + "MM-" + var[j] + "-threshold-3.arff"), map_geneId_interactions_mm);
            Instances d_sc = buildArff(DataUtils.loadDataSet(_pathGoDatasets + "SC-" + var[j] + "-threshold-3.arff"), map_geneId_interactions_sc);

            //Salva bases
            FileUtils.saveFile(d_ce.toString(), _path_output + (int) score_min + "/GO_PPI_CE-" + var[j] + ".arff");
            FileUtils.saveFile(d_dm.toString(), _path_output + (int) score_min + "/GO_PPI_DM-" + var[j] + ".arff");
            FileUtils.saveFile(d_mm.toString(), _path_output + (int) score_min + "/GO_PPI_MM-" + var[j] + ".arff");
            FileUtils.saveFile(d_sc.toString(), _path_output + (int) score_min + "/GO_PPI_SC-" + var[j] + ".arff");

            //Aplica o threshold = 3, ou seja cada termo deve aparecer em no minimo 3 instancias (genes)
            d_ce = thresholdFeatures(d_ce, 3);
            d_dm = thresholdFeatures(d_dm, 3);
            d_mm = thresholdFeatures(d_mm, 3);
            d_sc = thresholdFeatures(d_sc, 3);

            //salva as bases com threshold
            FileUtils.saveFile(d_ce.toString(), _path_output + (int) score_min + "/Threshold/GO_PPI_CE-" + var[j] + "-threshold-3.arff");
            FileUtils.saveFile(d_dm.toString(), _path_output + (int) score_min + "/Threshold/GO_PPI_DM-" + var[j] + "-threshold-3.arff");
            FileUtils.saveFile(d_mm.toString(), _path_output + (int) score_min + "/Threshold/GO_PPI_MM-" + var[j] + "-threshold-3.arff");
            FileUtils.saveFile(d_sc.toString(), _path_output + (int) score_min + "/Threshold/GO_PPI_SC-" + var[j] + "-threshold-3.arff");
        }
    }

    /***
     * Gera a base de dados em formato arff a partir da base GO e das interacoes das proteinas nessa base.
     *
     * @param dataGO base GO
     * @param map    interacoes de cada proteina (gene) nessa base com as outras proteinas na base string.
     * @return
     */
    private Instances buildArff(Instances dataGO, HashMap<String, ArrayList<String>> map) {

        //Monta a estrutura da base de dados
        FastVector vector = new FastVector();

        // All GO terms are binary.
        FastVector vals = new FastVector();
        vals.addElement("0");
        vals.addElement("1");

        // adiciona um atributo string que representa o entrezid (geneID)
        vector.addElement(new Attribute("entrez", (FastVector) null));

        // copia a estrutura da base GO
        for (int i = 1; i < dataGO.numAttributes() - 1; i++) {
            String go = dataGO.attribute(i).name();
            vector.addElement(new Attribute(go, vals));
        }

        // configura os atributos da ppi.
        ArrayList<String> ppis = new ArrayList<String>();
        for (String geneId : map.keySet()) {
            ArrayList<String> features = map.get(geneId);
            for (String ppi : features) {
                if (!ppis.contains(ppi)) {
                    ppis.add(ppi);
                    vector.addElement(new Attribute(ppi, vals));
                }
            }
        }

        // atributo classe
        vector.addElement(new Attribute("class", vals));

        //Preenche as instancias

        //cria uma base a partir da estrutura montada acima
        Instances data = new Instances("GO_PPI_" + dataGO.relationName(), vector, 0);

        // Adiciona a informacao dos PPIs na base GO
        for (Instance inst : dataGO) {
            double[] values = new double[vector.size()];

            String geneId = inst.stringValue(0);

            //entrez id
            values[0] = data.attribute(0).addStringValue(geneId);
            //adiciona os GOs
            for (int i = 1; i < dataGO.numAttributes() - 1; i++) {
                values[i] = inst.value(i);
            }

            ArrayList<String> inst_interactions = map.get(geneId);

            //se instancia tem PPI - adiciona.
            if (inst_interactions != null) {
                for (String interaction : inst_interactions) {
                    int index = data.attribute(interaction).index();
                    values[index] = 1;
                }
            } // senao os valores dos ppis para essa instancia serao missing-values (?).
            else {
                for (int i = dataGO.numAttributes() - 1; i < values.length - 1; i++) {
                    values[i] = Double.NaN;
                }
            }

            // repete a classe
            values[values.length - 1] = inst.classValue();

            // adiciona a instancia na base final
            data.add(new DenseInstance(1.0, values));
        }

        return data;
    }

    /**
     * Numero de vezes que cada GO term aparece na base.
     *
     * @param data
     * @return
     */
    private int[] numGOTerms(Instances data) {
        int[] numTotal = new int[data.numAttributes()];

        for (int i = 0; i < data.numInstances(); i++) {
            double[] vec = data.instance(i).toDoubleArray();
            for (int j = 0; j < vec.length; j++)
                if (!Double.isNaN(vec[j]))
                    numTotal[j] += vec[j];
        }

        return numTotal;
    }

    /***
     * Aplica o threshold
     *
     * @param data
     * @param threshold
     * @return
     */
    private Instances thresholdFeatures(Instances data, int threshold) {
        Instances dataAux = new Instances(data);
        dataAux.setClassIndex(data.numAttributes() - 1);
        int[] a = numGOTerms(data);

        for (int i = data.numAttributes() - 2; i > 0; i--) {
            System.out.println(data.attribute(i).name() + " - " + a[i]);
            if (a[i] < threshold) {
                System.out.println("Removed");
                dataAux.deleteAttributeAt(i);
            }
        }


        // Se uma instancia nao tiver PPIs deve-se configura os atributos como missing-values
        for (int j = data.numInstances() - 1; j >= 0; j--) {
            int numPPIs = 0;
            for (int i = data.numAttributes() - 2; i > 0; i--) {
                if (!data.attribute(i).name().contains("GO:")) {
                    numPPIs++;
                }
            }
            if (numPPIs == 0) {
                for (int i = 1; i < data.numAttributes() - 2; i++) {
                    if (!dataAux.attribute(i).name().contains("GO:"))
                        dataAux.instance(j).setValue(i, Double.NaN);
                }
            }
        }

        return dataAux;
    }

    /**
     * Obtem todas as interacoes das proteinas de um determinado organismo
     *
     * @param model     organismo
     * @param score_min score minimo
     * @return
     */
    private HashMap<String, ArrayList<String>> getInteractions(String model, double score_min) {

        HashMap<String, ArrayList<String>> map_geneId_interactions = new HashMap<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(_pathPPIData + model + ".protein.links.detailed.v10.txt"));
            String line = null;
            br.readLine();
            do {
                line = br.readLine();
                if (line != null) {
                    String[] splits = line.split(" ");
                    String protein1 = splits[0];
                    String protein2 = splits[1];
                    double score = Double.parseDouble(splits[splits.length - 1]);

                    // obs: esse score eh utilizado para montar a sua base.

                    if (score >= score_min) {
                        String geneId = string_gene.get(protein1);

                        if (geneId != null) {

                            ArrayList<String> ppi_list = map_geneId_interactions.get(geneId);
                            if (ppi_list == null)
                                ppi_list = new ArrayList<>();

                            if (!ppi_list.contains(protein2))
                                ppi_list.add(protein2);

                            //System.out.println(protein1 + " - " + protein2);

                            map_geneId_interactions.put(geneId, ppi_list);
                        }
                    }
                }
            } while (line != null);
            br.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

        return map_geneId_interactions;
    }
}
