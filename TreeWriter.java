import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 18/07/2013
 * Time: 16:33
 * To change this template use File | Settings | File Templates.
 */
public class TreeWriter {

    int n_lineages = 200;
    List<String> names = new ArrayList<String>();
    List<Double> sampleTimes = new ArrayList<Double>();
    int[][] coalEvents;
    double[] nodeTimes;

    public void setParams(String nodeTimes_filename, String coalEvent_filename, String names_filename) {

        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(nodeTimes_filename));
            BufferedReader reader2 = new BufferedReader(new FileReader(names_filename));
            BufferedReader reader3 = new BufferedReader(new FileReader(coalEvent_filename));

            nodeTimes = new double[(n_lineages*2)-1];
            coalEvents = new int[(n_lineages-1)][2];

            int i = 0;
            while(reader1.ready() && i < nodeTimes.length) {

                String line = reader1.readLine();
                line = line.trim();
                Double nodeTime = Double.parseDouble(line);
                nodeTimes[i] = nodeTime;
                i++;

            }

            int k = 0;
            while(reader2.ready()) {

                String line = reader2.readLine();
                line = line.trim();
                names.add(line);
                if(k < n_lineages) {

                    String [] parts = line.split("_");
                    String sampleTime_string = parts[parts.length-1].trim();
                    sampleTime_string = sampleTime_string.replaceAll("\"","");
                    sampleTimes.add(Double.parseDouble(sampleTime_string));
                    k++;

                }

            }

            int j = 0;
            while(reader3.ready() && j < coalEvents.length) {

                String line = reader3.readLine();
                line = line.trim();
                String[] parts = line.split(",");
                coalEvents[j][0] = Integer.parseInt(parts[0].trim());
                coalEvents[j][1] = Integer.parseInt(parts[1].trim());
                j++;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }
    public void writeNexusTree(int[][] coalescentEvents, double[] nodeTimes, List<String> names, String treeFileName) {
        //public void writeNexusTree(String nodeTime_file, String names_file, String coalescentEvents_file, String treeFile) {


        Map<Integer, String> nodeNames = new HashMap<Integer, String>();
        Map<Integer, Double> nodeHeights = new HashMap<Integer, Double>();
        fillNodeMaps(nodeNames, nodeHeights);

        FileWriter writer4 = null;
        //File treeOutputFile = new File("tree_"+n_lineages+"N_"+params.N+"_antigenicMu_"+params.antigenMu+"_s_"+params.dfe+"_psi_"+params.psi+"_psij_"+params.psi_j+"_samplePeriod_"+sampleStartTime/365.25+"_to_"+sampleEndTime/365.25+"_simNo_"+simNo+"_segment_"+(i+1)+".tre");

        try {
            writer4 = new FileWriter(new File(treeFileName));


            writer4.write("#NEXUS\n");
            writer4.write("begin taxa;\n");
            writer4.write("\tdimensions ntax="+n_lineages+";\n");
            writer4.write("\ttaxlabels\n");

            int n = 0;
            while(n < n_lineages) {
                writer4.write("\t'"+names.get(n)+"'\n");
                n++;
            }
            writer4.write(";\n");
            writer4.write("end;\n");
            writer4.write("begin trees;\n");
            int k = 0;
            String treeString = "";

            while(k < coalescentEvents.length) {

                int child1 = coalescentEvents[k][0];
                int child2 = coalescentEvents[k][1];


                String node1 = "";
                String node2 = "";

                //nodeNames should contain the children names! Just a check that everything is working as it should be!
                if(nodeNames.containsKey(child1)) {

                    node1 = nodeNames.get(child1);
                }
                if(nodeNames.containsKey(child2)) {

                    node2 = nodeNames.get(child2);
                }

                String[] partsNode1 = node1.split("_");
                String[] partsNode2 = node2.split("_");
                String[] partsParent = names.get(n_lineages+k).split("_");

                String patchInfo = "?";
                String reassInfo = "?";

                String patchNode1 = "?";
                String patchNode2 = "?";
                String reassNode1 = "?";
                String reassNode2 = "?";

                if(partsNode1.length > 2 && partsNode2.length > 2) {
                    patchNode1 = partsNode1[1];
                    patchNode2 = partsNode2[1];
                    reassNode1 = partsNode1[2];
                    reassNode2 = partsNode2[2];
                }

                if(partsParent.length > 5)  {
                    patchInfo = partsParent[3];
                    reassInfo = partsParent[5];
                }


                treeString = "("+node1+"[&patch=\""+patchNode1+"\",reassortant=\""+reassNode1+"\"]:"+nodeTimes[child1-1]+
                        ","+node2+"[&patch=\""+patchNode2+"\",reassortant=\""+reassNode2+"\"]:"+nodeTimes[child2-1]
                        +")[&parent="+names.get(n_lineages+k)+",parentReass=\""+reassInfo+"\"]";

                //System.out.println(nodeString);

                nodeNames.remove(child1);
                nodeNames.remove(child2);

                nodeNames.put(n_lineages+k+1, treeString);
                nodeHeights.put(n_lineages+k+1,nodeTimes[n_lineages+k-1]);

                k++;
            }

            //add the root node info:
            //treeString = treeString+"[&parent='root'"

            System.out.println("******");
            System.out.println(treeString);
            writer4.write("\t tree TREE1 = [&R] [&R=true]"+treeString+";\n");
            writer4.write("end;\n");
            writer4.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    private void fillNodeMaps(Map<Integer, String> nodeNames, Map<Integer, Double> nodeHeights){
        //, List<Integer> sampledLineages, List<Double> sampledTimes) {

        for(int i=0; i<n_lineages; i++) {
            nodeNames.put(i+1, "'"+names.get(i)+"'");
            nodeHeights.put(i+1, sampleTimes.get(i));

        }
    }

    public static void main(String [] args) {

        TreeWriter wTrees = new TreeWriter();
        wTrees.setParams
                ("/Users/jayna/Documents/MATLAB/Duke_Reassortment/2013-07-19/varying_s/psi_10.0_antigenicMu_0.005/s_0.05/nodeTimes_endemic_200N_7000000_antigenicMu_0.0050_s_0.05_psi_1.0_psij_1.0_samplePeriod_39.0_to_40.0_simNo_0_segment_2.csv",
                 "/Users/jayna/Documents/MATLAB/Duke_Reassortment/2013-07-19/varying_s/psi_10.0_antigenicMu_0.005/s_0.05/coalescentPairs_endemic_200_N_7000000_antigenicMu_0.0050_s_0.05_psi_1.0_psij_1.0_samplePeriod_39.0_to_40.0_simNo_0_segment_2.csv",
                 "/Users/jayna/Documents/MATLAB/Duke_Reassortment/2013-07-19/varying_s/psi_10.0_antigenicMu_0.005/s_0.05/names_endemic_200N_7000000_antigenicMu_0.0050_s_0.05_psi_10.0_psij_1.0_samplePeriod_39.0_to_40.0_simNo_0_segment_2.csv");
        //(args[1], args[2], args[3]);
        wTrees.writeNexusTree(wTrees.coalEvents, wTrees.nodeTimes, wTrees.names,
                "/Users/jayna/Documents/MATLAB/Duke_Reassortment/2013-07-19/varying_s/psi_10.0_antigenicMu_0.005/s_0.05/tree_200N_7000000_antigenicMu_0.0050_s_0.05_psi_10.0_psij_1.0_samplePeriod_39.0_to_40.0_simNo_0_segment_2.tre");
       // args[4]);

    }
}
