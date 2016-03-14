import Matrix.methods;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;




/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 9/18/12
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimulateTree {

    int n_lineages = 100;

    int[][] b;
    double[] d;
    List<Double> sampleTimes = new ArrayList<Double>();
    List<Integer> sampledLineages = new ArrayList<Integer>();
    List<String> names;
    methods tools = new methods();

    double sampleStartTime = 20*365.25; // in days
    double sampleEndTime = 21*365.25; // in days
    double daysBefore = 50; // in days

    boolean sampleBothPatches = false;
    boolean sampleSeasonalPatch = false;
    int samplingSchemeForTwoPatch = 0; //make these schemes interpretable or add notes
    int samplingSchemeForOnePatchModel = 1; // sample randomly from start time to end time;

    //int samplingScheme = 6;

    class coalescentEvent{

        double timeToCoalescence;
        int coalParent;
        int[] coalDaughters;

        coalescentEvent(double timeCoalescence, Integer coalParent, int[] coalDaughters) {

            this.timeToCoalescence = timeCoalescence;
            this.coalParent = coalParent;
            this.coalDaughters = coalDaughters;
        }

        private double getTimeToCoalescence(){

            return this.timeToCoalescence;
        }
        private Integer getCoalParent() {

            return this.coalParent;
        }
        private int[] getCoalDaughters() {

            return this.coalDaughters;
        }
    }

    public void getTransmissionTree(EpiModel model, int seg) {

        infectionHistory iMatrix = model.getIMatrix();

        System.out.println("iMat size "+iMatrix.birth.size());

        coinfectionHistory icoMatrix = model.getIcoMatrix();
        EpiParams params = new EpiParams();

        // setting the sampling scheme
        int samplingScheme = 0;
        if(sampleBothPatches) {
            samplingScheme = 3;
        }
        else if(sampleSeasonalPatch) {
            samplingScheme = 5;
        }
        else{
            samplingScheme = 4;
        }

        System.out.println("scheme "+samplingScheme);

        int n_seqs = n_lineages;

        if(sampledLineages.isEmpty()) {
            sampleLineages(iMatrix, samplingScheme);
        }

        //create names and label reassortants, patch numbers
        names = createNames(sampledLineages, sampleTimes);
        tools.labelReassortants(names, sampledLineages, iMatrix.reassortant, iMatrix.parentCo, iMatrix.patch, iMatrix.fitness);

        List<List<Integer>> parentLineages = null;

        if(seg==0) {

            parentLineages = getParentLineages(n_seqs, sampledLineages, iMatrix.seg1parent);
        }
        else {
            parentLineages = getParentLineages(n_seqs, sampledLineages, iMatrix.seg2parent);
        }

        b = new int[(n_lineages-1)][2];
        d = new double[2*n_lineages-1];

        List<Integer> currIndividuals = new ArrayList<Integer>(sampledLineages);
        List<Integer> completeIndividuals = new ArrayList<Integer>(sampledLineages);
        List<Double> completeSeqTimes = new ArrayList<Double>(sampleTimes);

        Map<Integer, String> nodeNames = new HashMap<Integer, String>();
        Map<Integer, Double> nodeHeights = new HashMap<Integer, Double>();
        fillNodeMaps(nodeNames, nodeHeights);

        boolean condition = true;

        int loc_b = 0;
        int x_event = 1;

        List<Integer> coalescedList = new ArrayList<Integer>();
        //List<Integer> coalescedList2 = new ArrayList<Integer>();
        while(condition) {

            System.out.println(".coal event "+x_event);

            coalescentEvent coalescentEvent = findMostRecentCoalescence(currIndividuals, parentLineages, iMatrix, seg);

            double timeOfCoalescence = coalescentEvent.getTimeToCoalescence();
            Integer coal_parent = coalescentEvent.getCoalParent();
            int[] coalDaughters = coalescentEvent.getCoalDaughters();

            if(timeOfCoalescence == Double.NEGATIVE_INFINITY) {
                break;
            }

            List<Integer> indiv1_index = tools.find(completeIndividuals, coalDaughters[0]);
            List<Integer> indiv2_index = tools.find(completeIndividuals, coalDaughters[1]);

            List<Integer> non_coalescedList1 = new ArrayList<Integer>();
            List<Integer> non_coalescedList2 = new ArrayList<Integer>();

            for(Integer i: indiv1_index) {

                Integer x = i + 1;
                if(!coalescedList.contains(x)) {

                    non_coalescedList1.add(i);
                }
            }

            for(Integer i: indiv2_index) {

                Integer x = i + 1;
                if(!coalescedList.contains(x)) {

                    non_coalescedList2.add(i);
                }
            }


            int index1 = non_coalescedList1.get(0);
            int index2 = non_coalescedList2.get(0);

            if(index1==index2) {

                index2 = non_coalescedList1.get(1);

            }


            d[index1] = Math.abs(completeSeqTimes.get(index1) - timeOfCoalescence)/365.25;
            d[index2] = Math.abs(completeSeqTimes.get(index2) - timeOfCoalescence)/365.25;


            int[] children = new int[2];

            b[loc_b] = new int[2];
            children[0] = index1+1;
            children[1] = index2+1;
            if(coalescedList.contains(children[0])|| coalescedList.contains(children[1])) {

                System.out.println("coalesced lineage already");


            }
            coalescedList.add(children[0]);
            coalescedList.add(children[1]);

            Arrays.sort(children);
            b[loc_b] = children;

            String nodeName = "";
            String reassState = "";
            if(iMatrix.getReassortant(coal_parent) == 0) {
                reassState = "nr";
            }
            else{
                reassState = "r";
            }

            if(iMatrix.patch.size()>0)  {
                nodeName = "'coal_"+x_event+"_patch_"+iMatrix.getPatch(coal_parent)+"_reass_"+reassState+"_fitness_"+iMatrix.getFitness(coal_parent)+"_node_"+(loc_b+1)+"_"+(timeOfCoalescence/365.25)+"'";
            }
            else{
                nodeName = "'coal_"+x_event+"_patch_0"+"_reass_"+reassState+"_fitness_"+iMatrix.getFitness(coal_parent)+"_node_"+(loc_b+1)+"_"+(timeOfCoalescence/365.25)+"'";

            }

            names.add(nodeName);

            loc_b += 1;

            //check if coalParent is in the complete individuals list - because it can only be there at max - right?
            completeIndividuals.add(coal_parent);


            completeSeqTimes.add(timeOfCoalescence);

            int loc1_inCurr = tools.find(currIndividuals, (Integer) coalDaughters[0]);
            currIndividuals.remove(loc1_inCurr);
            //currIndividuals.remove((Integer) coalDaughters[0]);
            int loc2_inCurr = tools.find(currIndividuals, (Integer) coalDaughters[1]);
            currIndividuals.remove(loc2_inCurr);
            //currIndividuals.remove((Integer) coalDaughters[1]);

            currIndividuals.add(coal_parent);
            n_seqs = currIndividuals.size();

            if(seg==0) {
                parentLineages = getParentLineages(n_seqs, currIndividuals, iMatrix.seg1parent);
            }
            else{
                parentLineages = getParentLineages(n_seqs, currIndividuals, iMatrix.seg2parent);
            }

            x_event++;
        }
        //Not all lineages coalesce - so need to figure out roughly lambda

        double lamba_approx = (double)n_lineages*(((double)n_lineages)-1.0)*(EpiParams.nu_s + EpiParams.nu_co + EpiParams.mu + EpiParams.mu)/(EpiParams.Is_init + EpiParams.Ico_init);
        double k_approx = lamba_approx/((n_seqs+1)*n_seqs);

        k_approx = k_approx*10;

        //if basal b and basal d are at the right size just append on the b and d...don't need to go through the
        //loop below!!
//
// ****commented out the code that fixes the basal part of the infection tree***
//
//        System.out.println("*******CHECKS*******");
//        System.out.println("names a: "+names.size());
//        int limit = loc_b;
//        int basal_limit = n_lineages-limit-1;
//
//        System.out.println("Basal_names limit: "+(n_lineages-limit-1));
//        System.out.println("loc_b "+limit);
//        System.out.println(">1 "+basal_b.size()+" : "+basal_d.size()+" : "+basal_names.size());
//
//        if(basal_b.size() < basal_limit && basal_d.size() < 2*basal_limit && basal_names.size() < basal_limit) {

        //basal_b = new ArrayList<int[]>();
        //basal_d = new HashMap<Integer, Double>();
        //basal_names = new ArrayList<String>();

        while(condition) {
            if(n_seqs==1) {
                break;
            }

            List<Integer> indices = tools.getIndices(currIndividuals);
            Collections.shuffle(indices);

            // these should be unique
            int index_1 = indices.get(0);
            int index_2 = indices.get(1);

            int[] coalDaughters = new int[2];
            //coalDaughters might be the same
            coalDaughters[0] = currIndividuals.get(index_1);
            coalDaughters[1] = currIndividuals.get(index_2);

            List<Integer> indiv1_index = tools.find(completeIndividuals, coalDaughters[0]);
            List<Integer> indiv2_index = tools.find(completeIndividuals, coalDaughters[1]);

            List<Integer> non_coalescedList1 = new ArrayList<Integer>();
            List<Integer> non_coalescedList2 = new ArrayList<Integer>();

            for(Integer i: indiv1_index) {

                Integer x = i + 1;
                if(!coalescedList.contains(x)) {

                    non_coalescedList1.add(i);
                }
            }

            for(Integer i: indiv2_index) {

                Integer x = i + 1;
                if(!coalescedList.contains(x)) {

                    non_coalescedList2.add(i);
                }
            }

            int index1 = non_coalescedList1.get(0);
            int index2 = non_coalescedList2.get(0);

            //this means non_coalesced lists 1 and 2 should be the same, and should contain at most
            if(index1==index2) {

                System.out.println("nc list1: "+non_coalescedList1);
                System.out.println("nc list2: "+non_coalescedList2);

                index2 = non_coalescedList1.get(1);//tools.lastIndex(non_coalescedList1));


            }

            int[] children = new int[2];
            b[loc_b] = new int[2];
            children[0] = index1+1;
            children[1] = index2+1;

            if(children[0]==children[1]) {

                System.out.println("Stop");
            }
            Arrays.sort(children);

            b[loc_b] = children;

            coalescedList.add(children[0]);
            coalescedList.add(children[1]);

            //if(basal_b.size() < basal_limit){
            //basal_b.add(children);
            //}

            double timeOfCoalescence = EpiParams.simulationStartTime -50.0;

            d[index1] = (completeSeqTimes.get(index1) - timeOfCoalescence)/365.25;
            d[index2] = (completeSeqTimes.get(index2) - timeOfCoalescence)/365.25;

            //if(basal_d.size() < (2*basal_limit)){


            //basal_d.put(index1,d[index1]);
            //basal_d.put(index2,d[index2]);

            //}
            int coal_parent = Collections.max(completeIndividuals)+1;
            String nodeName = "'coal_"+x_event+"_patch_0_reass_nr_fitness_"+iMatrix.getFitness(coal_parent)+"_node_"+(loc_b+1)+"_"+(timeOfCoalescence/365.25)+"'";

            x_event++;
            names.add(nodeName);



            loc_b += 1;



            completeIndividuals.add(coal_parent);
            completeSeqTimes.add(timeOfCoalescence);
            currIndividuals.add(coal_parent);

            int loc1_inCurr = tools.find(currIndividuals, (Integer) coalDaughters[0]);
            currIndividuals.remove(loc1_inCurr);

            int loc2_inCurr = tools.find(currIndividuals, (Integer) coalDaughters[1]);
            currIndividuals.remove(loc2_inCurr);

            n_seqs = currIndividuals.size();

        }
        System.out.println("names b: "+names.size());

        //}

        d[tools.lastIndex(d)] = 0.0;


        Collections.sort(coalescedList);
        //Collections.sort(coalescedList2);

        System.out.println(coalescedList.size()) ;
        //System.out.println(coalescedList2.size());


    }

    public void sampleLineages(infectionHistory iMatrix, int samplingScheme) {


        List<Integer> locsNotPrevUsed = tools.getIndices(iMatrix.birth);

//        double simStartTime = params.simulationStartTime;
//        double simEndTime = params.simulationEndTime;
//        double totalSimTime = simEndTime-simStartTime;

        double startTime = this.sampleStartTime;
        double endTime = this.sampleEndTime;
        double totalSampleTime = startTime-endTime;


        double sampleTime;
        Random random = new Random();
        int sampleIndex;
        Integer sampledLineage;

        Sample patchSamples;
        switch(samplingScheme) {

            case 0: //sample uniformly across the epidemic - only use this when sampling from day 0!!!

                int samplingFrequency = (int)Math.floor(totalSampleTime/n_lineages);

                sampleTime = Math.floor(Math.random());

                for(int i=0; i<n_lineages; i++) {

                    List<Integer> locsBornOK = tools.find(iMatrix.birth, 5, sampleTime);

                    List<Integer> locsDeathOK = tools.find(iMatrix.birth, 1, sampleTime);
                    List<Integer> locsOK = tools.intersect(tools.intersect(locsBornOK, locsDeathOK), locsNotPrevUsed);

                    double[] deaths_locOK = tools.copyArrayRangeOf(iMatrix.death, locsOK);
                    Arrays.sort(deaths_locOK);

                    double min_death = deaths_locOK[0];

                    List<Integer> loc_indiv = tools.find(deaths_locOK, 0, min_death);
                    Collections.shuffle(loc_indiv);

                    sampledLineage = locsOK.get(loc_indiv.get(0));

                    sampledLineages.add(sampledLineage);
                    sampleTimes.add(iMatrix.getDeath(sampledLineage));
                    locsNotPrevUsed.remove(sampledLineage);
                    sampleTime += samplingFrequency;

                }
                break;
            case 1: //sampling randomly across the epidemic - only use this method when sampling from day 0!!!

                for(int i=0; i<n_lineages; i++) {

                    System.out.println("> sample "+(i+1));
                    sampleTime = (Math.random()*totalSampleTime);
                    List<Integer> alive = tools.findAliveLocs(iMatrix.birth, iMatrix.death, sampleTime);
                    List<Integer> locsOK = tools.intersect(alive, locsNotPrevUsed);

                    if(locsOK.isEmpty()) {

                        while(locsOK.isEmpty()) {
                            sampleTime = (Math.random()*totalSampleTime);
                            alive = tools.findAliveLocs(iMatrix.birth, iMatrix.death, sampleTime);
                            locsOK = tools.intersect(alive, locsNotPrevUsed);
                        }
                    }
                    //what does this step do?
                    double[] deaths_locOK = tools.copyArrayRangeOf(iMatrix.death, locsOK);
                    Arrays.sort(deaths_locOK);

                    //sample with min_death
                    double min_death = deaths_locOK[0];

                    List<Integer> loc_indiv = tools.find(deaths_locOK, min_death);//find(deaths_locOK, 0, min_death);
                    Collections.shuffle(loc_indiv, random);

                    sampledLineage = locsOK.get(loc_indiv.get(0));

                    sampledLineages.add(sampledLineage);
                    sampleTimes.add(iMatrix.getDeath(sampledLineage)); //sampleTime is the recovery or death time;
                    locsNotPrevUsed.remove(sampledLineage); //= setDiff(locsNotPrevUsed, tempSampledLineages);


                }

                break;
            case 2:

                for(int i=0; i<n_lineages; i++) {


                    //sample only from the last 50 days of the simulation
                    sampleTime = (Math.random()*(totalSampleTime-daysBefore))+daysBefore;

                    List<Integer> locsBornOK = tools.find(iMatrix.birth, -1, sampleTime);
                    List<Integer> bornPrev = tools.find(iMatrix.birth, 2, Double.NaN);
                    locsBornOK = tools.union(locsBornOK, bornPrev);

                    List<Integer> locsDeathOK = tools.find(iMatrix.death, 1, sampleTime);
                    List<Integer> alive = tools.intersect(locsBornOK, locsDeathOK);
                    List<Integer> locsOK = tools.intersect(alive, locsNotPrevUsed);


                    double[] deaths_locOK = tools.copyArrayRangeOf(iMatrix.death, locsOK);
                    Arrays.sort(deaths_locOK);


                    double min_death = deaths_locOK[0];

                    List<Integer> loc_indiv = tools.find(deaths_locOK, 0, min_death);
                    Collections.shuffle(loc_indiv);

                    sampledLineage = locsOK.get(loc_indiv.get(0));//locsOK.get(sampleIndex);

                    sampledLineages.add(sampledLineage);
                    sampleTimes.add(iMatrix.getDeath(sampledLineage)); //sampleTime is the recovery or death time;
                    locsNotPrevUsed.remove(sampledLineage);
                }
                break;
            case 3:   //sample both patches

                patchSamples = getTwoPatchSamples(samplingSchemeForTwoPatch, iMatrix);

                sampledLineages.addAll(patchSamples.getLineages());
                sampleTimes.addAll(patchSamples.getTimes());

                break;

            case 4: // sample for onePatch model

                //samples randomly from sim startTime and sample endTime; from 0 to END;
                patchSamples = getSamples(n_lineages, samplingSchemeForOnePatchModel, iMatrix);

                sampledLineages.addAll(patchSamples.getLineages());
                sampleTimes.addAll(patchSamples.getTimes());

                break;

            case 5: // sample only seasonal patch in the twoPatch model

                patchSamples = getPatchSamples(2, n_lineages, samplingSchemeForOnePatchModel, iMatrix);

                sampledLineages.addAll(patchSamples.getLineages());
                sampleTimes.addAll(patchSamples.getTimes());



        }
        sampledLineages = tools.sortListsByTime(sampleTimes, sampledLineages);

        Set<Integer> set = new HashSet<Integer>();
        set.addAll(sampledLineages);
        System.out.println(sampledLineages.size() + ":" + set.size());


    }

    public List<List<Integer>> getParentLineages(int n_seq, List<Integer> sampledLineages, List<Integer> parents) {

        List<List<Integer>> parentLineages = new ArrayList<List<Integer>>();

        for(int i=0; i<n_seq; i++) {

            List<Integer> lineages = getSampleLineage(sampledLineages.get(i), parents);
            parentLineages.add(lineages);

        }

        return parentLineages;
    }

    public List<Integer> getSampleLineage(Integer indiv_lineage, List<Integer> parents) {

        int i = 0;

        List<Integer> this_lineage = new ArrayList<Integer>();

        this_lineage.add(indiv_lineage);

        boolean condition = true;
        while(condition) {

            Integer lineage_parent = parents.get(this_lineage.get(i));
            if(lineage_parent < 0) {

                return this_lineage;

            }
            else{

                this_lineage.add(lineage_parent);
                i++;
            }
        }

        return this_lineage;
    }

    //returns timeOfCoalescence

    public coalescentEvent findMostRecentCoalescence(List<Integer> current_indiv, List<List<Integer>> parentLineages, infectionHistory iMatrix, coinfectionHistory icoMatrix) {


        double mostRecentTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double coalTime1 = -1.0;
        double coalTime2 = -1.0;
        int n_individuals_sampled = current_indiv.size();
        int[] coalDaughters = new int[2];
        Integer coalParent = (int)Double.NEGATIVE_INFINITY;

        infectionHistory local_iMatrix = new infectionHistory();


        // System.out.println("> "+current_indiv.size());
        coalescentEvent coalescentEvent = null;
        List<Integer> currIndividuals = new ArrayList<Integer>(current_indiv);

        for(int i=0; i<n_individuals_sampled;i++){

            for(int j=(i+1); j<n_individuals_sampled; j++) {

                //System.out.println(i+ " :" + j);
                Integer oldestParent_i = parentLineages.get(i).get(tools.lastIndex(parentLineages.get(i)));
                Integer oldestParent_j = parentLineages.get(j).get(tools.lastIndex(parentLineages.get(j)));

                if(iMatrix.getParentCo(oldestParent_i) >=0 && iMatrix.getParent(oldestParent_i) < 0) {  // if there is a parent_co but no parent_s then this is a coinfected transmission

                    while(iMatrix.getParentCo(oldestParent_i) >=0 && iMatrix.getParent(oldestParent_i) < 0)  {  // there could multiple coinfected transmission from multiple reassortment events
                        //find Ico_matrix entry
                        //int[] coinfectedHistory = Ico_matrix.get(parent_co.get(oldestParent_i]);
                        int coinfected = iMatrix.getParentCo(oldestParent_i);
                        //choose parent randomly to be the genetic donor of the lineage i;

                        int parentIndex = (int)Math.round(Math.random());   //the parentIndex could be either 1 or 2
                        Integer donorParent = (int)Double.NEGATIVE_INFINITY;//coinfectedHistory[parentIndex];
                        if(parentIndex==0) {

                            donorParent = icoMatrix.getParent1(coinfected);
                        }
                        else{

                            donorParent = icoMatrix.getParent2(coinfected);
                        }

                        if(donorParent>=0) {   // could be -1 if it is a reassortant from the initial coinfecteds - these can't be detected!!
                            List<Integer> donorLineage = getSampleLineage(donorParent,iMatrix.parent);
                            iMatrix.setParent(oldestParent_i, donorParent);
                            //System.out.println(oldestParent_i+ " : "+donorParent);
                            parentLineages.get(i).addAll(donorLineage);

                            oldestParent_i = parentLineages.get(i).get(tools.lastIndex(parentLineages.get(i)));
                        }
                        else{
                            break;
                        }
                    }


                }
                else if(iMatrix.getParentCo(oldestParent_j) >=0 && iMatrix.getParent(oldestParent_j) < 0) {

                    // while it is a primary reassortant....
                    while(iMatrix.getParentCo(oldestParent_j)>=0 && iMatrix.getParent(oldestParent_j) < 0)  {

                        //find Ico_matrix entry
                        //int[] coinfectedHistory = Ico_matrix.get(parent_co[oldestParent_j]);
                        int coinfected = iMatrix.getParentCo(oldestParent_j);
                        //choose parent randomly to be the genetic donor of the lineage i;
                        int parentIndex = (int)Math.round(Math.random());//+2;
                        Integer donorParent = (int)Double.NEGATIVE_INFINITY;//coinfectedHistory[parentIndex];
                        if(parentIndex==0) {

                            donorParent = icoMatrix.getParent1(coinfected);
                        }
                        else{

                            donorParent = icoMatrix.getParent2(coinfected);
                        }

                        if(donorParent>=0) {

                            List<Integer> donorLineage = getSampleLineage(donorParent,iMatrix.parent);

                            iMatrix.setParent(oldestParent_j, donorParent);
                            //System.out.println("j pl "+parentLineages.get(j));
                            //System.out.println(oldestParent_j+ " : "+donorParent);

                            parentLineages.get(j).addAll(donorLineage);
                            //System.out.println("j pl "+parentLineages.get(j));

                            oldestParent_j = parentLineages.get(j).get(tools.lastIndex(parentLineages.get(j)));
                        }
                        else{
                            break;
                        }
                    }
                }

                //could make this faster?
                List<Integer> parentsInCommon = tools.intersect(parentLineages.get(i), parentLineages.get(j));
                //tools.intersect(parentLineages.get(i), parentLineages.get(j));

                if(!parentsInCommon.isEmpty()) {

                    Integer this_parent = Collections.max(parentsInCommon);

                    int loc1 = tools.find(parentLineages.get(i), this_parent);
                    int loc2 = tools.find(parentLineages.get(j), this_parent);

                    if(currIndividuals.get(i).equals(this_parent)) {

                        coalTime1 = iMatrix.getDeath(currIndividuals.get(i));   //this could be a problem is death is -1!!
                    }
                    else{
                        coalTime1 = iMatrix.getBirth(parentLineages.get(i).get(loc1-1)); //birth of curr individual i's parent
                    }


                    if(currIndividuals.get(j).equals(this_parent)) {

                        coalTime2 = iMatrix.getDeath(currIndividuals.get(j));

                    }
                    else{

                        coalTime2 = iMatrix.getBirth(parentLineages.get(j).get(loc2-1));
                    }

                    thisTimeOfCoalescence = Math.min(coalTime1, coalTime2);//tools.min(coalTime1, coalTime2);

                    if(thisTimeOfCoalescence > mostRecentTimeOfCoalescence) {

                        if(currIndividuals.get(i).equals(currIndividuals.get(j))) {

                            thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;
                        }
                        coalDaughters = new int[2];
                        coalDaughters[0] = currIndividuals.get(i);
                        coalDaughters[1] = currIndividuals.get(j);
                        coalParent = this_parent;
                        mostRecentTimeOfCoalescence = thisTimeOfCoalescence;

                    }

                }

            }

        }

        coalescentEvent = new coalescentEvent(mostRecentTimeOfCoalescence, coalParent, coalDaughters);

        return coalescentEvent;

    }

    public coalescentEvent findMostRecentCoalescence(List<Integer> current_indiv, List<List<Integer>> parentLineages, infectionHistory iMatrix, int seg) {


        double mostRecentTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double coalTime1 = -1.0;
        double coalTime2 = -1.0;
        int n_individuals_sampled = current_indiv.size();
        int[] coalDaughters = new int[2];
        Integer coalParent = (int)Double.NEGATIVE_INFINITY;

        infectionHistory local_iMatrix = new infectionHistory();


        // System.out.println("> "+current_indiv.size());
        coalescentEvent coalescentEvent = null;
        List<Integer> currIndividuals = new ArrayList<Integer>(current_indiv);


        for(int i=0; i<n_individuals_sampled;i++){

            Integer oldestParent_i = parentLineages.get(i).get(tools.lastIndex(parentLineages.get(i)));


            for(int j=(i+1); j<n_individuals_sampled; j++) {

                //System.out.println(i+ " :" + j);
                Integer oldestParent_j = parentLineages.get(j).get(tools.lastIndex(parentLineages.get(j)));

                if(iMatrix.getParentCo(oldestParent_i) >=0 && iMatrix.getParent(oldestParent_i) < 0) {  // if there is a parent_co but no parent_s then this is a coinfected transmission

                    while(iMatrix.getParentCo(oldestParent_i) >=0 && iMatrix.getParent(oldestParent_i) < 0)  {  // there could multiple coinfected transmission from multiple reassortment events

                        Integer donorParent = (int)Double.NEGATIVE_INFINITY;//coinfectedHistory[parentIndex];
                        if(seg==0) {

                            donorParent = iMatrix.getSeg1parent(i);
                        }
                        else{

                            donorParent = iMatrix.getSeg2parent(i);
                        }

                        if(donorParent>=0) {   // could be -1 if it is a reassortant from the initial coinfecteds - these can't be detected!!
                            List<Integer> donorLineage = getSampleLineage(donorParent,iMatrix.parent);
                            iMatrix.setParent(oldestParent_i, donorParent);
                            //System.out.println(oldestParent_i+ " : "+donorParent);
                            parentLineages.get(i).addAll(donorLineage);

                            oldestParent_i = parentLineages.get(i).get(tools.lastIndex(parentLineages.get(i)));
                        }
                        else{
                            break;
                        }
                    }


                }
                else if(iMatrix.getParentCo(oldestParent_j) >=0 && iMatrix.getParent(oldestParent_j) < 0) {

                    // while it is a primary reassortant....
                    while(iMatrix.getParentCo(oldestParent_j)>=0 && iMatrix.getParent(oldestParent_j) < 0)  {

                       Integer donorParent = (int)Double.NEGATIVE_INFINITY;//coinfectedHistory[parentIndex];
                        if(seg==0) {

                            donorParent = iMatrix.getSeg1parent(j);
                        }
                        else{

                            donorParent = iMatrix.getSeg2parent(j);
                        }

                        if(donorParent>=0) {

                            List<Integer> donorLineage = getSampleLineage(donorParent,iMatrix.parent);

                            iMatrix.setParent(oldestParent_j, donorParent);

                            parentLineages.get(j).addAll(donorLineage);

                            oldestParent_j = parentLineages.get(j).get(tools.lastIndex(parentLineages.get(j)));
                        }
                        else{
                            break;
                        }
                    }
                }

                //could make this faster?
                List<Integer> parentsInCommon = tools.intersect(parentLineages.get(i), parentLineages.get(j));
                //tools.intersect(parentLineages.get(i), parentLineages.get(j));

                if(!parentsInCommon.isEmpty()) {

                    Integer this_parent = Collections.max(parentsInCommon);

                    int loc1 = tools.find(parentLineages.get(i), this_parent);
                    int loc2 = tools.find(parentLineages.get(j), this_parent);

                    if(currIndividuals.get(i).equals(this_parent)) {

                        coalTime1 = iMatrix.getDeath(currIndividuals.get(i));   //this could be a problem is death is -1!!
                    }
                    else{
                        coalTime1 = iMatrix.getBirth(parentLineages.get(i).get(loc1-1)); //birth of curr individual i's parent
                    }


                    if(currIndividuals.get(j).equals(this_parent)) {

                        coalTime2 = iMatrix.getDeath(currIndividuals.get(j));

                    }
                    else{

                        coalTime2 = iMatrix.getBirth(parentLineages.get(j).get(loc2-1));
                    }

                    thisTimeOfCoalescence = Math.min(coalTime1, coalTime2);//tools.min(coalTime1, coalTime2);

                    if(thisTimeOfCoalescence > mostRecentTimeOfCoalescence) {

                        if(currIndividuals.get(i).equals(currIndividuals.get(j))) {

                            thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;
                        }
                        coalDaughters = new int[2];
                        coalDaughters[0] = currIndividuals.get(i);
                        coalDaughters[1] = currIndividuals.get(j);
                        coalParent = this_parent;
                        mostRecentTimeOfCoalescence = thisTimeOfCoalescence;

                    }

                }

            }

        }


        coalescentEvent = new coalescentEvent(mostRecentTimeOfCoalescence, coalParent, coalDaughters);

        return coalescentEvent;

    }

    public coalescentEvent findMostRecentCoalescence(List<Integer> current_indiv, List<List<Integer>> parentLineages, double[] births, double[] deaths) {


        double mostRecentTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double coalTime1 = -1.0;
        double coalTime2 = -1.0;
        int n_individuals_sampled = current_indiv.size();
        int[] coalDaughters = new int[2];
        Integer coalParent = (int)Double.NEGATIVE_INFINITY;

        coalescentEvent coalescentEvent = null;
        List<Integer> currIndividuals = new ArrayList<Integer>(current_indiv);

        for(int i=0; i<n_individuals_sampled;i++){

            for(int j=(i+1); j<n_individuals_sampled; j++){

                List<Integer> parentsInCommon = tools.intersect(parentLineages.get(i), parentLineages.get(j));

                if(!parentsInCommon.isEmpty()) {

                    Integer this_parent = Collections.max(parentsInCommon);

                    int loc1 = tools.find(parentLineages.get(i), this_parent);
                    int loc2 = tools.find(parentLineages.get(j), this_parent);

                    if(currIndividuals.get(i).equals(this_parent)) {

                        coalTime1 = deaths[currIndividuals.get(i)];
                    }
                    else{

                        coalTime1 = births[parentLineages.get(i).get(loc1-1)];
                    }

                    if(currIndividuals.get(j).equals(this_parent)) {

                        coalTime2 = deaths[currIndividuals.get(j)];
                    }
                    else{

                        coalTime2 = births[parentLineages.get(j).get(loc2-1)];
                    }

                    thisTimeOfCoalescence = Math.min(coalTime1, coalTime2);


                    if(currIndividuals.get(i).equals(currIndividuals.get(j))) {
                        //what do we here?
                        //well it shouldn't if we have sampled the lineages correctly - DUH

                    }
                    if(thisTimeOfCoalescence > mostRecentTimeOfCoalescence) {

                        if(currIndividuals.get(i).equals(currIndividuals.get(j))) {

                            thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;

                        }
                        coalDaughters = new int[2];
                        coalDaughters[0] = currIndividuals.get(i);
                        coalDaughters[1] = currIndividuals.get(j);
                        coalParent = this_parent;
                        mostRecentTimeOfCoalescence = thisTimeOfCoalescence;

                    }

                }

            }

        }

        coalescentEvent = new coalescentEvent(mostRecentTimeOfCoalescence, coalParent, coalDaughters);

        return coalescentEvent;

    }

    private Sample getTwoPatchSamples(int samplingScheme, infectionHistory iMatrix) {

        Sample patch1Samples = new Sample();
        Sample patch2Samples = new Sample();
        switch(samplingScheme) {

            case 0:

                patch1Samples = getPatchSamples(1, n_lineages/2, 0, iMatrix);
                patch2Samples = getPatchSamples(2, n_lineages/2, 1, iMatrix);

                break;

            case 1:

                patch1Samples = getPatchSamples(1, n_lineages/2, 1, iMatrix);
                patch2Samples = getPatchSamples(2, n_lineages/2, 1, iMatrix);

                break;

            case 2:

                patch1Samples = getPatchSamples(1, n_lineages/2, 0, iMatrix);
                patch2Samples = getPatchSamples(2, n_lineages/2, 0, iMatrix);

                break;

            case 3:
                patch1Samples = getPatchSamples(1, n_lineages/2, 1, iMatrix);
                patch2Samples = getPatchSamples(2, n_lineages/2, 0, iMatrix);
                break;
        }

        Sample samples = new Sample();
        samples.addSamples(patch1Samples);
        samples.addSamples(patch2Samples);

        return samples;

    }

    public Sample getSamples(int n_samples, int samplingScheme, infectionHistory iMatrix) {

        Sample samples = new Sample();

        Map<Integer, Double> births = new HashMap<Integer, Double>();
        Map<Integer, Double> deaths = new HashMap<Integer, Double>();

        double startTime = sampleStartTime;
        double endTime = sampleEndTime;
        double totalSampleTime = endTime-startTime;

        Set<Integer> locsNotPrevUsed = new HashSet<Integer>();
        List<Integer> indices = new ArrayList<Integer>();

        System.out.println("imat size "+iMatrix.birth.size());

        for(int i=0; i<iMatrix.birth.size(); i++) {

            if(iMatrix.getBirth(i) >= startTime && iMatrix.getBirth(i) < endTime && iMatrix.getDeath(i)>0) {


//                births.put(i, iMatrix.getBirth(i));
//                deaths.put(i, iMatrix.getDeath(i));
                    indices.add(i);

                //locsNotPrevUsed.add(i);
            }

        }


        System.out.println(indices.size());
        int samplingFrequency = (int)Math.floor(totalSampleTime/(n_samples));

//
        //locsNotPrevUsed.addAll(births.keySet());
        Random random = new Random();

        //indices.addAll(locsNotPrevUsed);
        System.out.println(">>" + locsNotPrevUsed.size());
        Collections.shuffle(indices, random);


        List<Integer> sampledLineages = new ArrayList<Integer>();
        List<Double> sampleTimes = new ArrayList<Double>();
        double sampleTime;
        switch(samplingScheme) {
            // uniform sampling could be better!!!
            case 0:  // samples uniformly from 0 to sampleEndTime;

                double intervalEnd = startTime+samplingFrequency;
                double intervalStart = startTime;


                for(int i=0; i<(n_samples); i++) {

                    double intervalPeriod = intervalEnd-samplingFrequency;
                    sampleTime = random.nextDouble()*intervalPeriod + intervalStart;

                    //System.out.println("> sample "+((n_samples/2)+i));
                    Set<Integer> alive = tools.findAliveLocs(births, deaths, sampleTime);
                    List<Integer> locsOK = tools.intersect(alive, locsNotPrevUsed);


                    if(locsOK.isEmpty()) {

                        while(locsOK.isEmpty()) {
                            sampleTime = random.nextDouble()*intervalPeriod + intervalStart;
                            alive = tools.findAliveLocs(births, deaths, sampleTime);
                            locsOK = tools.intersect(alive, locsNotPrevUsed);
                        }
                    }
                    //what does this step do?
                    double[] deaths_locOK = tools.copyArrayRangeOf(iMatrix.death, locsOK);
                    Arrays.sort(deaths_locOK);

                    //sample with min_death
                    double min_death = deaths_locOK[0];

                    List<Integer> loc_indiv = tools.find(deaths_locOK, min_death);//find(deaths_locOK, 0, min_death);
                    Collections.shuffle(loc_indiv, random);

                    Integer sampledLineage = locsOK.get(loc_indiv.get(0));

                    sampledLineages.add(sampledLineage);
                    sampleTimes.add(iMatrix.death.get(sampledLineage)); //sampleTime is the recovery or death time;
                    locsNotPrevUsed.remove(sampledLineage); //= setDiff(locsNotPrevUsed, tempSampledLineages);


                    intervalEnd += samplingFrequency;
                    intervalStart += samplingFrequency;


                }

                break;
            case 1:   // samples randomly from startTime to EndTime;

                System.out.println("indices "+indices.size());
                for(int i=0; i<(n_samples); i++) {

                    Integer sampledLineage = indices.get(i);
                    Double sampledTime = iMatrix.death.get(sampledLineage);

                    sampledLineages.add(sampledLineage);
                    sampleTimes.add(sampledTime);

                }
                break;

        }

        samples.addSamples(sampledLineages, sampleTimes);
        return samples;
    }

    public Sample getPatchSamples(int patchNo, int n_patchSamples, int samplingScheme, infectionHistory iMatrix) {

        Sample patchSamples = new Sample();

        double startTime = sampleStartTime;
        double endTime = sampleEndTime;
        double totalSampleTime = endTime-startTime;
        List<Integer> patchIndices = new ArrayList<Integer>();
        int samplingFrequency = (int)Math.floor(totalSampleTime/(n_patchSamples));
        for(int i=0; i<iMatrix.patch.size(); i++) {

            if(iMatrix.getPatch(i)==patchNo) {

                if(iMatrix.getBirth(i) >= startTime && iMatrix.getBirth(i) < endTime && iMatrix.getDeath(i)>startTime) {
                    patchIndices.add(i);
                }
            }

        }
        Random random = new Random();
        Collections.shuffle(patchIndices, random);

        //List<Integer> locsNotPrevUsed = new ArrayList<Integer>();
 //       Set<Integer> locsNotPrevUsed = new HashSet<Integer>();
//
//        locsNotPrevUsed = patchBirths.keySet();    // this is only for the old method (from Katia's matlab code). My method only uses patchindices

        List<Integer> patchSampledLineages = new ArrayList<Integer>();
        List<Double> patchSampleTimes = new ArrayList<Double>();

        switch(samplingScheme) {

            //I think case 0 and 1 is pointless - get rid of it...
//            case 0:  // samples uniformly from 0 to sampleEndTime;
//
//                //this doesn't actually uniformly sample between 0 and sample end time
//                //could samply uniformly along the locs as they're sorted by time of birth
//
//                double intervalEnd = startTime+samplingFrequency;
//                double intervalStart = startTime;
//
//                for(int i=0; i<(n_patchSamples); i++) {
//
//                    double intervalPeriod = intervalEnd-samplingFrequency;
//                    sampleTime = random.nextDouble()*intervalPeriod + intervalStart;
//
//                    //System.out.println("> sample "+((n_samples/2)+i));
//                    Set<Integer> alive = tools.findAliveLocs(patchBirths, patchDeaths, sampleTime);
//                    List<Integer> locsOK = tools.intersect(alive, locsNotPrevUsed);
//
//
//                    if(locsOK.isEmpty()) {
//
//                        while(locsOK.isEmpty()) {
//                            sampleTime = random.nextDouble()*intervalPeriod + intervalStart;
//                            alive = tools.findAliveLocs(patchBirths, patchDeaths, sampleTime);
//                            locsOK = tools.intersect(alive, locsNotPrevUsed);
//                        }
//                    }
//                    //what does this step do?
//                    double[] deaths_locOK = tools.copyArrayRangeOf(iMatrix.death, locsOK);
//                    Arrays.sort(deaths_locOK);
//
//                    //sample with min_death
//                    double min_death = deaths_locOK[0];
//
//                    List<Integer> loc_indiv = tools.find(deaths_locOK, min_death);//find(deaths_locOK, 0, min_death);
//                    Collections.shuffle(loc_indiv, random);
//
//                    Integer sampledLineage = locsOK.get(loc_indiv.get(0));
//
//                    patchSampledLineages.add(sampledLineage);
//                    patchSampleTimes.add(patchDeaths.get(sampledLineage)); //sampleTime is the recovery or death time;
//                    locsNotPrevUsed.remove(sampledLineage); //= setDiff(locsNotPrevUsed, tempSampledLineages);
//
//                    intervalEnd += samplingFrequency;
//                    intervalStart += samplingFrequency;
//
//                }
//                break;
            case 1:   // samples randomly from 0 to sampleEndTime;


                for(int i=0; i<(n_patchSamples); i++) {

//                    int index = (int)(Math.floor(Math.random()*patchIndices.size()));
                    Integer sampledLineage = patchIndices.get(i);

                    Double sampledTime = iMatrix.getDeath().get(sampledLineage);
                    patchSampledLineages.add(sampledLineage);
                    patchSampleTimes.add(sampledTime);

                    //patchIndices.remove(index);

                }
                break;


        }

        patchSamples.addSamples(patchSampledLineages, patchSampleTimes);
        return patchSamples;

    }

    private List<String> createNames(List<Integer> sampledLineages, List<Double> sampleTimes) {

        List<String> names = new ArrayList<String>();

        for(int i=0; i<sampledLineages.size(); i++) {

            String name = "sample_"+(i+1)+"_"+Double.toString((sampleTimes.get(i)) / 365.25);
            names.add(name);
        }

        return names;
    }

    public void appendList(int[][] list, List<int[]> listToAppend, int n_lineages) {

        int loc_b = (n_lineages-listToAppend.size()-1);
        System.out.println("loc_b "+loc_b);
        for(int i=0; i<listToAppend.size(); i++) {
            int[] children = listToAppend.get(i);
            list[loc_b+i] = children;
        }
    }

    public void appendList(List<String> list, List<String> listToAppend, int n_lineages) {

        int loc_b = (n_lineages-listToAppend.size()-1);
        System.out.println("loc_b "+loc_b);
        System.out.println("b_names "+listToAppend.size());
        System.out.println("l_names "+list.size());
        for(int i = 0; i< listToAppend.size(); i++) {

            //while((loc_b+i)<=(2*n_lineages-1)) {

            if((n_lineages+loc_b) < list.size()) {
                //if(list.get((loc_b+i))!= null) {

                list.set((n_lineages+loc_b), listToAppend.get(i));
                loc_b++;

            }
            else{
                list.add((n_lineages+loc_b), listToAppend.get(i));
                loc_b++;
            }
            //}
        }

    }

    public void appendList(double[] list, Map<Integer, Double> mapToAppend) {

        for(Integer i: mapToAppend.keySet()) {

            list[i] = mapToAppend.get(i);
        }
    }

    private void fillNodeMaps(Map<Integer, String> nodeNames, Map<Integer, Double> nodeHeights){
        //, List<Integer> sampledLineages, List<Double> sampledTimes) {

        for(int i=0; i<sampledLineages.size(); i++) {
            nodeNames.put(i+1, "'"+names.get(i)+"'");
            nodeHeights.put(i+1, sampleTimes.get(i));

        }
    }

    public void writeNexusTree(int[][] coalescentEvents, double[] nodeTimes, List<String> names, File treeFile) {

        Map<Integer, String> nodeNames = new HashMap<Integer, String>();
        Map<Integer, Double> nodeHeights = new HashMap<Integer, Double>();
        fillNodeMaps(nodeNames, nodeHeights);

        FileWriter writer4 = null;
        //File treeOutputFile = new File("tree_"+n_lineages+"N_"+params.N+"_antigenicMu_"+params.antigenMu+"_s_"+params.dfe+"_psi_"+params.psi+"_psij_"+params.psi_j+"_samplePeriod_"+sampleStartTime/365.25+"_to_"+sampleEndTime/365.25+"_simNo_"+simNo+"_segment_"+(i+1)+".tre");

        try {
            writer4 = new FileWriter(treeFile);


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
                //System.out.println(partsNode1.length);
                //System.out.println(partsNode1.length);
                String patchInfo = "?";
                String reassInfo = "?";
                String fitnessInfo = "?";

                String patchNode1 = "?";
                String patchNode2 = "?";
                String reassNode1 = "?";
                String reassNode2 = "?";
                String fitness1 = "?";
                String fitness2 = "?";


                if(partsNode1.length > 2 && partsNode2.length > 2) {
                    patchNode1 = partsNode1[1];
                    patchNode2 = partsNode2[1];
                    reassNode1 = partsNode1[2];
                    reassNode2 = partsNode2[2];
                    fitness1 = partsNode1[5];
                    fitness2 = partsNode2[5];

                }


                if(partsParent.length > 5)  {
                    patchInfo = partsParent[3];
                    reassInfo = partsParent[5];
                    fitnessInfo = partsParent[7];

                }


                treeString = "("+node1+"[&patch=\""+patchNode1+"\",reassortant=\""+reassNode1+"\",fitness="+fitness1+"]:"+nodeTimes[child1-1]+
                        ","+node2+"[&patch=\""+patchNode2+"\",reassortant=\""+reassNode2+"\",fitness="+fitness2+"]:"+nodeTimes[child2-1]
                        +")[&parent="+names.get(n_lineages+k)+",parentReass=\""+reassInfo+"\""+",parentPatch=\""+patchInfo+"\"]";

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

    public void getTransmissionTrees(int n_segments, EpiModel model, int simNo) {

        //EpiModel model = new reassortmentTwoPatch_new();
        EpiParams params = new EpiParams();
        //model.runSimulation(params);

        FileWriter writer1 = null;
        FileWriter writer2 = null;
        FileWriter writer3 = null;
        //FileWriter writer4 = null;

        EpiModel.iMatrix.segment1Fitness.clear();
        EpiModel.iMatrix.fitness.clear();

        for(int i=0; i<n_segments; i++) {

            //File outputfile_b = new File("coalescentPairs_epidemic_"+n_lineages+"N_"+params.N+"_segment_mij_"+params.m_ij+"_psii_"+params.psi_i+"_psij_"+params.psi_j+"_"+(i)+".csv");
            File outputfile_b = new File("coalescentPairs_endemic_"+n_lineages+"_N_"+ EpiParams.N +"_antigenicMu_"+ EpiParams.antigenicMu_a +"_s_"+ EpiParams.s_b1 +"_psi_"+ EpiParams.psi_i +"_psij_"+ EpiParams.psi_j +"_samplePeriod_"+sampleStartTime/365.25+"_to_"+sampleEndTime/365.25+"_simNo_"+simNo+"_segment_"+(i+1)+".csv");

            //File outputfile_d = new File("nodeTimes_epidemic_100_20e6_segment_mij_"+params.m_ij+"_psii_"+params.psi_i+"_psij_"+params.psi_j+"_"+(i)+".csv");
            File outputfile_d = new File("nodeTimes_endemic_"+n_lineages+"N_"+ EpiParams.N +"_antigenicMu_"+ EpiParams.antigenicMu_a +"_s_"+ EpiParams.s_b1 +"_psi_"+ EpiParams.psi_i +"_psij_"+ EpiParams.psi_j +"_samplePeriod_"+sampleStartTime/365.25+"_to_"+sampleEndTime/365.25+"_simNo_"+simNo+"_segment_"+(i+1)+".csv");

            //File outputfile_n = new File("names_epidemic_100_20e6_segment_mij_"+params.m_ij+"_psii_"+params.psi_i+"_psij_"+params.psi_j+"_"+(i)+".csv");
            File outputfile_n = new File("names_endemic_"+n_lineages+"N_"+ EpiParams.N +"_antigenicMu_"+ EpiParams.antigenicMu_a +"_s_"+ EpiParams.s_b1 +"_psi_"+ EpiParams.psi +"_psij_"+ EpiParams.psi_j +"_samplePeriod_"+sampleStartTime/365.25+"_to_"+sampleEndTime/365.25+"_simNo_"+simNo+"_segment_"+(i+1)+".csv");

            File treeOutputFile = new File("tree_"+n_lineages+"N_"+ EpiParams.N +"_p_"+ EpiParams.p +"_antigenicMu_"+ EpiParams.antigenicMu_a +"_s_"+ EpiParams.s_b1 +"_psi_"+ EpiParams.psi +"_samplePeriod_"+(Math.ceil(sampleStartTime/365.25))+"_to_"+(Math.ceil(sampleEndTime/365.25))+"_simNo_"+simNo+"_segment_"+(i+1)+".tre");

//            try {
//                writer1 = new FileWriter(outputfile_b);
//                writer2 = new FileWriter(outputfile_d);
//                writer3 = new FileWriter(outputfile_n);
//                //writer4 = new FileWriter(treeOutputFile);
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }


            System.out.println("Segment "+ (i+1));


            getTransmissionTree(model, i);


            int[][] coalescentEvents = b;
            double[] nodeTimes = d;

            System.out.println("cl "+coalescentEvents.length);
            System.out.println("names "+names.size());
            System.out.println("nodeTimes "+nodeTimes.length);

//            try {

//
//                int j=0;
//                while (j < coalescentEvents.length) {
//                    if(writer1 != null) {
//                        writer1.write(coalescentEvents[j][0]+","+coalescentEvents[j][1]+"\n");
//                        writer1.flush();
//                    }
//                    j++;
//                }
//                for (double nodeTime : nodeTimes) {
//                    if (writer2 != null) {
//                        writer2.write(nodeTime + "\n");
//                        writer2.flush();
//                    }
//                }
//
//                for (String name : names) {
//
//                    if (writer3 != null) {
//                        writer3.write("\"" + name + "\"\n");
//                        writer3.flush();
//                    }
//                }
//
//                Map<Integer, String> nodeNames = new HashMap<Integer, String>();
//                Map<Integer, Double> nodeHeights = new HashMap<Integer, Double>();
//                fillNodeMaps(nodeNames, nodeHeights);

                //writing trees:

                writeNexusTree(coalescentEvents, nodeTimes, names, treeOutputFile);

//                assert writer1 != null;
//                writer1.close();
//                assert writer2 != null;
//                writer2.close();
//                assert writer3 != null;
//                writer3.close();
                //writer4.close();
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }

        }
    }

    public static void runOnePatchSim(String [] args) {

        EpiModel model = new reassortmentTauLeap();
        EpiParams params = new EpiParams();

        int n_segments = 2;
        int n_seasons = 40;
        int n_lineages = 300;

        n_lineages = Integer.parseInt(args[0]);
        EpiParams.p = Double.parseDouble(args[1]);
        EpiParams.q = 1- EpiParams.p;
        EpiParams.psi = Double.parseDouble(args[2]);


        model.runSimulation(params);


        for(int j=0; j<n_seasons; j++) {

            System.out.println("season "+(j+1));
            SimulateTree tree = new SimulateTree();
            tree.n_lineages = n_lineages;

            tree.sampleStartTime += j*365.25;
            tree.sampleEndTime += j*365.25;
            tree.getTransmissionTrees(n_segments, model, 0);
        }


        //tree over 40 years

        SimulateTree tree = new SimulateTree();
        tree.n_lineages = 300;
        tree.sampleStartTime = 20*365.25;
        tree.sampleEndTime = 60*365.25;
        tree.getTransmissionTrees(n_segments, model, 0);


    }

    public static void runTwoPatchSim(String [] args)  {

        EpiModel model = new reassortmentTwoPatch_new();
        EpiParams params = new EpiParams();

        EpiParams.psi_i = Double.parseDouble(args[0].trim());
        EpiParams.psi_j =  Double.parseDouble(args[1].trim());
        EpiParams.s_b1 = Double.parseDouble(args[2].trim());
        EpiParams.antigenicMu_a = Double.parseDouble(args[3].trim());



        int n_segments = Integer.parseInt(args[4].trim());
        int n_seasons = Integer.parseInt(args[5].trim());
        int n_lineages = Integer.parseInt(args[6].trim());

        model.runSimulation(params);

        for(int j=0; j<n_seasons; j++) {

            System.out.println("season "+(j+1));
            SimulateTree tree = new SimulateTree();
            tree.sampleSeasonalPatch = true;
            tree.n_lineages = n_lineages;

            tree.sampleStartTime+=j*365.25;
            tree.sampleEndTime+=j*365.25;
            tree.getTransmissionTrees(n_segments, model, 0);
        }

    }

    public static void runMultiSamples(int n_samples, int n_segments, int n_seasons) {

        // make sure that n_seasons is =< than (simulationEndTime-sampleStartTime)years;
        // make sure also that sampleEndTime-sampleStartTime = 365.25 days or 1yr if n_seasons>1

        // when changing between single and two-patch model, make sure to change the sampling scheme
        EpiModel model = new reassortmentTwoPatch_new();
        EpiParams params = new EpiParams();
        model.runSimulation(params);

        for(int i=0; i<n_samples; i++) {

            System.out.println("pop sample no: "+ (i+1));
            for(int j=0; j<n_seasons; j++){
                SimulateTree tree = new SimulateTree();
                tree.sampleStartTime+=j*365.25;
                tree.sampleEndTime+=j*365.25;
                tree.getTransmissionTrees(n_segments, model, (i+1));

            }
        }

    }

    //args - psi, s, antigenicMu

    public static void main(String [] args){


//        EpiParams params = new EpiParams();
//
//        int no_of_sims = 1;
//        double[][] tmrca_seg1 = new double[no_of_sims][61];
//        double[][] tmrca_seg2 = new double[no_of_sims][61];
//        double[][] fitness = new double[no_of_sims][61];
//        double[][] fitnessv = new double[no_of_sims][61];
//
//        reassortmentTauLeap model = new reassortmentTauLeap();
//        model.runSimulation(params, 0, tmrca_seg1, tmrca_seg2, fitness, fitnessv);
//
//        SimulateTree tree = new SimulateTree();
//
//
//        tree.getTransmissionTrees(2, model, 1);

        //runOnePatchSim(args);

        //runOnePatchSim(args);

        //runTwoPatchSim(args);

        //SimulateTree simulateTree = new SimulateTree();

    }


}
