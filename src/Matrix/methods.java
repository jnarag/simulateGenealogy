package Matrix;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 1/15/13
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class methods {

    public int[] getColumnArray(List<int[]> I_matrix, int columnNo) {

        int[] columnVector = new int[I_matrix.size()];

        List<Integer> vector  = new ArrayList<Integer>();
        for(int i = 0; i < I_matrix.size(); i++) {

            int x = I_matrix.get(i)[columnNo];

            columnVector[i] = x;
            vector.add(x);
        }

        Collections.sort(vector);

        return columnVector;
    }

    public double[] getTimeVector(List<int[]> I_matrix, int columnNo) {

        int[] columnVector = getColumnArray(I_matrix, columnNo);

        double[] timeVector = new double[columnVector.length];

        int j = 0;

        double x = Double.NaN;
        for(int i: columnVector) {

            if(i >= 0) {

                x = (double)i/(60*60*24);
            }


            timeVector[j] = x;

            j++;
        }
        return timeVector;
    }

//    public int[] getIndices(List<int[]> matrix) {
//
//        int[] indices = new int[matrix.size()];
//        int i = 0;
//        while(i<matrix.size()){
//
//            indices[i] = i;
//            i++;
//
//        }
//        return indices;
//    }

//    public List<Integer> getIndices(List<Integer> list) {
//
//        List<Integer> indices = new ArrayList<Integer>();
//        int i = 0;
//
//        while(i<list.size()){
//
//            indices.add(i);
//            i++;
//
//        }
//        return indices;
//    }

//    public List<Integer> getIndices(double[] array) {
//
//        List<Integer> indices = new ArrayList<Integer>();
//        for(int i = 0; i<array.length; i++) {
//
//            indices.add(i);
//        }
//        return indices;
//    }

    public List<Integer> getIndices(List list) {

        List<Integer> indices = new ArrayList<Integer>();
        int i = 0;


        while(i<list.size()){

            indices.add(i);
            i++;

        }
        return indices;
    }

    public List<Integer> find(double[] array, int condition, double value) {


        List<Integer> subList = new ArrayList<Integer>();

        switch(condition) {


            case -2: //find elements that are equal to -1;

                for(int i=0; i<array.length; i++) {

                    if(array[i]<0.0){ //(-1.0/(60*60*24))){                 //why have I got this line here? why is it -1.0/(60*60*24)      //check if this is used AT ALL

                        subList.add(i);
                    }
                }

                break;
            case -1: //find elements less than the value

                for(int i=0; i<array.length; i++) {

                    if(array[i]<value){

                        //System.out.println("lessthan"+array[i]);
                        subList.add(i);
                    }
                }
                break;
            case 0: //find elements equal to the value

                for(int i=0; i<array.length; i++) {

                    if(array[i]==value){

                        subList.add(i);
                    }
                }
                break;
            case 1: //find elements larger than the value

                for(int i=0; i<array.length; i++) {

                    if(array[i]>value){
                        subList.add(i);
                    }
                }
                break;

            case 2: //find elements isNaN

                for(int i=0; i<array.length; i++) {

                    if(Double.isInfinite(array[i])){
                        subList.add(i);
                    }
                }
            case 3: //find elements !isNaN

                for(int i=0; i<array.length; i++) {

                    if(!Double.isInfinite(array[i])){
                        subList.add(i);
                    }
                }
            case 4:

                for(int i=0; i<array.length; i++) {

                    if(array[i]<value){
                        subList.add(i);
                    }

                    if(Double.isInfinite(array[i])){
                        subList.add(i);
                    }
                }

            case 5:

                for(int i=0; i<array.length; i++) {

                    if(array[i]>value){
                        subList.add(i);
                    }

                    if(Double.isInfinite(array[i])){
                        subList.add(i);
                    }
                }
        }

        return subList;

    }

    public List<Integer> find(List<Double> array, int condition, double value) {


        List<Integer> subList = new ArrayList<Integer>();

        switch(condition) {


            case -2: //find elements that are equal to -1;

                for(int i=0; i<array.size(); i++) {

                    if(array.get(i)<0.0){ //(-1.0/(60*60*24))){                 //why have I got this line here? why is it -1.0/(60*60*24)      //check if this is used AT ALL

                        subList.add(i);
                    }
                }

                break;
            case -1: //find elements less than the value

                for(int i=0; i<array.size(); i++) {

                    if(array.get(i)<value){

                        //System.out.println("lessthan"+array[i]);
                        subList.add(i);
                    }
                }
                break;
            case 0: //find elements equal to the value

                for(int i=0; i<array.size(); i++) {

                    if(array.get(i)==value){

                        subList.add(i);
                    }
                }
                break;
            case 1: //find elements larger than the value

                for(int i=0; i<array.size(); i++) {

                    if((Double)array.get(i)>value){
                        subList.add(i);
                    }
                }
                break;

            case 2: //find elements isNaN

                for(int i=0; i<array.size(); i++) {

                    if(Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }
            case 3: //find elements !isNaN

                for(int i=0; i<array.size(); i++) {

                    if(!Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }
            case 4:

                for(int i=0; i<array.size(); i++) {

                    if(array.get(i)<value){
                        subList.add(i);
                    }

                    if(Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }

            case 5:

                for(int i=0; i<array.size(); i++) {

                    if(array.get(i)>value){
                        subList.add(i);
                    }

                    if(Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }
        }

        return subList;

    }

    public List<Integer> find(Map<Integer, Double> array, int condition, double value) {


        List<Integer> subList = new ArrayList<Integer>();

        switch(condition) {


            case -2: //find elements that are equal to -1;

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {

                    if(array.get(i)<0.0){ //(-1.0/(60*60*24))){                 //why have I got this line here? why is it -1.0/(60*60*24)      //check if this is used AT ALL

                        subList.add(i);
                    }
                }

                break;
            case -1: //find elements less than the value

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if(array.get(i)<value){

                        //System.out.println("lessthan"+array[i]);
                        subList.add(i);
                    }
                }
                break;
            case 0: //find elements equal to the value

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if(array.get(i)==value){

                        subList.add(i);
                    }
                }
                break;
            case 1: //find elements larger than the value

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if((Double)array.get(i)>value){
                        subList.add(i);
                    }
                }
                break;

            case 2: //find elements isNaN

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if(Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }
            case 3: //find elements !isNaN

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if(!Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }
            case 4:

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if(array.get(i)<value){
                        subList.add(i);
                    }

                    if(Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }

            case 5:

                //for(int i=0; i<array.size(); i++) {

                for(Integer i: array.keySet()) {


                    if(array.get(i)>value){
                        subList.add(i);
                    }

                    if(Double.isInfinite(array.get(i))){
                        subList.add(i);
                    }
                }
        }

        return subList;

    }

    public List<Integer> findAliveLocs(double[] array1, double[] array2, double value) {

        List<Integer> locs = new ArrayList<Integer>();

        int arrayLength = array1.length;

        for(int i=0; i < arrayLength; i++) {

            if(array1[i] < value && array2[i] > value ) {


                locs.add(i);

            }
            if(Double.isInfinite(array1[i]) && array2[i] > value ) {

                locs.add(i);

            }
        }

        return locs;
    }

    public List<Integer> findAliveLocs(List<Double> array1, List<Double> array2, double value) {

        Set<Integer> locs = new HashSet<Integer>();

        int arrayLength = array1.size();

        for(int i=0; i < arrayLength; i++) {

            if(array1.get(i) < value && array2.get(i) > value ) {


                locs.add(i);

            }
            if(Double.isInfinite(array1.get(i)) && array2.get(i) > value ) {

                locs.add(i);

            }
        }

        List<Integer> locsList = new ArrayList<Integer>();
        locsList.addAll(locs);
        return locsList;
    }

    public List<Integer> findAliveLocsX(Map<Integer, Double> array1, Map<Integer, Double> array2, double value) {

        List<Integer> locs = new ArrayList<Integer>();

        //int arrayLength = array1.size();

        Set<Integer> keySet = array1.keySet();
        //for(int i=0; i < arrayLength; i++) {
        for(Integer i: keySet) {
            if(array1.get(i) < value && array2.get(i) > value ) {
                //System.out.println(array1.get(i)+ " : "+array2.get(i));

                locs.add(i);

            }


            if(Double.isInfinite(array1.get(i)) && array2.get(i) > value ) {

                locs.add(i);

            }
        }

        return locs;
    }

    public Set<Integer> findAliveLocs(Map<Integer, Double> array1, Map<Integer, Double> array2, double value) {

        Set<Integer> locs = new HashSet<Integer>();

        //int arrayLength = array1.size();

        Set<Integer> keySet = array1.keySet();
        //for(int i=0; i < arrayLength; i++) {
        for(Integer i: keySet) {
            if(array1.get(i) < value && array2.get(i) > value ) {
                //System.out.println(array1.get(i)+ " : "+array2.get(i));

                locs.add(i);

            }
            if(Double.isInfinite(array1.get(i)) && array2.get(i) > value ) {

                locs.add(i);

            }
        }

        return locs;
    }

    public List<Integer> find(int[] array, int condition, int value) {

        //retains the indices whose elements certify the condition
        List<Integer> subList = new ArrayList<Integer>();

        switch(condition) {


            case -2: //find elements that are equal to -1;
                for(int i=0; i<array.length; i++) {

                    if(array[i]==-1){

                        subList.add(i);
                    }
                }

                break;
            case -1: //find elements less than the value

                for(int i=0; i<array.length; i++) {

                    if(array[i]<value){

                        subList.add(i);
                    }
                }
                break;
            case 0: //find elements equal to the value
                for(int i=0; i<array.length; i++) {

                    if(array[i]==value){

                        subList.add(i);
                    }
                }
                break;
            case 1: //find elements larger than the value

                for(int i=0; i<array.length; i++) {

                    if(array[i]>value){

                        subList.add(i);
                    }
                }
                break;


        }

        return subList;
    }

    public List<Integer> find(List<Integer> list, int value) {

        List<Integer> indices = new ArrayList<Integer>();
        List<Integer> copyList = new ArrayList<Integer>();
        copyList.addAll(list);
        Integer intObject = Integer.valueOf(value);
        while(copyList.contains(intObject)) {

            int index = copyList.indexOf(intObject);
            indices.add(index);
            copyList.set(index,(int)Double.NEGATIVE_INFINITY);
        }

        return indices;

    }

    public List<Integer> find(double[] array, double value) {


        List<Integer> indices = new ArrayList<Integer>();
        double[] arrayCopy = array;
        while(Arrays.binarySearch(arrayCopy, value) >= 0) {

            int index = Arrays.binarySearch(arrayCopy, value);
            indices.add(index);
            arrayCopy[index] = Double.NaN;

        }

        return indices;

    }

    public int find(List<Integer> list, Integer value) {

        return list.indexOf(value);

    }

    public List<Integer> find(List<Double> list, Double value) {

        List<Integer> indices = new ArrayList<Integer>();

        List<Double> copyList = new ArrayList<Double>(list);
        while(copyList.contains(value)) {

            int index = copyList.indexOf(value);
            indices.add(index);
            copyList.set(index,Double.NaN);
            //copyList.remove(index);

        }

        return indices;

    }

    public List<Integer> union(List<Integer> list1, List<Integer> list2) {

        Set<Integer> lists1and2 = new HashSet<Integer>();
        lists1and2.addAll(list1);
        lists1and2.addAll(list2);

        List<Integer> unionList = new ArrayList<Integer>();
        unionList.addAll(lists1and2);

        return unionList;

    }

//    public List<Integer> intersect(List<Integer> list1, List<Integer> list2) {
//
//        BitSet set1 = listToBitSet(list1);
//        BitSet set2 = listToBitSet(list2);
//        BitSet onlySet1 = (BitSet)set1.clone();
//
//        onlySet1.andNot(set2);
//        set1.andNot(onlySet1);
//
//        List<Integer> intersectList = bitSetToList(set1);
//
//        return intersectList;
//
//
//    }

    public List<Integer> intersect(List<Integer> list1, List<Integer> list2) {


        Set<Integer> set1 = new HashSet<Integer>();
        Set<Integer> set2 = new HashSet<Integer>();

        set1.addAll(list1);
        set2.addAll(list2);
        //Collections.copy(newList, list1);

        set1.retainAll(set2);

        List<Integer> intersectList = new ArrayList<Integer>();
        intersectList.addAll(set1);
        return intersectList;
    }

    public List<Integer> intersect(Set<Integer> set1, Set<Integer> set2) {

        List<Integer> list = new ArrayList<Integer>();
        set1.retainAll(set2);
        list.addAll(set1);
        return list;
    }

    public BitSet listToBitSet(List<Integer> list) {

        BitSet set = new BitSet();
        for(Integer i: list) {

            set.set(i);
        }
        return set;
    }

    public List<Integer> bitSetToList(BitSet set) {

        List<Integer> list = new ArrayList<Integer>();
        for(int i=0; i<set.size(); i++) {

            if(set.get(i)){

                list.add(i);

            }
        }
        return list;

    }

    public List<Integer> setDiff(List<Integer> list, List<Integer> sampleList) {

        List<Integer> tempList = new ArrayList<Integer>();

        BitSet set1 = listToBitSet(list);
        BitSet set2 = listToBitSet(sampleList);

        set1.andNot(set2);
        tempList = bitSetToList(set1);

        return tempList;

    }

    public List<Integer> sortListsByTime(List<Double> sampleTimes, List<Integer> sampledLineages) {


        List<Double> originalTimeArray = new ArrayList<Double>(sampleTimes);

        Collections.sort(sampleTimes);

        List<Integer> sortedSampledLineages = new ArrayList<Integer>();

        int i = 0;
        while(sortedSampledLineages.size()<sampledLineages.size()) {


            List<Integer> indicesInOriginal = find(originalTimeArray, sampleTimes.get(i));

            if(indicesInOriginal.size()==1) {
                Integer sampleLineage = sampledLineages.get(indicesInOriginal.get(0));
                sortedSampledLineages.add(sampleLineage);
                i++;
            }
            else{
                int j= 0;
                while(j < indicesInOriginal.size()) {

                    Integer sampleLineage = sampledLineages.get(indicesInOriginal.get(j));
                    sortedSampledLineages.add(sampleLineage);
                    j++;
                    i++;

                }

            }

        }

        return sortedSampledLineages;

    }


    public List<String> createNames(List<Integer> sampledLineages, List<Double> sampleTimes) {


        List<String> names = new ArrayList<String>();

        for(int i=0; i<sampledLineages.size(); i++) {

            String name = "sample_"+(i+1)+"_"+sampleTimes.get(i).toString();
            names.add(name);
        }

        return names;
    }

    public void labelReassortants(List<String> names, List<Integer> sampledLineages, List<Integer> reassortants, List<Integer> parent_co, List<Integer> patch) {

        //L//ist<String> newNames = new ArrayList<String>();

        for(int i=0; i<sampledLineages.size(); i++) {


            Integer sample = sampledLineages.get(i);


            int patch_loc;
            if(patch==null) {
                patch_loc = (int)Double.NEGATIVE_INFINITY;
            }
            else{
                patch_loc = patch.get(sample);
            }

            //System.out.println("x "+sample+" : "+patch_loc+" : "+reassortants[sample]+ " : "+parent_co[sample]);
            String name = "patch_"+patch_loc+"_"+names.get(i);
            String newName = name;
            if(reassortants.get(sample)==1) {

                //newName = newName+"_sample_r";
                if(parent_co.get(sample)>0) {
                    newName = name.replaceAll("sample_","r_prim_sample_");
                }
                else{
                    newName = name.replaceAll("sample_","r_seco_sample_");
                }
            }
            else{
                newName = name.replaceAll("sample_","nr_sample_");
            }

            names.set(i, newName);


        }

    }

    //copies the array for the corresponding indices
    public double[] copyArrayRangeOf(double[] array, List<Integer> indices) {

        double[] subArray = new double[indices.size()];

        for(int i=0; i<indices.size(); i++) {

            subArray[i] = array[indices.get(i)];

        }

        return subArray;
    }

    public double[] copyArrayRangeOf(List<Double> array, List<Integer> indices) {

        double[] subArray = new double[indices.size()];

        for(int i=0; i<indices.size(); i++) {

            subArray[i] = array.get(indices.get(i));

        }

        return subArray;
    }

    public double[] copyArrayRangeOf(List<Double> array, Set<Integer> indices) {

        double[] subArray = new double[indices.size()];

        int j = 0;
        for(Integer i: indices) {

            subArray[j] = array.get(i);
            j++;
        }

        return subArray;
    }

    public double[] copyArrayRangeOf(Map<Integer, Double> array, List<Integer> indices) {

        double[] subArray = new double[indices.size()];

        for(int i=0; i<indices.size(); i++) {


            subArray[i] = array.get(indices.get(i));

        }

        return subArray;
    }

    public int lastIndex(List<Integer> list) {

        int lastIndex = list.size()-1;


        return lastIndex;
    }

    public int lastIndex(double[] array) {

        int lastIndex = array.length-1;
        return lastIndex;
    }

    public void isReassortant(List<Integer> sampledLineages, List<Integer> reassortants, List<Integer> parent_co, List<List<Integer>> parentLineages)  {


        for(int i = 0; i< sampledLineages.size(); i++) {

            if(reassortants.get(sampledLineages.get(i))==1 && parent_co.get(sampledLineages.get(i))>=0) {

                System.out.println("1 reassortment "+sampledLineages.get(i));

            }
            else if(reassortants.get(sampledLineages.get(i))==1 && parent_co.get(sampledLineages.get(i))<0) {


                System.out.println("2 reassortment "+sampledLineages.get(i));

            }

        }

    }

    public double[] listToArray(List<Double> list) {

        double[] array = new double[list.size()];
        for(int i=0; i<list.size(); i++) {

            array[i] = list.get(i);
        }
        return array;
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

                list.set((n_lineages+loc_b), ">"+listToAppend.get(i));
                loc_b++;

            }
            else{
                list.add((n_lineages+loc_b), "<"+listToAppend.get(i));
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
}
