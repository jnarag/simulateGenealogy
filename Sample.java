import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 11/28/12
 * Time: 4:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class Sample {

    List<Integer> lineages;
    List<Double> times;
    List<Double> sampleBirths;
    List<Double> sampleDeaths;

    public Sample() {

        lineages = new ArrayList<Integer>();
        times = new ArrayList<Double>();

    }

    public void addSamples(List<Integer> lineages, List<Double> times) {

        this.lineages.addAll(lineages);
        this.times.addAll(times);

    }

    public void addSamples(Sample sample) {

        this.lineages.addAll(sample.getLineages());
        this.times.addAll(sample.getTimes());
    }

    public List<Integer> getLineages() {

        return lineages;
    }

    public List<Double> getTimes() {

        return times;
    }
}
