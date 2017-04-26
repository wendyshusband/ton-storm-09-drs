package resa.optimize;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Created by Tom.fu on 23/4/2014.
 */
public class GeneralServiceNode {

    private static final Logger LOG = LoggerFactory.getLogger(GeneralServiceNode.class);

    protected String componentID;
    protected int executorNumber;
    protected double compSampleRate;

    protected double avgSendQueueLength;
    protected double avgRecvQueueLength;

    protected double avgServTimeHis;
    protected double scvServTimeHis;
    protected double mu;

    protected double numCompleteTuples;
    protected double sumDurationSeconds;
    protected double tupleCompleteRate;

    /*metrics on recv_queue*/
    protected double lambda;
    protected double lambdaByInterArrival;
    protected double interArrivalScv;

    protected double exArrivalRate;
    protected double exArrivalRateByInterArrival;

    protected double ratio;
    protected double ratioByInterArrival;

    protected double rho = lambda / (executorNumber * mu);
    protected double rhoBIA = lambdaByInterArrival / (executorNumber * mu);

    public GeneralServiceNode(
            String componentID,
            int executorNumber,
            double compSampleRate,
            double avgSendQueueLength,
            double avgRecvQueueLength,
            double avgServTimeHis,
            double scvServTimeHis,
            double numCompleteTuples,
            double sumDurationSeconds,
            double tupleCompleteRate,
            double lambda,
            double lambdaByInterArrival,
            double interArrivalScv,
            double exArrivalRate,
            double exArrivalRateByInterArrival) {
        this.componentID = componentID;
        this.executorNumber = executorNumber;
        this.compSampleRate = compSampleRate;
        this.avgSendQueueLength = avgSendQueueLength;
        this.avgRecvQueueLength = avgRecvQueueLength;
        this.avgServTimeHis = avgServTimeHis;
        this.scvServTimeHis = scvServTimeHis;
        this.numCompleteTuples = numCompleteTuples;
        this.sumDurationSeconds = sumDurationSeconds;
        this.tupleCompleteRate = tupleCompleteRate;
        this.lambda = lambda;
        this.lambdaByInterArrival = lambdaByInterArrival;
        this.interArrivalScv = interArrivalScv;
        this.exArrivalRate = exArrivalRate;
        this.exArrivalRateByInterArrival = exArrivalRateByInterArrival;

        this.mu = this.avgServTimeHis > 0.0 ? (1000.0 / this.avgServTimeHis) : Double.MAX_VALUE;
        this.ratio = this.exArrivalRate > 0.0 ? (this.lambda / this.exArrivalRate) : 0;
        this.ratioByInterArrival = this.exArrivalRateByInterArrival > 0.0 ? (this.lambdaByInterArrival / this.exArrivalRateByInterArrival) : 0;

        rho = this.lambda / (this.executorNumber * mu);
        rhoBIA = this.lambdaByInterArrival / (this.executorNumber * mu);

        LOG.info(toString());
    }

    @Override
    public String toString() {
//        return String.format(
//                "Component(ID, eNum):(%s,%d), tupleProcCnt: %.1f, sumMeasuredDur: %.1f, sampleRate: %.1f, tupleProcRate: %.3f, " +
//                        "avgSendQLen: %.1f, avgRecvQLen: %.1f, avgServTimeMS: %.3f, scvServTime: %.3f, mu: %.3f, " +
//                        "arrRateHis: %.3f, arrRateBIA: %.3f, interArrivalScv: %.3f, " +
//                        "ratio: %.3f, ratioBIA: %.3f, rho: %.3f, rhoBIA: %.3f",
//                componentID, executorNumber, numCompleteTuples, sumDurationSeconds, compSampleRate, tupleCompleteRate,
//                avgSendQueueLength, avgRecvQueueLength, avgServTimeHis, scvServTimeHis, mu,
//                lambda, lambdaByInterArrival, interArrivalScv,
//                ratio, ratioByInterArrival, rho, rhoBIA);

        return String.format(
                "(ID, eNum):(%s,%d), ProcRate: %.3f, avgSTime: %.3f, scvSTime: %.3f, mu: %.3f, ProcCnt: %.1f, Dur: %.1f, sample: %.1f, SQLen: %.1f, RQLen: %.1f, " +
                "-----> arrRate: %.3f, arrRateBIA: %.3f, arrScv: %.3f, ratio: %.3f, ratioBIA: %.3f, rho: %.3f, rhoBIA: %.3f",
                componentID, executorNumber, tupleCompleteRate, avgServTimeHis, scvServTimeHis, mu,
                numCompleteTuples, sumDurationSeconds, compSampleRate, avgSendQueueLength, avgRecvQueueLength,
                lambda, lambdaByInterArrival, interArrivalScv, ratio, ratioByInterArrival, rho, rhoBIA);
    }

    public String getComponentID() {
        return componentID;
    }

    public int getExecutorNumber() {
        return executorNumber;
    }

    public double getCompSampleRate() {
        return compSampleRate;
    }

    public double getAvgSendQueueLength() {
        return avgSendQueueLength;
    }

    public double getAvgRecvQueueLength() {
        return avgRecvQueueLength;
    }

    public double getAvgServTimeHis() {
        return avgServTimeHis;
    }

    public double getScvServTimeHis() {
        return scvServTimeHis;
    }

    public double getNumCompleteTuples() {
        return numCompleteTuples;
    }

    public double getSumDurationSeconds() {
        return sumDurationSeconds;
    }

    public double getTupleCompleteRate() {
        return tupleCompleteRate;
    }

    public double getLambda() {
        return lambda;
    }

    public double getLambdaByInterArrival() {
        return lambdaByInterArrival;
    }

    public double getInterArrivalScv() {
        return interArrivalScv;
    }

    public double getExArrivalRate() {
        return exArrivalRate;
    }

    public double getExArrivalRateByInterArrival() {
        return exArrivalRateByInterArrival;
    }

    public double getMu() {
        return mu;
    }

    public double getRatio() {
        return ratio;
    }

    public double getRatioByInterArrival() {
        return ratioByInterArrival;
    }

    public double getRho() {
        return rho;
    }

    public double getRhoBIA() {
        return rhoBIA;
    }
}
