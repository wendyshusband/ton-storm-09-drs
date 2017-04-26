package resa.optimize;

import backtype.storm.Config;
import backtype.storm.generated.StormTopology;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import resa.util.ConfigUtil;
import resa.util.ResaConfig;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by Tom Fu on Feb-15-2016
 */
public class MMKAllocCalculator extends AllocCalculator {
    private static final Logger LOG = LoggerFactory.getLogger(MMKAllocCalculator.class);
    private HistoricalCollectedData spoutHistoricalData;
    private HistoricalCollectedData boltHistoricalData;
    private int historySize;
    private int currHistoryCursor;

    @Override
    public void init(Map<String, Object> conf, Map<String, Integer> currAllocation, StormTopology rawTopology) {
        super.init(conf, currAllocation, rawTopology);
        historySize = ConfigUtil.getInt(conf, ResaConfig.OPTIMIZE_WIN_HISTORY_SIZE, 1);
        currHistoryCursor = ConfigUtil.getInt(conf, ResaConfig.OPTIMIZE_WIN_HISTORY_SIZE_IGNORE, 0);
        spoutHistoricalData = new HistoricalCollectedData(rawTopology, historySize);
        boltHistoricalData = new HistoricalCollectedData(rawTopology, historySize);
    }

    @Override
    public AllocResult calc(Map<String, AggResult[]> executorAggResults, int maxAvailableExecutors) {
        executorAggResults.entrySet().stream().filter(e -> rawTopology.get_spouts().containsKey(e.getKey()))
                .forEach(e -> spoutHistoricalData.putResult(e.getKey(), e.getValue()));
        executorAggResults.entrySet().stream().filter(e -> rawTopology.get_bolts().containsKey(e.getKey()))
                .forEach(e -> boltHistoricalData.putResult(e.getKey(), e.getValue()));
        // check history size. Ensure we have enough history data before we run the optimize function
        currHistoryCursor++;
        if (currHistoryCursor < historySize) {
            LOG.info("currHistoryCursor < historySize, curr: " + currHistoryCursor + ", Size: " + historySize
                    + ", DataHistorySize: "
                    + spoutHistoricalData.compHistoryResults.entrySet().stream().findFirst().get().getValue().size());
            return null;
        } else {
            currHistoryCursor = historySize;
        }

        ///TODO: Here we assume only one spout, how to extend to multiple spouts?
        ///TODO: here we assume only one running topology, how to extend to multiple running topologies?
        double targetQoSMs = ConfigUtil.getDouble(conf, ResaConfig.OPTIMIZE_SMD_QOS_MS, 5000.0);
        double completeTimeMilliSecUpper = ConfigUtil.getDouble(conf, ResaConfig.OPTIMIZE_SMD_QOS_UPPER_MS, 2000.0);
        double completeTimeMilliSecLower = ConfigUtil.getDouble(conf, ResaConfig.OPTIMIZE_SMD_QOS_LOWER_MS, 500.0);
        int maxSendQSize = ConfigUtil.getInt(conf, Config.TOPOLOGY_EXECUTOR_SEND_BUFFER_SIZE, 1024);
        int maxRecvQSize = ConfigUtil.getInt(conf, Config.TOPOLOGY_EXECUTOR_RECEIVE_BUFFER_SIZE, 1024);
        double sendQSizeThresh = ConfigUtil.getDouble(conf, ResaConfig.OPTIMIZE_SMD_SEND_QUEUE_THRESH, 5.0);
        double recvQSizeThreshRatio = ConfigUtil.getDouble(conf, ResaConfig.OPTIMIZE_SMD_RECV_QUEUE_THRESH_RATIO, 0.6);
        double recvQSizeThresh = recvQSizeThreshRatio * maxRecvQSize;
        int resourceUnit = ConfigUtil.getInt(conf, ResaConfig.OPTIMIZE_SMD_RESOURCE_UNIT,1);

        double compSampleRate = ConfigUtil.getDouble(conf, ResaConfig.COMP_SAMPLE_RATE, 1.0);

        Map<String, GeneralSourceNode> spInfos = spoutHistoricalData.compHistoryResults.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> {
                    SpoutAggResult hisCar = AggResult.getHorizontalCombinedResult(new SpoutAggResult(), e.getValue());
                    int numberExecutor = currAllocation.get(e.getKey());
                    double avgSendQLenHis = hisCar.getSendQueueResult().getAvgQueueLength();
                    double avgRecvQLenHis = hisCar.getRecvQueueResult().getAvgQueueLength();

                    double departRateHis = hisCar.getDepartureRatePerSec();
                    //if (ackerIsEnabled)
                    ///TODO: there we multiply 1/2 for this particular implementation
                    ///TODO: this shall be adjusted and configurable for ackering mechanism
                    /// the calculation on tupleEmitRate depends on whether we enable ackering mechanism.
                    double tupleEmitRate = departRateHis * numberExecutor / 2.0;
                    double arrivalRateHis = hisCar.getArrivalRatePerSec();
                    ///For DRS calculation, we instead use externalTupleArrivalRate as the \lambda_0, this will not be affected by ackering issue.
                    double externalTupleArrivalRate = arrivalRateHis * numberExecutor;
                    double tupleEmitRateByInterArrival = hisCar.getSendQueueResult().getAvgArrivalRatePerSecond()* numberExecutor;
                    double tupleEmitInterArrivalScv = hisCar.getSendQueueResult().getScvInterArrivalTimes();
                    double externalRateByInterArrival = hisCar.getRecvQueueResult().getAvgArrivalRatePerSecond()* numberExecutor;
                    double externalTupleInterArrivalScv = hisCar.getRecvQueueResult().getScvInterArrivalTimes();

                    double avgCompleteLatencyHis = hisCar.getCombinedCompletedLatency().getAvg();///unit is millisecond
                    double scvCompleteLatencyHis = hisCar.getCombinedCompletedLatency().getScv();

                    double totalCompleteTupleCnt = hisCar.getCombinedCompletedLatency().getCount();
                    double totalDurationSeconds  = hisCar.getDurationSeconds();
                    double tupleCompleteRate = totalCompleteTupleCnt * numberExecutor / (totalDurationSeconds * compSampleRate);

                    return new GeneralSourceNode(
                            e.getKey(), numberExecutor, compSampleRate, avgSendQLenHis, avgRecvQLenHis, avgCompleteLatencyHis, scvCompleteLatencyHis,
                            totalCompleteTupleCnt, totalDurationSeconds, tupleCompleteRate,
                            tupleEmitRate, tupleEmitRateByInterArrival, tupleEmitInterArrivalScv,
                            externalTupleArrivalRate, externalRateByInterArrival, externalTupleInterArrivalScv);
                }));

        GeneralSourceNode spInfo = spInfos.entrySet().stream().findFirst().get().getValue();

        Map<String, GeneralServiceNode> queueingNetwork = boltHistoricalData.compHistoryResults.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> {
                    BoltAggResult hisCar = AggResult.getHorizontalCombinedResult(new BoltAggResult(), e.getValue());
                    int numberExecutor = currAllocation.get(e.getKey());

                    double avgSendQLenHis = hisCar.getSendQueueResult().getAvgQueueLength();
                    double avgRecvQLenHis = hisCar.getRecvQueueResult().getAvgQueueLength();

                    double avgServTimeHis = hisCar.getCombinedProcessedResult().getAvg();///unit is millisecond
                    double scvServTimeHis = hisCar.getCombinedProcessedResult().getScv();

                    double arrivalRateHis = hisCar.getArrivalRatePerSec();
                    double lambdaHis = arrivalRateHis * numberExecutor;
                    double arrivalByInterArrival = hisCar.getRecvQueueResult().getAvgArrivalRatePerSecond() * numberExecutor;
                    double interArrivalScv = hisCar.getRecvQueueResult().getScvInterArrivalTimes();

                    double totalProcessTupleCnt = hisCar.getCombinedProcessedResult().getCount();
                    double totalDurationSecond = hisCar.getDurationSeconds();
                    double tupleProcessRate = totalProcessTupleCnt * numberExecutor / (totalDurationSecond * compSampleRate);

                    ///TODO: here i2oRatio can be INFINITY, when there is no data sent from Spout.
                    ///TODO: here we shall decide whether to use external Arrival rate, or tupleLeaveRateOnSQ!!
                    ///TODO: major differences 1) when there is max-pending control, tupleLeaveRateOnSQ becomes the
                    ///TODO: the tupleEmit Rate, rather than the external tuple arrival rate (implicit load shading)
                    ///TODO: if use tupleLeaveRateOnSQ(), be careful to check if ACKing mechanism is on, i.e.,
                    ///TODO: there are ack tuple. otherwise, divided by tow becomes meaningless.
                    ///TODO: shall we put this i2oRatio calculation here, or later to inside ServiceModel?
                    ///double i2oRatio = lambdaHis / spInfo.getTupleLeaveRateOnSQ();

                    return new GeneralServiceNode(
                            e.getKey(), numberExecutor, compSampleRate, avgSendQLenHis, avgRecvQLenHis,
                            avgServTimeHis, scvServTimeHis,
                            totalProcessTupleCnt, totalDurationSecond, tupleProcessRate,
                            lambdaHis, arrivalByInterArrival, interArrivalScv,
                            spInfo.getExArrivalRate(), spInfo.getExArrivalRateByInterArrival());
                }));

        Map<String, Integer> boltAllocation = currAllocation.entrySet().stream()
                .filter(e -> rawTopology.get_bolts().containsKey(e.getKey()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        ///totalAvailableExecutors - spoutExecutors, currently, it is assumed that there is only one spout
        ///Kmax
        int maxThreadAvailable4Bolt = maxAvailableExecutors - currAllocation.entrySet().stream()
                .filter(e -> rawTopology.get_spouts().containsKey(e.getKey()))
                .mapToInt(Map.Entry::getValue).sum();

        int currentUsedThreadByBolts = currAllocation.entrySet().stream()
                .filter(e -> rawTopology.get_bolts().containsKey(e.getKey())).mapToInt(Map.Entry::getValue).sum();

        LOG.info("Run Optimization, tQos: " + targetQoSMs + ", currUsed: " + currentUsedThreadByBolts + ", kMax: " + maxThreadAvailable4Bolt + ", currAllo: " + currAllocation);
        AllocResult allocResult = MMKServiceModel.checkOptimized(
                spInfo, queueingNetwork, completeTimeMilliSecUpper, completeTimeMilliSecLower, boltAllocation, maxThreadAvailable4Bolt, currentUsedThreadByBolts, resourceUnit);


        Map<String, Integer> retCurrAllocation = null;
        if (allocResult.currOptAllocation != null) {
            retCurrAllocation = new HashMap<>(currAllocation);
            retCurrAllocation.putAll(allocResult.currOptAllocation);
        }
        Map<String, Integer> retKMaxAllocation = null;
        if (allocResult.kMaxOptAllocation != null) {
            retKMaxAllocation = new HashMap<>(currAllocation);
            retKMaxAllocation.putAll(allocResult.kMaxOptAllocation);
        }
        Map<String, Integer> retMinReqAllocation = null;
        if (allocResult.minReqOptAllocation != null) {
            retMinReqAllocation = new HashMap<>(currAllocation);
            retMinReqAllocation.putAll(allocResult.minReqOptAllocation);
        }
        Map<String, Object> ctx = new HashMap<>();
        ctx.put("latency", allocResult.getContext());
        ctx.put("spout", spInfo);
        ctx.put("bolt", queueingNetwork);
        return new AllocResult(allocResult.status, retMinReqAllocation, retCurrAllocation, retKMaxAllocation).setContext(ctx);
    }

    @Override
    public void allocationChanged(Map<String, Integer> newAllocation) {
        super.allocationChanged(newAllocation);
        spoutHistoricalData.clear();
        boltHistoricalData.clear();
        currHistoryCursor = ConfigUtil.getInt(conf, ResaConfig.OPTIMIZE_WIN_HISTORY_SIZE_IGNORE, 0);
    }
}
