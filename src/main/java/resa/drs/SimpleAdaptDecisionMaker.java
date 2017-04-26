package resa.drs;

import backtype.storm.generated.StormTopology;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import resa.optimize.AllocResult;
import resa.util.ConfigUtil;
import resa.util.ResaConfig;

import java.util.Map;

import static resa.util.ResaConfig.OPTIMIZE_INTERVAL;

/**
 * Created by ding on 14-6-20.
 */
public class SimpleAdaptDecisionMaker implements DecisionMaker {

    private static final Logger LOG = LoggerFactory.getLogger(SimpleAdaptDecisionMaker.class);

    long startTimeMillis;
    long minExpectedIntervalMillis;
    int rbTypeValue;

    @Override
    public void init(Map<String, Object> conf, StormTopology rawTopology) {
        startTimeMillis = System.currentTimeMillis();
        long calcIntervalSec = ConfigUtil.getInt(conf, OPTIMIZE_INTERVAL, 30);
        /** if OPTIMIZE_MIN_EXPECTED_REBALANCE_INTERVAL is not found in configuration files,
         * we use twice of OPTIMIZE_INTERVAL as default
         * here we -50 for synchronization purpose, this needs to be tested **/
        minExpectedIntervalMillis = ConfigUtil.getLong(conf, ResaConfig.OPTIMIZE_MIN_EXPECTED_REBALANCE_INTERVAL, calcIntervalSec * 2) * 1000 - 50;
//        rbTypeValue = ConfigUtil.getInt(conf, ResaConfig.OPTIMIZE_REBALANCE_TYPE, RebalanceType.CurrentOpt.getValue());

        LOG.info("SimpleAdaptDecisionMaker.init(), stTime: " + startTimeMillis + ", minExpInteval: " + minExpectedIntervalMillis);
    }

    @Override
    public Map<String, Integer> make(AllocResult newAllocResult, Map<String, Integer> currAlloc) {

        long timeSpan = Math.max(0, System.currentTimeMillis() - startTimeMillis);
        if (newAllocResult == null) {
            LOG.info("SimpleAdaptDecisionMaker.make(), newAllocResult == null");
            return null;
        }


        if (timeSpan < minExpectedIntervalMillis) {
            /** if  timeSpan is not large enough, no rebalance will be triggered **/
            LOG.info("SimpleAdaptDecisionMaker.make(), timeSpan (" + timeSpan + ") < minExpectedIntervalMillis (" + minExpectedIntervalMillis + ")");
            return null;
        } else {

            if (newAllocResult.status.equals(AllocResult.Status.OVERPROVISIONING)) {

                LOG.info("SimpleAdaptDecisionMaker.make(), ewAllocResult.status == OVERPROVISIONING, rebalance is triggered with removing existing resources");
                return newAllocResult.minReqOptAllocation;

            } else if (newAllocResult.status.equals(AllocResult.Status.SHORTAGE)) {

                LOG.info("SimpleAdaptDecisionMaker.make(), ewAllocResult.status == OVERPROVISIONING, rebalance is triggered with adding new resources");
                return newAllocResult.minReqOptAllocation;

            } else if (newAllocResult.status.equals(AllocResult.Status.INFEASIBLE)) {

                LOG.info("SimpleAdaptDecisionMaker.make(), ewAllocResult.status == INFEASIBLE, rebalance is triggered with using maximum available resources");
                return newAllocResult.kMaxOptAllocation;

            } else if (newAllocResult.status.equals(AllocResult.Status.FEASIBLE)) {

                LOG.info("SimpleAdaptDecisionMaker.make(), " +
                        "ewAllocResult.status == FEASIBLE, rebalance is triggered without adjusting current used resources, but re-allocation to optimal may be possible");
                ///return newAllocResult.currOptAllocation;
                return null;

            } else {
                throw new IllegalArgumentException("Illegal status: " + newAllocResult.status);
            }
        }
    }
}
