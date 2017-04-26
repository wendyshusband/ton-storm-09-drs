package topology;

import backtype.storm.spout.SpoutOutputCollector;
import backtype.storm.task.TopologyContext;
import backtype.storm.topology.OutputFieldsDeclarer;
import backtype.storm.tuple.Fields;
import backtype.storm.tuple.Values;
import resa.topology.RedisQueueSpout;

import java.util.Map;
import static topology.Constant.*;

/**
 * Created by ding on 14-7-3.
 */
public class ImageSource extends RedisQueueSpout {

    private long frameId;
    private String idPrefix;
    private long acked;

    public ImageSource(String host, int port, String queue) {
        super(host, port, queue, true);
    }

    @Override
    public void open(Map conf, TopologyContext context, SpoutOutputCollector collector) {
        super.open(conf, context, collector);
        this.collector = collector;
        frameId = 0;
        acked = 0;
        idPrefix = String.format("s-%02d-", context.getThisTaskIndex() + 1);
    }

    @Override
    protected void emitData(Object data) {
        String id = idPrefix + frameId++;
        collector.emit(STREAM_IMG_OUTPUT, new Values(id, data), id);
    }

    @Override
    public void declareOutputFields(OutputFieldsDeclarer declarer) {
        declarer.declareStream(STREAM_IMG_OUTPUT, new Fields(FIELD_FRAME_ID, FIELD_IMG_BYTES));
    }

    @Override
    public void ack(Object msgId) {
        acked ++;
//        long unacked = frameId - acked;
        System.out.println(String.format("ImageSource.ack(), sent: %d, acked: %d", frameId, acked));
    }
}
