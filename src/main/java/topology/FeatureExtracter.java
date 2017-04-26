package topology;

import backtype.storm.task.OutputCollector;
import backtype.storm.task.TopologyContext;
import backtype.storm.topology.OutputFieldsDeclarer;
import backtype.storm.topology.base.BaseRichBolt;
import backtype.storm.tuple.Fields;
import backtype.storm.tuple.Tuple;
import backtype.storm.tuple.Values;
import org.bytedeco.javacpp.BytePointer;
import org.bytedeco.javacpp.opencv_core.*;
import org.bytedeco.javacpp.opencv_features2d.KeyPoint;
import org.bytedeco.javacpp.opencv_nonfree.SIFT;
import resa.util.ConfigUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_highgui.cvDecodeImage;
import static topology.Constant.*;

/**
 * Created by ding on 14-7-3.
 */
public class FeatureExtracter extends BaseRichBolt {

    private SIFT sift;
    private double[] buf;
    private OutputCollector collector;

    @Override
    public void prepare(Map stormConf, TopologyContext context, OutputCollector collector) {

        int nfeatures = ConfigUtil.getInt(stormConf, "sift-nfeatures", 0);
        double contrastThreshold = ConfigUtil.getDouble(stormConf, "sift-contrastThreshold", 0.05);
        int edgeThreshold = ConfigUtil.getInt(stormConf, "sift-edgeThreshold", 12);
        sift = new SIFT(nfeatures, 3, contrastThreshold, edgeThreshold, 1.6);
        //sift = new SIFT(500, 3, 0.05, 12, 1.6);
        buf = new double[128];
        this.collector = collector;
    }

    @Override
    public void execute(Tuple input) {
        byte[] imgBytes = (byte[]) input.getValueByField(FIELD_IMG_BYTES);
        IplImage image = cvDecodeImage(cvMat(1, imgBytes.length, CV_8UC1, new BytePointer(imgBytes)));
        KeyPoint points = new KeyPoint();
        Mat featureDesc = new Mat();
        Mat matImg = new Mat(image);
        sift.detect(matImg, points);
        sift.compute(matImg, points, featureDesc);
        try {
            cvReleaseImage(image);
        } catch (Exception e) {
        }
        int rows = featureDesc.rows();
        List<byte[]> selected = new ArrayList<>(rows);
        for (int i = 0; i < rows; i++) {
            featureDesc.rows(i).asCvMat().get(buf);
            // compress data
            byte[] siftFeat = new byte[buf.length];
            for (int j = 0; j < buf.length; j++) {
                siftFeat[j] = (byte) (((int) buf[j]) & 0xFF);
            }
            selected.add(siftFeat);
        }
        String frameId = input.getStringByField(FIELD_FRAME_ID);
        collector.emit(STREAM_FEATURE_DESC, input, new Values(frameId, selected));
        collector.emit(STREAM_FEATURE_COUNT, input, new Values(frameId, selected.size()));
        System.out.println("FrameID: " + frameId + ", rows: " + rows);
        collector.ack(input);
    }

    @Override
    public void declareOutputFields(OutputFieldsDeclarer declarer) {
        declarer.declareStream(STREAM_FEATURE_DESC, new Fields(FIELD_FRAME_ID, FIELD_FEATURE_DESC));
        declarer.declareStream(STREAM_FEATURE_COUNT, new Fields(FIELD_FRAME_ID, FIELD_FEATURE_CNT));
    }
}
