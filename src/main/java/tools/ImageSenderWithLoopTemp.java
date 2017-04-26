package tools;

import redis.clients.jedis.Jedis;
import resa.util.ConfigUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Map;

/**
 * Created by ding on 14-3-18.
 */
public class ImageSenderWithLoopTemp {

    private static final File END = new File("END");

    private String host;
    private int port;
    private byte[] queueName;
    private String path;
    private String imageFolder;
    private String filePrefix;

    public ImageSenderWithLoopTemp(Map<String, Object> conf) {
        this.host = (String) conf.get("redis.host");
        this.port = ((Number) conf.get("redis.port")).intValue();
        this.queueName = ((String) conf.get("redis.queue")).getBytes();
        this.path = (String) conf.get("sourceFilePath");
        this.imageFolder = (String) conf.get("imageFolder");
        this.filePrefix = (String) conf.get("filePrefix");

    }

    public void send2Queue(int st, int end, int fps1, int fps2, double p, int safeLen, int stopTime, int totalDurationInSeconds) throws IOException {
        Jedis jedis = new Jedis(host, port);
        int generatedFrames = 0;
        int fileIndex = st;

        try {
            final long startTime = System.currentTimeMillis();
            long last = startTime;
            long qLen = 0;
            int target = Math.random() > p ? fps1 : fps2;
            int finished = 0;

            while (true) {

                String fileName = path + imageFolder + System.getProperty("file.separator")
                        + String.format("%s%06d.jpg", filePrefix, (fileIndex++));
                File f = new File(fileName);
                if (f.exists() == false) {
                    System.out.println("File not exist: " + fileName);
                    continue;
                }

                jedis.rpush(this.queueName, Files.readAllBytes(new File(fileName).toPath()));

                if (fileIndex > end) {
                    fileIndex = st;
                }

                generatedFrames++;
                if (++finished == target) {
                    long current = System.currentTimeMillis();
                    long elapse = current - last;
                    long remain = 1000 - elapse;
                    if (remain > 0) {
                        Thread.sleep(remain);
                    }
                    last = System.currentTimeMillis();
                    qLen = jedis.llen(this.queueName);
                    System.out.println("Target: " + target + ", elapsed: " + (last - startTime) / 1000
                            + ",totalSend: " + generatedFrames + ", remain: " + remain + ", sendQLen: " + qLen);
                    while (qLen > safeLen) {
                        Thread.sleep(Math.max(100, stopTime));
                        System.out.println("qLen > safeLen, stop sending....");
                        qLen = jedis.llen(this.queueName);
                    }
                    target = Math.random() > p ? fps1 : fps2;
                    finished = 0;
                }

                if (totalDurationInSeconds > 0) {
                    int totalLast = (int) ((System.currentTimeMillis() - startTime) / 1000);
                    if (totalLast > totalDurationInSeconds) {
                        System.out.println(totalDurationInSeconds + " seconds passed and exit");
                        return;
                    }
                }
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 9) {
            System.out.println("usage: ImageSender <confFile> <fileSt> <fileEnd> <fps> <range> <safeLen> <stopTime> <mode> <totalDuration, -1 means infinity>");
            return;
        }
        ImageSenderWithLoopTemp sender = new ImageSenderWithLoopTemp(ConfigUtil.readConfig(new File(args[0])));
        System.out.println(String.format("start sender, st: %d, end: %d, fps1: %d, fps2: %d, p: %.3f, slen: %d, sTime: %d, totalDuration: %d",
                Integer.parseInt(args[1]), Integer.parseInt(args[2]),
                Integer.parseInt(args[3]), Integer.parseInt(args[4]), Double.parseDouble(args[5]),
                Integer.parseInt(args[6]), Integer.parseInt(args[7]),
                Integer.parseInt(args[8])));
        sender.send2Queue(Integer.parseInt(args[1]), Integer.parseInt(args[2]),
                Integer.parseInt(args[3]), Integer.parseInt(args[4]), Double.parseDouble(args[5]),
                Integer.parseInt(args[6]), Integer.parseInt(args[7]),
                Integer.parseInt(args[8]));
        System.out.println("end sender");
    }

}
