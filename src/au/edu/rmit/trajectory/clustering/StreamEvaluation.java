package au.edu.rmit.trajectory.clustering;

import java.io.IOException;
import java.io.PrintStream;

import au.edu.rmit.trajectory.clustering.kpaths.KPathsOptimization;
import au.edu.rmit.trajectory.clustering.streaming.Streaming;

public class StreamEvaluation {
    public static void main(String[] args) throws IOException, InterruptedException {
        runOnStreamingDataset(args);
        //runOnStaticDataset();
    }

    /**
     * test clustering over streaming data
     */
    private static void runOnStreamingDataset(String[] args) throws IOException, InterruptedException {
        KPathsOptimization<Object> kpaths = new KPathsOptimization<>(args);
        //window: 15min, 30min, 60min, 120min, 180min
        int[] windowSize = {900, 1800, 3600, 7200, 10800, 14600};
        int counter = 0;
        for (int k = 0; k < 5; k++) {
            Streaming stream = new Streaming(windowSize[k], windowSize[k], 1, args[6]);//args[6] is the streaming record.
            stream.start();
            for (int i = 0; i < 60; i++) {//the time stamp
                Thread.sleep(1000);//sleep one second
                stream.updateIndexes();
                String filename = "./logs/streaming/"+args[5]+"_"+windowSize[k]+"_"+args[1]+"_"+counter;
        		PrintStream fileOut = new PrintStream(filename);
        		System.setOut(fileOut);
                kpaths.streamingClustering(stream.trajectories, stream.edgeInvertedIndex, args, counter++);
            }
            System.err.println();
            stream.close();// this method must be called to shutdown the other thread which keeps reading from record file,
            // and this will release the record file.
        }
    }
}
