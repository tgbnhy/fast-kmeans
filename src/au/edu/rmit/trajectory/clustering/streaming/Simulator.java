package au.edu.rmit.trajectory.clustering.streaming;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

class Simulator {

    ConcurrentLinkedDeque<Tuple> stream;
    private BufferedReader reader;
    private int window;
    private AtomicInteger cur;
    private int startFrom;
    private Tuple preT;
    private int speedupFactor;
    private ExecutorService worker;
    private boolean stop = false;

    /**
     * (startFrom+1) should be dividable by speedupFactor
     *
     * e.g.
     * 3600, 3600, 60: read in an hour data before testing, the window size is an hour. When passed 1 second in real world, 60 second has passed in simulator.
     * 7200, 7200, 60: read in two hours data before testing, the window size is two hour. When passed 1 second in real world, 60 second has passed in simulator
     *
     */
    Simulator(int startFrom, int window, int speedupFactor, String datafile) throws IOException {
        this.speedupFactor = speedupFactor;
        stream = new ConcurrentLinkedDeque<>();
        reader = new BufferedReader(new FileReader(datafile));
        this.window = window;
        cur = new AtomicInteger(0);

        //preload the list
        String line;
        while ((line = reader.readLine()) !=null){
            Tuple t = Tuple.genFromCSV(line);
            if (t.timestemp >= startFrom){
                preT = t;
                break;
            }else{
                stream.addLast(t);
            }
        }

        this.startFrom = startFrom / speedupFactor;
        cur.set(this.startFrom);
    }

    /**
     * Start the simulation from the parameter specified in constructor.
     * Suppose the parameters in constructor is 3600, 3600, 60,
     * now we start from 3600 second, we will read data from 3600 second to 3660 second in current second, which depends on speedupFactor
     */
    public void start(){
        worker = Executors.newSingleThreadExecutor();
        worker.execute(()->{
            long startTime = System.currentTimeMillis();

            while (!stop){

                while ((System.currentTimeMillis() - startTime) / 1000 < cur.intValue() - startFrom + 1){
                    try {
                        Thread.sleep(50);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }

                try {
                    stop = update();
                } catch (IOException e) {
                    e.printStackTrace();
                }

                cur.incrementAndGet();
            }
            worker.shutdown();
            try {
                reader.close();
                worker.awaitTermination(2, TimeUnit.SECONDS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        });
    }

    private boolean update() throws IOException {

        if (cur.intValue()*speedupFactor - window + 1 > 0) {
            int discardTo = cur.intValue()*speedupFactor - window;
            while (stream.getFirst()!=null && stream.getFirst().timestemp <= discardTo) {
                stream.removeFirst();
            }
        }

        if (preT != null)
            stream.addLast(preT);

        String line ;
        while (true) {
            if ((line = reader.readLine()) == null)
                return true;

            Tuple t = Tuple.genFromCSV(line);
            if (t.timestemp > cur.intValue() * speedupFactor) {
                preT = t;
                return false;
            } else {
                stream.addLast(t);
            }
        }


    }

    public void close() throws IOException {
        stop = true;
        worker.shutdown();
        try {
            worker.awaitTermination(2, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }
}
