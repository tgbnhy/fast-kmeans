package jeffgriffith.mnist;

import static java.lang.Math.min;
import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

import org.junit.Test;

import jeffgriffith.mnist.MnistReader;

public class TestMnistReader {

	@Test
	public void test() throws FileNotFoundException {
		String LABEL_FILE = "E:/dataset/kmeans/minist/t10k-labels.idx1-ubyte";
		String IMAGE_FILE = "E:/dataset/kmeans/minist/t10k-images.idx3-ubyte";
		
		IMAGE_FILE = "E:/dataset/kmeans/minist/train-images.idx3-ubyte";
		LABEL_FILE = "E:/dataset/kmeans/minist/train-labels.idx1-ubyte";		

		int[] labels = MnistReader.getLabels(LABEL_FILE);
		List<int[][]> images = MnistReader.getImages(IMAGE_FILE);
		
		assertEquals(labels.length, images.size());
		assertEquals(28, images.get(0).length);
		assertEquals(28, images.get(0)[0].length);

		for (int i = 0; i < min(10, labels.length); i++) {
			printf("================= LABEL %d\n", labels[i]);
			printf("%s", MnistReader.renderImage(images.get(i)));
		}
		//write into the files to store the vectors
		String LOG_DIR = "E:/dataset/kmeans/minist_60000_784";
		PrintStream fileOut = new PrintStream(LOG_DIR);
		System.setOut(fileOut);	
		for (int i = 0; i < 60000; i++) {
			for(int j=0; j<28; j++) {
				for(int z=0; z<28; z++) {
					System.out.print(images.get(i)[j][z]+" ");
				}
			}
			System.out.println();
		}
	}

	public static void printf(String format, Object... args) {
		System.out.printf(format, args);
	}
}
