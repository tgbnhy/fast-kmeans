package au.edu.rmit.trajectory.clustering.validity;

/**
 * Created by Josem on 21/04/2017.
 */
public class Indice {

    private double resultado;
    private long time;

    public Indice(double resultado, long time) {
        this.resultado = resultado;
        this.time = time;
    }

    public double getResultado() {
        return resultado;
    }

    public void setResultado(double resultado) {
        this.resultado = resultado;
    }

    public long getTime() {
        return time;
    }

    public void setTime(long time) {
        this.time = time;
    }
}
