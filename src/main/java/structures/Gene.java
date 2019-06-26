package structures;

public class Gene {
    private final String name;

    public Gene(String name) {
        this.name = name;
    }

    @Override
    public String toString() {
        return  name;
    }
}
