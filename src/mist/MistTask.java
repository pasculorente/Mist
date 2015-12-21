package mist;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleObjectProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.concurrent.Task;

import java.io.*;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

/**
 * MistTask allows to perform a Mist analysis over a BAM file. Given an input BAM file, a threshold and a length,
 * MistTask will create a .mist output file with all the mist regions.
 * <p>
 * Given the following exon:
 * <p>
 * Pos DP<br>
 * 1000 5<br>
 * 1001 4<br>
 * 1002 8<br>
 * 1003 9<br>
 * 1004 10<br>
 * 1005 11<br>
 * 1006 9<br>
 * <p>
 * With threshold=10 and length=1 the mist regions would be [1000-1003] and [1006,1006].
 * With threshold=10 and length=3 the mist regions would be [1000-1003].
 * With threshold=9 and length=2 the mist regions would be [1000-1002].
 *
 * @author Lorente Arencibia, Pascual (pasculorente@gmail.com)
 */
public class MistTask extends Task<Integer> {


    private final File input, output;
    private final int threshold, length;

    private final static int WINDOW_SIZE = 10;
    private final static String INSIDE = "inside";
    private final static String OVERLAP = "overlap";
    private final static String LEFT = "left";
    private final static String RIGHT = "right";

    // chrom | start | end | gene_id | gene_name | exon_number | transcript_id | transcript_name |
    // transcript_info | gene_biotype
    private final static int EXON_CHR = 0;
    private final static int EXON_START = 1;
    private final static int EXON_END = 2;
    private final static int GENE_ID = 3;
    private final static int GENE_NAME = 4;
    private final static int EXON_N = 5;
    private final static int EXON_ID = 6;
    private final static int TRANS_NAME = 7;
    private final static int TRANS_INFO = 8;
    private final static int GENE_BIO = 9;

    private long genomeLength;

    private long startTime;
    private List<Chromosome> chromosomes;

    final AtomicReference<String> currentChromosome = new AtomicReference<>("0");
    private int[] depths;

//    AtomicInteger matches = new AtomicInteger();

    private Property<Integer> matches = new SimpleObjectProperty<>();
    private Property<String> coordinate = new SimpleStringProperty();
    private Property<Long> elapsed = new SimpleObjectProperty<>();
    private Property<Long> remaining = new SimpleObjectProperty<>();

    /**
     * Parameters are not checked inside MIST, please, be sure all of them are legal.
     *
     * @param input     the input BAM
     * @param output    the output MIST
     * @param threshold the DP threshold
     * @param length    the minimum length
     */
    public MistTask(File input, File output, int threshold, int length) {
        this.input = input;
        this.output = output;
        this.threshold = threshold;
        this.length = length;
    }

    final String[] headers = {"chrom", "exon_start", "exon_end", "mist_start", "mist_end",
            "gene_id", "gene_name", "exon_number", "exon_id", "transcript_name", "biotype", "match"};


    /*
     * IMPORTANT NOTE FOR DEVELOPERS. Genomic positions start at 1, Java array positions start at 0.
     * To avoid confusions, all Java arrays will have length incremented in 1, and I won't use
     * position 0 in them. So any time there is an array access (depths[i]) it is accessing
     * to genomic position.
     * NOTE 2: Firs implementation had a high cost: read each exon from Ensembl and call 'samtools
     * mpileup' for each. Even with parallelization, its estimated time was 2 or 3 days for a sample.
     * Current implementation piles up chromosome by chromosome by requesting.
     */
    @Override
    protected Integer call() throws Exception {
        int ret = startMIST();
        updateProgress(1, 1);
        updateMessage(ret == 0
                ? Texts.getString("finished.successfully") + ": " + Texts.getFormattedString("regions.found", getMatches())
                : Texts.getString("finished.with.errors") + ": " + Texts.getFormattedString("regions.found", getMatches()));
        return ret;

    }

    /*
     * 1: write headers
     * 2: Read exons
     * 3: if exon's chr is not loaded, mpileup in memory chr
     * 4: Locate mist regions.
     * 5: save mist regions
     */
    private int startMIST() {
        matches.setValue(0);
        updateProgress(0, 1);
        readChromosomes();
        writeHeader(output);
        startTime = System.currentTimeMillis();
        // Read the exons file
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(Mist.class.getResourceAsStream("HomoSapiens_v75_protein_coding.tsv.gz"))))) {
            // Skip first line
            reader.readLine();
            String line;
            while ((line = reader.readLine()) != null && !isCancelled()) {
                String[] row = line.split("\t");
                final String chr = row[0];
                if (!currentChromosome.get().equals(chr)) nextChromosome(chr);
                if (depths != null) computeMistRegions(row);
            }
        } catch (Exception e) {
            // Dont do anything, someone canceled the stream
            return 1;
        }
        return 0;
    }

    private void readChromosomes() {
        chromosomes = readBamHeaders(input);
        chromosomes.forEach(chromosome -> genomeLength += chromosome.length);
    }

    private void nextChromosome(String chr) {
        depths = readBamContent(chr);
        chromosomes.stream().
                filter(chromosome -> chromosome.name.equals(chr)).
                forEach(chromosome -> chromosome.processed = true);
        currentChromosome.set(chr);
    }

    /**
     * Writes all the args in the file using a tab as separator and insert a newLine mark at the
     * end.
     *
     * @param output output file
     * @param values list of values
     */
    private void writeLine(File output, String... values) {
        try (BufferedWriter out = new BufferedWriter(new FileWriter(output, true))) {
            out.write(asString("\t", values));
            out.newLine();
        } catch (IOException ex) {
            Logger.getLogger(Mist.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Writes the first line of the output line, which contains the headers for the columns.
     *
     * @param output outputFile
     */
    private void writeHeader(File output) {
        if (output.exists()) output.delete();
        writeLine(output, headers);
    }

    private void computeMistRegions(String[] exon) {
        int start = getStart(exon[EXON_START]);
        int end = getEnd(exon[EXON_END]);
        findMistRegions(exon, start, end);

    }

    private int getStart(String s) {
        int start = Integer.parseInt(s) - WINDOW_SIZE;
        return start < 1 ? 1 : start;
    }

    private int getEnd(String s) {
        int end = Integer.parseInt(s) + WINDOW_SIZE;
        return end <= depths.length ? end : depths.length;
    }

    private void findMistRegions(String[] exon, int start, int end) {
        boolean inMist = false;
        int mistStart = 0;
        for (int i = start; i < end; i++) {
            if (depths[i] < threshold) {
                if (!inMist) {
                    inMist = true;
                    mistStart = i;
                }
            } else {
                if (inMist) {
                    inMist = false;
                    if (checkRegionLength(mistStart, i - 1)) {
                        printMist(exon, mistStart, i - 1);
                        matches.setValue(getMatches() + 1);
                    }
                }
            }
        }
    }

    /**
     * Stores a MIST region only if its length is greater than the length parameter.
     *
     * @param line      the TSV exon
     * @param mistStart the start of the mist region
     * @param mistEnd   the end of the mist region
     */
    private void printMist(String[] line, int mistStart, int mistEnd) {
        final int exonStart = Integer.valueOf(line[EXON_START]);
        final int exonEnd = Integer.valueOf(line[EXON_END]);
        // Determine type of match
        String match = determineMatch(exonStart, exonEnd, mistStart, mistEnd);
        // chrom, exon_start, exon_end, mist_start, mist_end, gene_id, gene_name, exon_id,
        // transcript_name, biotype, match
        writeLine(output, line[EXON_CHR], line[EXON_START], line[EXON_END], mistStart + "", mistEnd + "", line[GENE_ID],
                line[GENE_NAME], line[EXON_N], line[EXON_ID], line[TRANS_NAME], line[GENE_BIO], match);

    }

    /**
     * True if the region is longer than length
     *
     * @param mistStart start of mist region
     * @param mistEnd   end of mist region
     * @return mistEnd - mistStart + 1 >= length
     */
    private boolean checkRegionLength(int mistStart, int mistEnd) {
        return mistEnd - mistStart + 1 >= length;
    }

    /**
     * Given an exon coordinates and a MIST region coordinates determines if the MIST region if
     * left, right, inside or overlapping the exon.
     *
     * @param exonStart start of the exon
     * @param exonEnd   end of the exon
     * @param mistStart start of the mist region
     * @param mistEnd   end of the mist region
     * @return left, rigth, inside or overlap
     */
    private String determineMatch(int exonStart, int exonEnd, int mistStart, int mistEnd) {
        return mistStart < exonStart
                ? mistEnd > exonEnd ? OVERLAP : LEFT
                : mistEnd > exonEnd ? RIGHT : INSIDE;
    }

    private void setProgress(String chr, int pos) {
        final long gpos = chromosomes.stream().
                filter(chromosome -> chromosome.processed).
                mapToLong(chromosome -> chromosome.length).
                sum() + pos;
        final double percentage = gpos * 100.0 / genomeLength;
        final long time = System.currentTimeMillis() - startTime;
        final long remaining = genomeLength * time / gpos - time;
//        final String elapsed = humanReadableTime(time);
//        final String rem = humanReadableTime(remaining);
//        updateMessage(String.format("%s (%s:%,d) %d matches (%s)", elapsed, chr, pos, getMatches(), rem));
        updateProgress(percentage, 100.0);
        this.coordinate.setValue(String.format("%s:%,d", chr, pos));
        this.elapsed.setValue(time);
        this.remaining.setValue(remaining);
    }

    public String humanReadableTime(long millis) {
        long days = TimeUnit.MILLISECONDS.toDays(millis);
        millis -= TimeUnit.DAYS.toMillis(days);
        long hours = TimeUnit.MILLISECONDS.toHours(millis);
        millis -= TimeUnit.HOURS.toMillis(hours);
        long minutes = TimeUnit.MILLISECONDS.toMinutes(millis);
        millis -= TimeUnit.MINUTES.toMillis(minutes);
        long seconds = TimeUnit.MILLISECONDS.toSeconds(millis);
        String ret = days > 0 ? days + " d " : "";
        return ret + String.format("%02d:%02d:%02d", hours, minutes, seconds);
    }

    /**
     * Creates a int[] with the length of chr + 1, using 1-based coordinates, son ret[0] is empty.
     *
     * @param chr the name of the chromosome
     * @return the array of the depths of the chromosome
     */
    private int[] readBamContent(String chr) {
        Chromosome next = getChromosome(chr);
        if (next == null) {
            System.out.println(Texts.getFormattedString("missed.chromosome", chr));
            return null;
        }
        return loadChromosome(next);
    }

    private Chromosome getChromosome(String chr) {
        Optional<Chromosome> first = chromosomes.stream()
                .filter(chromosome -> chromosome.name.equals(chr))
                .findFirst();
        return first.isPresent() ? first.get() : null;
    }

    private int[] loadChromosome(Chromosome chromosome) {
        final int le = (int) chromosome.length;
        final int[] depths = new int[le + 1];
        final ProcessBuilder pb = new ProcessBuilder("samtools", "mpileup", "-r", chromosome.name, input.getAbsolutePath());
        final AtomicInteger iterations = new AtomicInteger();
        try {
            Process process = pb.start();
            try (BufferedReader command = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                command.lines().parallel().filter(s -> !isCancelled()).map(line -> line.split("\t")).forEach(pileup -> {
                    final int pos = Integer.valueOf(pileup[1]);
                    final int depth = Integer.valueOf(pileup[3]);
                    depths[pos] = depth;
                    if (iterations.incrementAndGet() % 1000000 == 0) setProgress(chromosome.name, pos);
                });
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return depths;
    }

    /**
     * Reads input bam headers and returns a List with chromosome names and lengths or an empty list when no chromosomes
     * are found.
     *
     * @param input the sam or bam file.
     * @return a list of chromosomes in the BAM file. Collections.emptyList() if no chromosomes in the file
     */
    private List<Chromosome> readBamHeaders(File input) {
        // samtools view -H input.bam
        // @SQ	SN:1	LN:249250621
        // @SQ	SN:GL000249.1	LN:38502
        final ProcessBuilder pb = new ProcessBuilder("samtools", "view", "-H", input.getAbsolutePath());
        try {
            return mapToChromosomes(pb.start().getInputStream());
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return Collections.emptyList();
    }

    private List<Chromosome> mapToChromosomes(InputStream inputStream) {
        try (BufferedReader in = new BufferedReader(new InputStreamReader(inputStream))) {
            return in.lines()
                    .filter(line -> line.startsWith("@SQ"))
                    .map(line -> line.split("\t"))
                    .map(this::createChromosome)
                    .collect(Collectors.toList());
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return Collections.emptyList();
    }

    private Chromosome createChromosome(String[] line) {
        final String chr = line[1].substring(3);
        final int size = Integer.valueOf(line[2].substring(3));
        return new Chromosome(chr, size);
    }

    /**
     * Converts an Array to String using the separator. Omits the last separator. [value1 value2
     * value3] to value1,value2,value3
     *
     * @param separator something like "\t" or ","
     * @param values    a list of values
     * @return the stringified list
     */
    private static String asString(String separator, String[] values) {
        if (values.length == 0) return "";
        final StringBuilder builder = new StringBuilder(values[0]);
        for (int i = 1; i < values.length; i++) builder.append(separator).append(values[i]);
        return builder.toString();
    }

    @Override
    protected void cancelled() {
        super.cancelled();
        System.out.println("Cancelled");
    }

    public Property<Integer> matchesProperty() {
        return matches;
    }

    public int getMatches() {
        return matches.getValue();
    }

    public Property<Long> elapsedProperty() {
        return elapsed;
    }

    public Property<Long> remainingProperty() {
        return remaining;
    }

    public Property<String> coordinateProperty() {
        return coordinate;
    }

    /**
     * Tiny class to store together a chrom with its length and an already processed flag. This is
     * only used for progress purpose.
     */
    private class Chromosome {
        String name;
        long length;

        boolean processed = false;

        public Chromosome(String name, long length) {
            this.name = name;
            this.length = length;
        }

    }
}
