package mist;

import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.TextField;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;

import java.io.File;
import java.util.concurrent.TimeUnit;

public class ViewController {
    public Label positionLabel;
    public Label elapsedTimeLabel;
    public Label remainingTimeLabel;
    public Label matchesLabel;
    public VBox progressInfo;
    @FXML
    private TextField fileTextField;
    @FXML
    private Label message;
    @FXML
    private ProgressBar progressBar;
    @FXML
    private TextField thresholdTextField;
    @FXML
    private TextField lengthTextField;
    @FXML
    private Button startButton;
    @FXML
    private Button browseButton;

    private MistTask task = null;

    @FXML
    private void initialize() {
        browseButton.requestFocus();
    }

    public void selectFile(ActionEvent event) {
        final File file = openBam();
        if (file != null) {
            fileTextField.setText(file.getAbsolutePath());
            startButton.setDisable(false);
        }
    }

    private File openBam() {
        final FileChooser chooser = new FileChooser();
        chooser.getExtensionFilters().add(new FileChooser.ExtensionFilter("Binary sequence alignment map file (.bam)", "*.bam"));
        chooser.setTitle(Texts.getString("select.input"));
        if (fileTextField.getText() != null && !fileTextField.getText().isEmpty())
            chooser.setInitialDirectory(new File(fileTextField.getText()).getParentFile());
        else chooser.setInitialDirectory(new File(System.getProperty("user.home")));
        return chooser.showOpenDialog(Mist.getStage());
    }

    public void start(ActionEvent event) {
        int threshold, length;
        try {
            threshold = Integer.valueOf(thresholdTextField.getText());
        } catch (NumberFormatException ex) {
            message.setText(Texts.getString("invalid.threshold"));
            return;
        }
        try {
            length = Integer.valueOf(lengthTextField.getText());
        } catch (NumberFormatException ex) {
            message.setText(Texts.getString("invalid.length"));
            return;
        }
        File output = saveMist();
        if (output != null) {
            if (!output.getName().endsWith(".mist")) output = new File(output.getAbsolutePath() + ".mist");
            final File input = new File(fileTextField.getText());
            task = new MistTask(input, output, threshold, length);
            prepareGUI(task);
            task.setOnSucceeded(event1 -> restoreGUI());
            task.setOnCancelled(event1 -> restoreGUI());
            task.setOnFailed(event1 -> restoreGUI());
            final Thread th = new Thread(task);
            th.setDaemon(false);
            th.start();
        }
    }

    private void prepareGUI(MistTask mistTask) {
        progressBar.progressProperty().bind(mistTask.progressProperty());
        message.textProperty().bind(mistTask.messageProperty());
        progressInfo.setVisible(true);
        startButton.setDisable(true);
        browseButton.setDisable(true);
        task.matchesProperty().addListener((observable, oldValue, newValue) -> Platform.runLater(() -> matchesLabel.setText(String.format("%,d", newValue))));
        task.elapsedProperty().addListener((observable, oldValue, newValue) -> Platform.runLater(() -> elapsedTimeLabel.setText(humanReadableTime(newValue))));
        task.remainingProperty().addListener((observable, oldValue, newValue) -> Platform.runLater(() -> remainingTimeLabel.setText(humanReadableTime(newValue))));
        task.coordinateProperty().addListener((observable, oldValue, newValue) -> Platform.runLater(() -> positionLabel.setText(newValue)));
    }

    private void restoreGUI() {
        progressBar.progressProperty().unbind();
        message.textProperty().unbind();
        progressInfo.setVisible(false);
        startButton.setDisable(false);
        browseButton.setDisable(false);
    }

    private File saveMist() {
        final FileChooser chooser = new FileChooser();
        chooser.getExtensionFilters().add(new FileChooser.ExtensionFilter("Mist file (.mist)", "*.mist"));
        chooser.setTitle(Texts.getString("select.output"));
        if (fileTextField.getText() != null && !fileTextField.getText().isEmpty())
            chooser.setInitialDirectory(new File(fileTextField.getText()).getParentFile());
        else chooser.setInitialDirectory(new File(System.getProperty("user.home")));
        return chooser.showSaveDialog(Mist.getStage());
    }

    public void stopTask(ActionEvent event) {
        if (task != null) task.cancel(true);
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
}
