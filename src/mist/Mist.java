package mist;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;

import java.util.ResourceBundle;

public class Mist extends Application {

    private  static Stage stage;

    @Override
    public void start(Stage primaryStage) throws Exception{
        final ResourceBundle resourceBundle = Texts.getResources();
        Parent root = FXMLLoader.load(getClass().getResource("view.fxml"), resourceBundle);
        primaryStage.setTitle("Mist");
        primaryStage.setScene(new Scene(root));
        Mist.stage = primaryStage;
        primaryStage.show();
    }


    public static void main(String[] args) {
        launch(args);
    }

    public static Stage getStage(){
        return stage;

    }
}
