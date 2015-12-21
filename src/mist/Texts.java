package mist;

import java.text.MessageFormat;
import java.util.ResourceBundle;

/**
 * @author Lorente Arencibia, Pascual (pasculorente@gmail.com)
 */
public class Texts {

    private final static ResourceBundle resources = ResourceBundle.getBundle("mist.texts");

    public static String getString (String string) {
        return resources.getString(string);
    }

    public static String getFormattedString(String key, Object... args) {
        return new MessageFormat(resources.getString(key), resources.getLocale()).format(args);

    }

    public static ResourceBundle getResources() {
        return resources;
    }
}
