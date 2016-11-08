import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * This class provides methods which are used to manipulate files.
 *
 * @author pablonsilva
 * @version 20151203
 */
public class FileUtils {
    /**
     * Saves a file in disk given the path and the content.
     *
     * @param content
     * @param path
     */
    public static void saveFile(String content, String path) {
        try {
            FileWriter fw = new FileWriter(path, true);
            fw.write(content);
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Method responsible to create a file in a given directory
     *
     * @param content
     * @param path
     * @param append
     */
    public static void saveFile(String content, String path, boolean append) {
        try {
            FileWriter fw = new FileWriter(path, append);
            fw.write(content);
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * This method creates a directory and return if it was created.
     *
     * @param path
     * @return if the directory was created
     */
    public static boolean createDir(String path) {
        File f = new File(path);

        if (!f.exists())
            return f.mkdir();
        return false;
    }

    /**
     * This method removes a file and return if it was removed.
     *
     * @param path
     * @return if the file was removed
     */
    public static boolean removeFile(String path) {
        File f = new File(path);

        if (f.exists())
            return f.delete();
        return false;
    }

    /**
     * This method verifies if a file exists in the system.
     *
     * @param path
     * @return if the file exists
     */
    public static boolean fileExist(String path)
    {
        File f = new File(path);
        boolean exists = false;
        try {
            exists = f.exists();
        }catch (Exception e)
        {

        }

        return  exists;
    }


/**
     * This method verifies if a directory exists in the system.
     *
     * @param path
     * @return if the directory exists
     */
    public static boolean directoryExists(String path)
    {
        File f = new File(path);
        boolean exists = false;
        try {
            exists = f.isDirectory();
        }catch (Exception e)
        {

        }

        return  exists;
    }
}
