import com.comsol.model.*;
import com.comsol.model.util.*;
public class ComsolBatchSmoke {
  private static final String OUTPUT_DIR = "D:\\AI agent\\codex\\DDM_Schur\\data\\generated\\two_by_three_comsol61\\smoke";
  public static Model run() {
    Model model=ModelUtil.create("smoke");
    model.component().create("comp1",true);
    model.component("comp1").geom().create("geom1",3);
    model.component("comp1").geom("geom1").create("blk1","Block");
    model.component("comp1").geom("geom1").feature("blk1").set("size",new String[]{"1","1","1"});
    model.component("comp1").geom("geom1").run();
    return model;
  }
  public static void main(String[] args) throws Exception {
    ModelUtil.loadPreferences();
    System.out.println("Preferences loaded");
    System.out.println("File system access key = security.external.filepermission");
    System.out.println("File system access = " + ModelUtil.getPreference("security.external.filepermission"));
    System.out.println("MPHTXT export path = " + OUTPUT_DIR + "\\smoke_cube.mphtxt");
    Model model=run();
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").run();
    model.component("comp1").mesh("mesh1").export(OUTPUT_DIR+"\\smoke_cube.mphtxt");
    System.out.println("COMSOL_SMOKE_OK");
  }
}
