import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

/** COMSOL 5.6 batch model builder. Entity IDs are queried after geom.run(). */
public class BuildTwoByThreeCubes {
  static final String OUT="D:\\AI agent\\codex\\DDM_Schur\\data\\generated\\two_by_three_comsol";
  static String join(int[] a) { StringBuilder s=new StringBuilder(); for(int x:a)s.append(x).append(';'); return s.toString(); }
  static int[] toArray(ArrayList<Integer> a) { int[] r=new int[a.size()]; for(int i=0;i<r.length;i++)r[i]=a.get(i).intValue(); return r; }
  static int[] entities(Model m,String tag) { return m.component("comp1").selection(tag).entities(); }
  static void setBox(Model m,String tag,int dim,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
    m.component("comp1").geom("geom1").create(tag,"BoxSelection");
    m.component("comp1").geom("geom1").feature(tag).set("entitydim",dim);
    m.component("comp1").geom("geom1").feature(tag).set("xmin",Double.toString(xmin)); m.component("comp1").geom("geom1").feature(tag).set("xmax",Double.toString(xmax));
    m.component("comp1").geom("geom1").feature(tag).set("ymin",Double.toString(ymin)); m.component("comp1").geom("geom1").feature(tag).set("ymax",Double.toString(ymax));
    m.component("comp1").geom("geom1").feature(tag).set("zmin",Double.toString(zmin)); m.component("comp1").geom("geom1").feature(tag).set("zmax",Double.toString(zmax));
    m.component("comp1").geom("geom1").feature(tag).set("condition","inside");
  }
  static void material(Model m,int i,int did,double k) {
    String t="mat"+i; m.component("comp1").material().create(t,"Common"); m.component("comp1").material(t).selection().set(new int[]{did});
    m.component("comp1").material(t).propertyGroup("def").set("thermalconductivity",new String[]{k+"[W/(m*K)]","0","0","0",k+"[W/(m*K)]","0","0","0",k+"[W/(m*K)]"});
    m.component("comp1").material(t).propertyGroup("def").set("density","1000[kg/m^3]"); m.component("comp1").material(t).propertyGroup("def").set("heatcapacity","1000[J/(kg*K)]");
  }
  public static void main(String[] args) {
    try {
      ModelUtil.clear(); Model m=ModelUtil.create("two_by_three_cubes"); m.modelPath(OUT);
      m.component().create("comp1",true); m.component("comp1").geom().create("geom1",3); m.component("comp1").geom("geom1").lengthUnit("mm");
      String[] n={"C1","C2","C3","C4","C5","C6"}; String[] t={"c1","c2","c3","c4","c5","c6"}; int[][] xy={{0,0},{10,0},{20,0},{0,10},{10,10},{20,10}};
      for(int i=0;i<6;i++){ m.component("comp1").geom("geom1").create(t[i],"Block"); m.component("comp1").geom("geom1").feature(t[i]).set("size",new String[]{"10","10","10"}); m.component("comp1").geom("geom1").feature(t[i]).set("pos",new String[]{""+xy[i][0],""+xy[i][1],"0"}); m.component("comp1").geom("geom1").feature(t[i]).set("selresult","on"); }
      m.component("comp1").geom("geom1").create("uni1","Union"); m.component("comp1").geom("geom1").feature("uni1").selection("input").set(t); m.component("comp1").geom("geom1").feature("uni1").set("intbnd",true);
      // Thin geometric boxes return actual boundary IDs after the union is finalized.
      setBox(m,"bleft",2,-.001,.001,-.001,20.001,-.001,10.001); setBox(m,"bright",2,29.999,30.001,-.001,20.001,-.001,10.001);
      setBox(m,"btop",2,-.001,30.001,19.999,20.001,-.001,10.001); setBox(m,"bbot",2,-.001,30.001,-.001,.001,-.001,10.001);
      setBox(m,"bfront",2,-.001,30.001,-.001,20.001,-.001,.001); setBox(m,"bback",2,-.001,30.001,-.001,20.001,9.999,10.001);
      setBox(m,"bi12",2,9.999,10.001,-.001,10.001,-.001,10.001); setBox(m,"bi13",2,9.999,10.001,9.999,20.001,-.001,10.001); setBox(m,"bi23",2,9.999,30.001,9.999,10.001,-.001,10.001);
      m.component("comp1").geom("geom1").run();
      int[] did=new int[6]; for(int i=0;i<6;i++) did[i]=entities(m,"geom1_"+t[i]+"_dom")[0];
      int[] left=entities(m,"geom1_bleft_bnd"), right=entities(m,"geom1_bright_bnd"), top=entities(m,"geom1_btop_bnd"), bot=entities(m,"geom1_bbot_bnd"), front=entities(m,"geom1_bfront_bnd"), back=entities(m,"geom1_bback_bnd");
      int[] I12=entities(m,"geom1_bi12_bnd"); int[] I13=entities(m,"geom1_bi13_bnd"); int[] I23=entities(m,"geom1_bi23_bnd");
      double[] k={10,1,20,2,15,.5}; for(int i=0;i<6;i++)material(m,i+1,did[i],k[i]);
      m.component("comp1").physics().create("ht","HeatTransferInSolids","geom1"); m.component("comp1").physics("ht").feature("init1").set("Tinit","300[K]");
      m.component("comp1").physics("ht").create("temp_left","TemperatureBoundary",2); m.component("comp1").physics("ht").feature("temp_left").selection().set(left); m.component("comp1").physics("ht").feature("temp_left").set("T0","300[K]");
      m.component("comp1").physics("ht").create("temp_right","TemperatureBoundary",2); m.component("comp1").physics("ht").feature("temp_right").selection().set(right); m.component("comp1").physics("ht").feature("temp_right").set("T0","350[K]");
      m.component("comp1").physics("ht").create("conv_top","HeatFluxBoundary",2); m.component("comp1").physics("ht").feature("conv_top").selection().set(top); m.component("comp1").physics("ht").feature("conv_top").set("HeatFluxType","ConvectiveHeatFlux"); m.component("comp1").physics("ht").feature("conv_top").set("h","10[W/(m^2*K)]"); m.component("comp1").physics("ht").feature("conv_top").set("Text","293.15[K]");
      m.component("comp1").physics("ht").create("flux_bottom","HeatFluxBoundary",2); m.component("comp1").physics("ht").feature("flux_bottom").selection().set(bot); m.component("comp1").physics("ht").feature("flux_bottom").set("HeatFluxType","GeneralInwardHeatFlux"); m.component("comp1").physics("ht").feature("flux_bottom").set("q0","1000[W/m^2]");
      m.component("comp1").physics("ht").create("ins_z","ThermalInsulation",2); int[] z=new int[front.length+back.length]; System.arraycopy(front,0,z,0,front.length);System.arraycopy(back,0,z,front.length,back.length);m.component("comp1").physics("ht").feature("ins_z").selection().set(z);
      m.component("comp1").physics("ht").create("src_c5","HeatSource",3);m.component("comp1").physics("ht").feature("src_c5").selection().set(new int[]{did[4]});m.component("comp1").physics("ht").feature("src_c5").set("Q0","1e6[W/m^3]");
      m.component("comp1").mesh().create("mesh_ref");m.component("comp1").mesh("mesh_ref").feature("size").set("hauto",2);m.component("comp1").mesh("mesh_ref").run(); m.study().create("std_stationary");m.study("std_stationary").create("stat","Stationary");m.study("std_stationary").run();m.study().create("std_time");m.study("std_time").create("time","Transient");m.study("std_time").feature("time").set("tlist","range(0,1,100)");m.study("std_time").run();m.component("comp1").mesh("mesh_ref").export(OUT+"\\reference_full.mphtxt");
      int[][] sd={{did[0],did[3]},{did[1],did[2]},{did[4],did[5]}};double[] hm={2.5,1.8,2.2}; for(int s=0;s<3;s++){String mt="mesh_sd"+(s+1);m.component("comp1").mesh().create(mt);m.component("comp1").mesh(mt).create("ftet1","FreeTet");m.component("comp1").mesh(mt).feature("ftet1").selection().geom("geom1",3);m.component("comp1").mesh(mt).feature("ftet1").selection().set(sd[s]);m.component("comp1").mesh(mt).create("size1","Size");m.component("comp1").mesh(mt).feature("size1").selection().geom("geom1",3);m.component("comp1").mesh(mt).feature("size1").selection().set(sd[s]);m.component("comp1").mesh(mt).feature("size1").set("hmax",hm[s]+"[mm]");m.component("comp1").mesh(mt).feature("size1").set("hmin",(hm[s]/4)+"[mm]");m.component("comp1").mesh(mt).run();m.component("comp1").mesh(mt).export(OUT+"\\subdomain_"+(s+1)+".mphtxt");}
      m.save(OUT+"\\two_by_three_cubes.mph"); System.out.println("COMSOL_BUILD_COMPLETE domains="+join(did)+" I12="+join(I12)+" I13="+join(I13)+" I23="+join(I23));
    } catch(Exception e) { e.printStackTrace(); throw new RuntimeException(e); }
  }
}
