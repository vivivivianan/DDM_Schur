// COMSOL 6.1 Model Method body. Paste into a Method named build_two_by_three.
// The supplied `model` object is used directly; no ModelUtil or external class.
String out = "D:\\AI agent\\codex\\DDM_Schur\\data\\generated\\two_by_three_comsol61";
model.component().clear();
model.component().create("comp1", true);
model.component("comp1").geom().create("geom1", 3);
model.component("comp1").geom("geom1").lengthUnit("mm");
String[] cubeTag={"c1","c2","c3","c4","c5","c6"};
String[] cubeName={"C1","C2","C3","C4","C5","C6"};
String[][] pos={{"0","0","0"},{"10","0","0"},{"20","0","0"},{"0","10","0"},{"10","10","0"},{"20","10","0"}};
for(int i=0;i<6;i++) {
  model.component("comp1").geom("geom1").create(cubeTag[i],"Block");
  model.component("comp1").geom("geom1").feature(cubeTag[i]).set("size",new String[]{"10","10","10"});
  model.component("comp1").geom("geom1").feature(cubeTag[i]).set("pos",pos[i]);
  model.component("comp1").geom("geom1").feature(cubeTag[i]).set("selresult","on");
}
model.component("comp1").geom("geom1").create("uni1","Union");
model.component("comp1").geom("geom1").feature("uni1").selection("input").set(cubeTag);
model.component("comp1").geom("geom1").feature("uni1").set("intbnd",true);
String[] btag={"bleft","bright","btop","bbot","bfront","bback","bi12","bi13","bi23"};
double[][] box={{-0.001,0.001,-0.001,20.001,-0.001,10.001},{29.999,30.001,-0.001,20.001,-0.001,10.001},{-0.001,30.001,19.999,20.001,-0.001,10.001},{-0.001,30.001,-0.001,0.001,-0.001,10.001},{-0.001,30.001,-0.001,20.001,-0.001,0.001},{-0.001,30.001,-0.001,20.001,9.999,10.001},{9.999,10.001,-0.001,10.001,-0.001,10.001},{9.999,10.001,9.999,20.001,-0.001,10.001},{9.999,30.001,9.999,10.001,-0.001,10.001}};
for(int i=0;i<btag.length;i++) {
  model.component("comp1").geom("geom1").create(btag[i],"BoxSelection");
  model.component("comp1").geom("geom1").feature(btag[i]).set("entitydim",2);
  model.component("comp1").geom("geom1").feature(btag[i]).set("xmin",Double.toString(box[i][0])); model.component("comp1").geom("geom1").feature(btag[i]).set("xmax",Double.toString(box[i][1]));
  model.component("comp1").geom("geom1").feature(btag[i]).set("ymin",Double.toString(box[i][2])); model.component("comp1").geom("geom1").feature(btag[i]).set("ymax",Double.toString(box[i][3]));
  model.component("comp1").geom("geom1").feature(btag[i]).set("zmin",Double.toString(box[i][4])); model.component("comp1").geom("geom1").feature(btag[i]).set("zmax",Double.toString(box[i][5]));
  model.component("comp1").geom("geom1").feature(btag[i]).set("condition","inside");
}
model.component("comp1").geom("geom1").run();
int[] did=new int[6];
for(int i=0;i<6;i++) did[i]=model.component("comp1").selection("geom1_"+cubeTag[i]+"_dom").entities()[0];
int[][] bid=new int[9][];
for(int i=0;i<btag.length;i++) bid[i]=model.component("comp1").selection("geom1_"+btag[i]+"_bnd").entities();
double[] kval={10,1,20,2,15,0.5};
for(int i=0;i<6;i++) {
  String mt="mat"+(i+1); model.component("comp1").material().create(mt,"Common"); model.component("comp1").material(mt).selection().set(new int[]{did[i]});
  String k=Double.toString(kval[i])+"[W/(m*K)]";
  model.component("comp1").material(mt).propertyGroup("def").set("thermalconductivity",new String[]{k,"0","0","0",k,"0","0","0",k});
  model.component("comp1").material(mt).propertyGroup("def").set("density","1000[kg/m^3]"); model.component("comp1").material(mt).propertyGroup("def").set("heatcapacity","1000[J/(kg*K)]");
}
model.component("comp1").physics().create("ht","HeatTransferInSolids","geom1");
model.component("comp1").physics("ht").feature("init1").set("Tinit","300[K]");
model.component("comp1").physics("ht").create("temp_left","TemperatureBoundary",2); model.component("comp1").physics("ht").feature("temp_left").selection().set(bid[0]); model.component("comp1").physics("ht").feature("temp_left").set("T0","300[K]");
model.component("comp1").physics("ht").create("temp_right","TemperatureBoundary",2); model.component("comp1").physics("ht").feature("temp_right").selection().set(bid[1]); model.component("comp1").physics("ht").feature("temp_right").set("T0","350[K]");
model.component("comp1").physics("ht").create("conv_top","HeatFluxBoundary",2); model.component("comp1").physics("ht").feature("conv_top").selection().set(bid[2]); model.component("comp1").physics("ht").feature("conv_top").set("HeatFluxType","ConvectiveHeatFlux"); model.component("comp1").physics("ht").feature("conv_top").set("h","10[W/(m^2*K)]"); model.component("comp1").physics("ht").feature("conv_top").set("Text","293.15[K]");
model.component("comp1").physics("ht").create("flux_bottom","HeatFluxBoundary",2); model.component("comp1").physics("ht").feature("flux_bottom").selection().set(bid[3]); model.component("comp1").physics("ht").feature("flux_bottom").set("HeatFluxType","GeneralInwardHeatFlux"); model.component("comp1").physics("ht").feature("flux_bottom").set("q0","1000[W/m^2]");
int[] z=new int[bid[4].length+bid[5].length]; System.arraycopy(bid[4],0,z,0,bid[4].length); System.arraycopy(bid[5],0,z,bid[4].length,bid[5].length);
model.component("comp1").physics("ht").create("ins_z","ThermalInsulation",2); model.component("comp1").physics("ht").feature("ins_z").selection().set(z);
model.component("comp1").physics("ht").create("src_c5","HeatSource",3); model.component("comp1").physics("ht").feature("src_c5").selection().set(new int[]{did[4]}); model.component("comp1").physics("ht").feature("src_c5").set("Q0","1e6[W/m^3]");
model.component("comp1").mesh().create("mesh_ref"); model.component("comp1").mesh("mesh_ref").feature("size").set("hauto",2); model.component("comp1").mesh("mesh_ref").run();
model.study().create("std_stationary"); model.study("std_stationary").create("stat","Stationary"); model.study("std_stationary").run();
model.study().create("std_time"); model.study("std_time").create("time","Transient"); model.study("std_time").feature("time").set("tlist","range(0,1,100)"); model.study("std_time").run();
model.component("comp1").mesh("mesh_ref").export(out+"\\reference_full.mphtxt");
int[][] sd={{did[0],did[3]},{did[1],did[2]},{did[4],did[5]}}; double[] hmax={2.5,1.8,2.2};
for(int s=0;s<3;s++) { String mt="mesh_sd"+(s+1); model.component("comp1").mesh().create(mt); model.component("comp1").mesh(mt).create("ftet1","FreeTet"); model.component("comp1").mesh(mt).feature("ftet1").selection().geom("geom1",3); model.component("comp1").mesh(mt).feature("ftet1").selection().set(sd[s]); model.component("comp1").mesh(mt).create("size1","Size"); model.component("comp1").mesh(mt).feature("size1").selection().geom("geom1",3); model.component("comp1").mesh(mt).feature("size1").selection().set(sd[s]); model.component("comp1").mesh(mt).feature("size1").set("hmax",Double.toString(hmax[s])+"[mm]"); model.component("comp1").mesh(mt).feature("size1").set("hmin",Double.toString(hmax[s]/4)+"[mm]"); model.component("comp1").mesh(mt).run(); model.component("comp1").mesh(mt).export(out+"\\subdomain_"+(s+1)+".mphtxt"); }
java.io.PrintWriter csv=new java.io.PrintWriter(new java.io.File(out+"\\comsol_entity_map.csv"));
csv.println("kind,name,comsol_ids,subdomain,material_or_role");
for(int i=0;i<6;i++) csv.println("domain,"+cubeName[i]+","+did[i]+","+((i==0||i==3)?1:((i==1||i==2)?2:3))+",k="+kval[i]);
String[] role={"left_temperature","right_temperature","top_convection","bottom_inward_flux","front_insulation","back_insulation","interface_1_2","interface_1_3","interface_2_3"};
for(int i=0;i<9;i++) { String v=""; for(int j=0;j<bid[i].length;j++) v+=((j==0?"":";")+bid[i][j]); csv.println("boundary,"+btag[i]+","+v+",,"+role[i]); }
csv.close();
System.out.println("COMSOL_2X3_METHOD_OK domains="+java.util.Arrays.toString(did));
