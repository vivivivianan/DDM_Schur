// Paste this source into the COMSOL 6.1 Application Builder > Methods editor.
// It runs inside an MPH model method: `model` is supplied by COMSOL.
// No ModelUtil calls and no external Java class are used.

// Make batch replays deterministic when the saved template already contains
// the model built during the Desktop smoke run.
model.component().clear();
model.component().create("comp1", true);
model.component("comp1").geom().create("geom1", 3);
model.component("comp1").geom("geom1").create("blk1", "Block");
model.component("comp1").geom("geom1").feature("blk1").set("size", new String[]{"1", "1", "1"});
model.component("comp1").geom("geom1").run();

model.component("comp1").mesh().create("mesh1");
model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
model.component("comp1").mesh("mesh1").run();
model.component("comp1").mesh("mesh1").export(
  "D:\\AI agent\\codex\\DDM_Schur\\data\\generated\\two_by_three_comsol61\\smoke\\method_smoke.mphtxt"
);
