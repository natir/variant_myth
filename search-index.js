var searchIndex = new Map(JSON.parse('[\
["variant_myth",{"t":"CCCCCCCCCHHFFPGPPGPPPNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNFNNNNNNNNNFNNNNNNNNNNNNNNNNNNNNNNNNNNPGPPPPPPIPPNNNNNNNNNNNNNNNFFIFNNNNNNNNNNNNNNONNNNNNNNNNNNNNNNNNNNNNONNNNNNNNNNNNFFGPPPPPPPGPPPPPPPPFPPPPPPPPPPPPPPPPPPNNNNNNNNNNNNNNNNNNNNONNNNNNNNNNNNNNNNNNNNNNNONNNNONNNNNNNNNNNNNNNFNNNNNNNHNNNFNNNNNNNNHHNNNNFFONNNNNNNNNNNNNNNNOOONNNNNNNN","n":["annotation","annotations_db","cli","error","interval_tree","myth","sequences_db","translate","variant","variants2myth","vcf2json","Annotation","Attribute","Forward","Frame","One","Reverse","Strand","Two","Unknow","Zero","borrow","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","borrow_mut","clone","clone","clone","clone","clone_into","clone_into","clone_into","clone_into","default","fmt","fmt","fmt","fmt","fmt","fmt","fmt","fmt","from","from","from","from","from_annotation","from_byte_record","from_u8_slice","get_attribute","get_exon_number","get_feature","get_interval","get_seqname","get_source","get_strand","get_transcript_id","get_transcript_id","into","into","into","into","to_owned","to_owned","to_owned","to_owned","to_string","to_string","to_string","to_string","try_from","try_from","try_from","try_from","try_into","try_into","try_into","try_into","type_id","type_id","type_id","type_id","AnnotationsDataBase","borrow","borrow_mut","from","from_reader","get_annotation","into","try_from","try_into","type_id","Command","annotations","augment_args","augment_args_for_update","borrow","borrow_mut","command","command_for_update","fmt","from","from_arg_matches","from_arg_matches_mut","group_id","into","output","quiet","reference","timestamp","translate","try_from","try_into","type_id","update_from_arg_matches","update_from_arg_matches_mut","updown_distance","variant","verbosity","Err","Error","FloatParsing","GffBadFrame","GffBadStrand","IntParsing","Log","Ok","Result","StdIO","VcfBadRecord","borrow","borrow_mut","fmt","fmt","from","from","from","from","from","into","source","to_string","try_from","try_into","type_id","Entry","InternalEntry","Interval","IntervalTree","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","clone","clone","clone","clone_into","clone_into","clone_into","data","default","end","eq","eq","eq","find","find_into","fmt","fmt","fmt","from","from","from","from_iter","hash","hash","hash","index","insert","interval","into","into","into","new","start","to_owned","to_owned","to_owned","try_from","try_from","try_from","try_into","try_into","try_into","type_id","type_id","type_id","AnnotationMyth","AnnotationMythBuilder","AnnotationMythBuilderError","Cds","CodonChange","CodonChangePlusCodonDeletion","CodonChangePlusCodonInsertion","CodonDeletion","CodonInsertion","Downstream","Effect","Exon","ExonDeleted","FrameShift","Gene","Intergenic","IntergenicConserved","Intron","IntronConserved","Myth","NonSynonymousCoding","SpliceSiteAcceptor","SplitceSiteDonor","StartGained","StartLost","StopGained","StopLost","SynonymousCoding","SynonymousStart","SynonymousStop","Transcript","UninitializedField","Upstream","Utr3Deleted","Utr3Prime","Utr5Deleted","Utr5Prime","ValidationError","add_annotation","add_effect","borrow","borrow","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","build","builder","clone","clone","clone_into","clone_into","default","effects","effects","extend_effect","fmt","fmt","fmt","fmt","fmt","from","from","from","from","from","from","from","from_variant","into","into","into","into","into","serialize","serialize","serialize","source","source","to_owned","to_owned","to_string","transcript_id","transcript_id","try_from","try_from","try_from","try_from","try_from","try_into","try_into","try_into","try_into","try_into","type_id","type_id","type_id","type_id","type_id","SequencesDataBase","borrow","borrow_mut","from","from_reader","get_interval","get_transcript","into","rev_comp","try_from","try_into","type_id","Translate","borrow","borrow_mut","from","from_reader","get_aa","into","is_start","is_stop","nuc2bit","seq2bit","translate","try_from","try_into","type_id","Variant","VcfReader","alt_seq","borrow","borrow","borrow_mut","borrow_mut","clone","clone_into","fmt","from","from","from_byte_record","from_reader","get_interval","into","into","into_iter","next","position","ref_seq","seqname","serialize","to_owned","try_from","try_from","try_into","try_into","type_id","type_id"],"q":[[0,"variant_myth"],[11,"variant_myth::annotation"],[86,"variant_myth::annotations_db"],[96,"variant_myth::cli"],[123,"variant_myth::error"],[149,"variant_myth::interval_tree"],[203,"variant_myth::myth"],[306,"variant_myth::sequences_db"],[318,"variant_myth::translate"],[333,"variant_myth::variant"],[363,"std::io"],[364,"core::fmt"],[365,"core::result"],[366,"csv::byte_record"],[367,"core::ops::range"],[368,"alloc::string"],[369,"core::any"],[370,"alloc::boxed"],[371,"std::io::buffered::bufreader"],[372,"alloc::vec"],[373,"clap_builder::builder::command"],[374,"clap_builder::parser::matches::arg_matches"],[375,"clap_builder"],[376,"clap_builder::util::id"],[377,"core::option"],[378,"std::io::buffered::bufwriter"],[379,"stderrlog"],[380,"log"],[381,"std::io::error"],[382,"core::num::error"],[383,"core::num::dec2flt"],[384,"core::error"],[385,"core::clone"],[386,"core::cmp"],[387,"core::marker"],[388,"core::convert"],[389,"core::iter::traits::collect"],[390,"core::hash"],[391,"derive_builder::error"],[392,"serde::ser"]],"i":[0,0,0,0,0,0,0,0,0,0,0,0,0,11,0,12,11,0,12,12,12,11,12,13,14,11,12,13,14,11,12,13,14,11,12,13,14,13,11,11,12,12,13,13,14,14,11,12,13,14,14,14,13,14,13,14,14,14,14,14,13,14,11,12,13,14,11,12,13,14,11,12,13,14,11,12,13,14,11,12,13,14,11,12,13,14,0,1,1,1,1,1,1,1,1,1,0,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,8,0,41,41,41,41,41,8,0,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,0,0,0,0,47,51,52,47,51,52,47,51,52,47,51,52,51,52,29,47,51,52,52,52,47,51,52,47,51,52,52,47,51,52,52,52,51,47,51,52,52,29,47,51,52,47,51,52,47,51,52,47,51,52,0,0,0,63,63,63,63,63,63,63,0,63,63,63,63,63,63,63,63,0,63,63,63,63,63,63,63,63,63,63,63,64,63,63,63,63,63,64,5,62,63,62,64,61,5,63,62,64,61,5,62,61,63,61,63,61,62,62,61,62,63,64,64,61,5,63,62,64,64,64,61,5,5,63,62,64,61,5,63,61,5,62,61,63,61,64,62,61,63,62,64,61,5,63,62,64,61,5,63,62,64,61,5,0,2,2,2,2,2,2,2,0,2,2,2,0,3,3,3,3,3,3,3,3,0,0,3,3,3,3,0,0,4,6,4,6,4,4,4,4,6,4,4,6,4,6,4,6,6,4,4,4,4,4,6,4,6,4,6,4],"f":"`````````{{bdfh}j}{{bdf{l{c}}e}{{A`{n}}}AbAd}``````````{ce{}{}}0000000{AfAf}{AhAh}{AjAj}{AlAl}{{ce}n{}{}}000{{}Aj}{{AfAn}{{Bb{nB`}}}}{{AfAn}Bd}{{AhAn}Bd}{{AhAn}{{Bb{nB`}}}}{{AjAn}Bd}{{AjAn}{{Bb{nB`}}}}{{AlAn}Bd}{{AlAn}{{Bb{nB`}}}}{cc{}}000{{Al{Bh{Bf}}}Al}{Bj{{A`{Al}}}}{{{Bh{Bf}}}{{A`{Aj}}}}{AlAj}{AjBl}{Al{{Bh{Bf}}}}{Al{{Bn{Bl}}}}11{AlAf}{Aj{{Bh{Bf}}}}3{ce{}{}}0000000{cC`{}}000{c{{Bb{e}}}{}{}}0000000{cCb{}}000`33={{{Ch{{Cf{Cd}}}}Bl}{{A`{b}}}}{{b{Bh{Bf}}{Cj{Bl}}}{{Cl{Al}}}}5332`{Cn{{A`{{Ch{{Cf{Cd}}}}}}}}{D`D`}077{{}D`}0{{CnAn}Bd}{cc{}}{Db{{Bb{CnDd}}}}0{{}{{Dh{Df}}}}<{Cn{{A`{{Dj{{Cf{Ad}}}}}}}}{CnDl}8{CnDn}9==<{{CnDb}{{Bb{nDd}}}}0{CnBl};{CnE`}```````````{ce{}{}}0{{EbAn}Bd}0:{EdEb}{EfEb}{EhEb}{EjEb}5{Eb{{Dh{El}}}}{cC`{}}{c{{Bb{e}}}{}{}}0{cCb{}}````999999{{{En{ce}}}{{En{ce}}}{F`FbF`Fd}F`}{{{Ff{ce}}}{{Ff{ce}}}{F`FbF`}F`}{{{Fh{ce}}}{{Fh{ce}}}{F`FbF`Fd}F`}{{ce}n{}{}}00{{{Ff{ce}}}e{FbF`}{}}{{}{{Fh{ce}}}{FbF`Fd}{}}`{{{En{ce}}{En{ce}}}Dl{FjFbF`Fd}Fj}{{{Ff{ce}}{Ff{ce}}}Dl{FjFbF`}Fj}{{{Fh{ce}}{Fh{ce}}}Dl{FjFbF`Fd}Fj}{{{Fh{ce}}g}{{Cl{{Ff{ce}}}}}{FbF`Fd}F`{{Fl{{Cj{c}}}}}}{{{Fh{ce}}g{Cl{{Ff{ce}}}}}n{FbF`Fd}F`{{Fl{{Cj{c}}}}}}{{{En{ce}}An}Bd{FnFbF`Fd}Fn}{{{Ff{ce}}An}Bd{FnFbF`}Fn}{{{Fh{ce}}An}Bd{FnFbF`Fd}Fn}{cc{}}00{i{{Fh{cg}}}{FbF`Fd}{{Fl{{Cj{c}}}}}F`{{Gd{}{{G`{{Gb{eg}}}}}}}}{{{En{ce}}g}n{GfFbF`Fd}GfGh}{{{Ff{ce}}g}n{GfFbF`}GfGh}{{{Fh{ce}}g}n{GfFbF`Fd}GfGh}{{{Fh{ce}}}n{FbF`Fd}F`}{{{Fh{ce}}ge}n{FbF`Fd}F`{{Fl{{Cj{c}}}}}}{{{Ff{ce}}}{{Cj{c}}}{FbF`}{}}{ce{}{}}00{{}{{Fh{ce}}}{FbF`Fd}F`}`111{c{{Bb{e}}}{}{}}00000{cCb{}}00``````````````````````````````````````{{jGj}n}{{GlGn}n}5555555555{Gl{{Bb{GjH`}}}}{{}Gl}{GnGn}{GjGj}{{ce}n{}{}}03{{Gl{Cl{Gn}}}Gl}`{{Gl{Bh{Gn}}}n}{{GnAn}Bd}{{H`An}Bd}0{{GjAn}Bd}{{jAn}Bd}{cc{}}00{HbH`}{C`H`}22{hj}{ce{}{}}0000{{Gnc}BbHd}{{Gjc}BbHd}{{jc}BbHd}{{Gl{Cl{Bf}}}Gl}`44{cC`{}}1`{c{{Bb{e}}}{}{}}000000000{cCb{}}0000`77;{{{Ch{{Cf{Cd}}}}}{{A`{d}}}}{{d{Bh{Bf}}{Cj{Bl}}}{{Dh{{Bh{Bf}}}}}}{{d{Bh{Bf}}{Bh{{Gb{{Cj{Bl}}Af}}}}}{{Cl{Bf}}}}:{{{Bh{Bf}}}n}554`;;?{{{Ch{{Cf{Cd}}}}}{{A`{f}}}}{{f{Bh{Bf}}}Bf}={{f{Bh{Bf}}}Dl}0{BfBf}{{{Bh{Bf}}}Bf}{{f{Bh{Bf}}}{{Cl{Bf}}}};;:```{ce{}{}}000{hh}{{ce}n{}{}}{{hAn}{{Bb{nB`}}}}{cc{}}0{Bj{{A`{h}}}}{c{{l{c}}}Ab}{h{{Gb{{Bh{Bf}}{Bn{Bl}}}}}}777{{{l{c}}}{{Dh{e}}}Ab{}}```{{hc}BbHd}9{c{{Bb{e}}}{}{}}000{cCb{}}0","D":"AFj","p":[[5,"AnnotationsDataBase",86],[5,"SequencesDataBase",306],[5,"Translate",318],[5,"Variant",333],[5,"Myth",203],[5,"VcfReader",333],[1,"unit"],[8,"Result",123],[10,"BufRead",363],[10,"Write",363],[6,"Strand",11],[6,"Frame",11],[5,"Attribute",11],[5,"Annotation",11],[5,"Formatter",364],[5,"Error",364],[6,"Result",365],[8,"Result",364],[1,"u8"],[1,"slice"],[5,"ByteRecord",366],[1,"u64"],[5,"Range",367],[5,"String",368],[5,"TypeId",369],[10,"Read",363],[5,"Box",370],[5,"BufReader",371],[8,"Interval",149],[5,"Vec",372],[5,"Command",96],[5,"Command",373],[5,"ArgMatches",374],[8,"Error",375],[5,"Id",376],[6,"Option",377],[5,"BufWriter",378],[1,"bool"],[6,"Timestamp",379],[1,"usize"],[6,"Error",123],[5,"SetLoggerError",380],[5,"Error",381],[5,"ParseIntError",382],[5,"ParseFloatError",383],[10,"Error",384],[5,"InternalEntry",149],[10,"Clone",385],[10,"Ord",386],[10,"Copy",387],[5,"Entry",149],[5,"IntervalTree",149],[10,"PartialEq",386],[10,"Into",388],[10,"Debug",364],[17,"Item"],[1,"tuple"],[10,"IntoIterator",389],[10,"Hash",390],[10,"Hasher",390],[5,"AnnotationMyth",203],[5,"AnnotationMythBuilder",203],[6,"Effect",203],[6,"AnnotationMythBuilderError",203],[5,"UninitializedFieldError",391],[10,"Serializer",392]],"r":[],"b":[[38,"impl-Display-for-Strand"],[39,"impl-Debug-for-Strand"],[40,"impl-Debug-for-Frame"],[41,"impl-Display-for-Frame"],[42,"impl-Debug-for-Attribute"],[43,"impl-Display-for-Attribute"],[44,"impl-Debug-for-Annotation"],[45,"impl-Display-for-Annotation"],[136,"impl-Debug-for-Error"],[137,"impl-Display-for-Error"],[139,"impl-From%3CSetLoggerError%3E-for-Error"],[140,"impl-From%3CError%3E-for-Error"],[141,"impl-From%3CParseIntError%3E-for-Error"],[142,"impl-From%3CParseFloatError%3E-for-Error"],[264,"impl-Debug-for-AnnotationMythBuilderError"],[265,"impl-Display-for-AnnotationMythBuilderError"],[271,"impl-From%3CUninitializedFieldError%3E-for-AnnotationMythBuilderError"],[272,"impl-From%3CString%3E-for-AnnotationMythBuilderError"]],"c":"OjAAAAAAAAA=","e":"OzAAAAEAALYAHgAWABgAQwATAFgAAQBeAAIAYwAGAGsAAgB0AAQAhwADAIwAAwCRAAQAmgALAKcAAACpAAIArgACALQAAwDAAAsA9AAJAAABBAAIAQQAEAEBABoBAgAfAQIAJAEOADQBAQA8AQIAQAEBAEsBAgBRAQYAXwEBAGQBBwA="}]\
]'));
if (typeof exports !== 'undefined') exports.searchIndex = searchIndex;
else if (window.initSearch) window.initSearch(searchIndex);
