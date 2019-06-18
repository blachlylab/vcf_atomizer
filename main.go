package main

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"flag"
	"io"
	"math"
	"os"
	"reflect"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
)

func annfield(anns []string) []map[string]interface{} {
	var field_list = []string{"allele",
		"effect",
		"impact",
		"gene_name",
		"gene_id",
		"feature_type",
		"feature_id",
		"transcript_biotype",
		"rank_total",
		"hgvs_c",
		"hgvs_p",
		"cdna_position",
		"cds_position",
		"protein_position",
		"distance_to_feature",
		"errors_warnings_info"}

	//Decodes one or more annotations from VCF ANN field standard 1.0
	//http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
	//input: annotation string (either with or without ANN= prefix)
	//output: (Generator)
	//one or more dictionaries representing annotations
	//absent elements in the ANN field will not have
	//a corresponding key in the dict

	// if empty string
	if len(anns) < 1 {
		return []map[string]interface{}{}
	}
	// specs page 2:
	// Multiple effects / consequences are separated by comma.
	var ann_string string
	var ann_arr []string
	var res = make([]map[string]interface{}, 0)
	var effects []string
	for _, ann_string = range anns {
		ann_arr = strings.Split(ann_string, "|")
		effects = strings.Split(ann_arr[1], "&")
		var i int
		var val, effect string
		for _, effect = range effects {
			var eff_map = make(map[string]interface{})
			for i, val = range ann_arr {
				if i == 1 {
					eff_map["ANN_"+field_list[i]] = effect
				} else if val != "" {
					eff_map["ANN_"+field_list[i]] = val
				}
			}
			res = append(res, eff_map)
		}
	}
	return res
}

//unpacks maps into main map
//similar to ** operator in python
func unpack(main map[string]interface{}, maps ...map[string]interface{}) {
	var temp map[string]interface{}
	for _, temp = range maps {
		var k string
		var v interface{}
		for k, v = range temp {
			main[k] = v
		}
	}
}

func vcf_transform(filename string, mapping string, meta string, sr bool) {
	//Opens vcf and loops over rows
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	var r io.Reader
	r, err = gzip.NewReader(f)
	if err != nil {
		err=f.Close()
		if err!=nil{
			panic(err)
		}
		f, err = os.Open(filename)
		if err!=nil{
			panic(err)
		}
		r = io.Reader(f)
	}
	vr, err := vcfgo.NewReader(r, false)
	if err != nil {
		panic(err)
	}
	var variant *vcfgo.Variant
	var out = bufio.NewWriter(os.Stdout)
	var encoder = json.NewEncoder(out)
	if mapping != "" {
		j, _ := os.Create(mapping)
		var map_out = bufio.NewWriter(j)
		var map_encoder = json.NewEncoder(map_out)
		json_map := make(map[string]map[string]map[string]string)
		json_map["properties"] = make(map[string]map[string]string)
		for key, val := range vr.Header.Infos {
			key = "INFO_" + key
			json_map["properties"][key] = make(map[string]string)
			switch val.Type {
			case "Integer":
				json_map["properties"][key]["type"] = "long"
			case "Float":
				json_map["properties"][key]["type"] = "float"
			case "Flag":
				json_map["properties"][key]["type"] = "boolean"
			default:
				delete(json_map["properties"], key)
			}

		}
		err=map_encoder.Encode(json_map)
		if err!=nil{
			panic(err)
		}
	}
	if meta != "" {
		fi, _ := os.Create(meta)
		var meta_out = bufio.NewWriter(fi)
		var meta_encoder = json.NewEncoder(meta_out)
		for _, s := range vr.Header.SampleNames {
			meta_map := make(map[string]interface{})
			meta_map["sample"] = s
			meta_map["type"] = "variant_meta"
			meta_map["extras"] = vr.Header.Extras
			meta_map["file_format"] = vr.Header.FileFormat
			meta_map["filters"] = vr.Header.Filters
			meta_map["contigs"] = vr.Header.Contigs
			meta_map["infos"] = vr.Header.Infos
			err=meta_encoder.Encode(meta_map)
			if err!=nil{
				panic(err)
			}
		}
		err=meta_out.Flush()
		if err!=nil{
			panic(err)
		}
		err=fi.Close()
		if err!=nil{
			panic(err)
		}
	}
	for {
		variant = vr.Read()
		if variant == nil {
			break
		}
		parse_vcf_record(variant, encoder, sr)

	}
	err=out.Flush()
	if err!=nil{
		panic(err)
	}
}
func parse_vcf_field_array(val string, val_type string) []interface{} {
	var str_arr = strings.Split(val, ",")
	var ret []interface{}
	var str string
	for _, str = range str_arr {
		switch val_type {
		case "Integer":
			var i, err = strconv.Atoi(str)
			if err != nil {
				panic(err)
			} else {
				ret = append(ret, i)
			}
		case "Float":
			var i, err = strconv.ParseFloat(str, 64)
			if err != nil {
				panic(err)
			}
			if math.IsInf(i,0) || math.IsNaN(i) {
				ret = append(ret, nil)
			}else{
				ret = append(ret, i)
			}
		case "String":
			ret = append(ret, str)
		default:
			panic(val_type)
		}
	}
	return ret
}

func parse_vcf_record(variant *vcfgo.Variant, encoder *json.Encoder, sr bool) {
	//assign fields common to the row
	var common_fields = make(map[string]interface{})
	common_fields["type"] = "variant_vcf"
	common_fields["CHROM"] = variant.Chromosome
	common_fields["POS"] = variant.Pos
	common_fields["REF"] = variant.Reference
	common_fields["QUAL"] = variant.Quality
	common_fields["ID"] = variant.Id_
	common_fields["FILTER"] = strings.Split(variant.Filter, ";")
	var info_key string
	var anns []map[string]interface{}
	var ann_strings []string
	//parse the info fields
	var err error
	if !(len(variant.Info().Keys()) == 1 && variant.Info().Keys()[0] == "") {
		for _, info_key = range variant.Info().Keys() {
			var res, _ = variant.Info().Get(info_key)
			var s = reflect.ValueOf(res)
			if info_key == "ANN" {
				ann_strings = make([]string, s.Len())
				for i := 0; i < s.Len(); i++ {
					ann_strings[i] = s.Index(i).String()
				}
				if s.Type().Name() == "string" {
					anns = annfield([]string{s.String()})
				} else {
					anns = annfield(ann_strings)
				}
			} else {
				if variant.Header.Infos[info_key].Number != "1" {
					common_fields["INFO_"+info_key] = res
				} else {
					switch key_type := variant.Header.Infos[info_key].Type; key_type {
					case "Integer":
						common_fields["INFO_"+info_key] = s.Int()
					case "Float":
						if math.IsInf(s.Float(),0) || math.IsNaN(s.Float()) {
							common_fields["INFO_"+info_key] = nil
						}else {
							common_fields["INFO_"+info_key] = s.Float()
						}
					case "String":
						common_fields["INFO_"+info_key] = s.String()
					case "Flag":
						common_fields["INFO_"+info_key] = true
					default:
						panic(res)
					}
				}
			}
		}
	}
	var alt string
	var index int
	//loop over alts
	for index, alt = range variant.Alt() {
		var alt_fields = make(map[string]interface{})
		alt_fields["ALT"] = alt
		if len(variant.Samples) < 1 {
			if len(anns) < 1 {
				unpack(alt_fields, common_fields)
				err=encoder.Encode(alt_fields)
				if err!=nil{
					panic(err)
				}
				continue
			} else {
				var ann map[string]interface{}
				for _, ann = range anns {
					if reflect.ValueOf(ann["ANN_allele"]).String() == alt {
						unpack(ann, common_fields, alt_fields)
						err=encoder.Encode(ann)
						if err!=nil{
							panic(err)
						}
					}
				}
				continue
			}

		}
		//loop over samples
		for samidx, sample := range variant.Samples {
			var found = false
			var gt int
			for _, gt = range sample.GT {
				if index+1 == gt {
					found = true
				}
			}
			if !found && len(sample.GT)!=0{
				continue
			}
			var sample_fields = make(map[string]interface{})
			sample_fields["sample"] = variant.Header.SampleNames[samidx]
			var key, val string
			var err error
			for key, val = range sample.Fields {
				if key != "DP" && key != "GT" && key != "MQ" && key != "GL" && key != "GQ" && key != "AD" {
					if variant.Header.SampleFormats[key].Number != "1" {
						sample_fields[key] = parse_vcf_field_array(val, variant.Header.SampleFormats[key].Type)
					} else {
						switch key_type := variant.Header.SampleFormats[key].Type; key_type {
						case "Integer":
							sample_fields[key], err = strconv.Atoi(val)
						case "Float":
							if val != "." {
								var i float64
								i, err = strconv.ParseFloat(val, 64)
								if math.IsInf(i,0) || math.IsNaN(i) {
									sample_fields[key]=nil
								}else{
									sample_fields[key]=i
								}
							}
						case "String":
							sample_fields[key] = val
						default:
							panic(key_type)
						}
					}
					if err != nil {
						panic(err)
					}
				}
			}
			sample_fields["DP"] = sample.DP
			if len(sample.GT) > 0 {
				sample_fields["GT"] = sample.GT
			}
			sample_fields["MQ"] = sample.MQ
			if len(sample.GL) > 0 {
				sample_fields["GL"] = sample.GL
			}
			sample_fields["GQ"] = sample.GQ
			sample_fields["Ref_Depth"], _ = sample.RefDepth()
			sample_fields["Alt_depths"], _ = sample.AltDepths()
			var ann map[string]interface{}
			//for each sample loop over anns and link
			if len(anns) < 1 {
				unpack(sample_fields, common_fields, alt_fields)
				err=encoder.Encode(sample_fields)
			}
			if sr {
				unpack(sample_fields, common_fields, alt_fields)
				sample_fields["ANN"] = anns
				err=encoder.Encode(sample_fields)
			} else {
				for _, ann = range anns {
					if reflect.ValueOf(ann["ANN_allele"]).String() == alt {
						unpack(ann, sample_fields, common_fields, alt_fields)
						err=encoder.Encode(ann)
					}
					if err !=nil {
						panic(err)
					}
				}

			}
			if err !=nil {
				panic(err)
			}
		}
	}
}

func main() {
	var mapping = flag.String("mapping", "", "print predicted elasticsearch mapping")
	var meta = flag.String("meta", "", "write metadata to file")
	var sr = flag.Bool("one", false, "do not permute rows by annotations")
	flag.Parse()
	var filename = flag.Arg(0)
	if(filename!=""){
		vcf_transform(filename, *mapping, *meta, *sr)
	}else {
		flag.PrintDefaults()
	}
}
