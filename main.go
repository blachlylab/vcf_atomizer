package main

import (
	"github.com/brentp/vcfgo"
	"os"
	"io"
	"encoding/json"
	"bufio"
	"reflect"
	"flag"
	"strings"
	"fmt"
)

func annfield(anns []string)[]map[string]interface{}{
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
		"errors_warnings_info" }

	//Decodes one or more annotations from VCF ANN field standard 1.0
	//http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
	//input: annotation string (either with or without ANN= prefix)
	//output: (Generator)
	//one or more dictionaries representing annotations
	//absent elements in the ANN field will not have
	//a corresponding key in the dict

	// if empty string
	if len(anns)<1{
		return []map[string]interface{}{}
	}
	// specs page 2:
	// Multiple effects / consequences are separated by comma.
	var eff_string string
	var effs []string
	var res= make([]map[string]interface{},len(anns))
	var index int
	for index,eff_string= range anns{
		var eff_map=make(map[string]interface{})
		effs=strings.Split(eff_string,"|")
		var i int
		var val string
		for i,val=range effs{
			if i>=len(effs){
				fmt.Println()
				fmt.Println(effs[0])
				os.Exit(0)
			}
			eff_map["ANN_"+field_list[i]]=val
		}
		res[index]=eff_map
	}
	return res
}

//unpacks maps into main map
//similar to ** operator in python
func unpack(main map[string]interface{}, maps ...map[string]interface{}){
	var temp map[string]interface{}
	for _,temp=range maps{
		var k string
		var v interface{}
		for k,v=range temp{
			main[k]=v
		}
	}
}

func vcf_transform(filename string)  {
	//Opens vcf and loops over rows
	f, err := os.Open(filename)

	r := io.Reader(f)
	vr, err := vcfgo.NewReader(r, false)
	if err != nil {
		panic(err)
	}
	var variant *vcfgo.Variant
	var encoder=json.NewEncoder(bufio.NewWriter(os.Stdout))
	for {
		variant = vr.Read()
		if variant==nil{
			break
		}
		//if ANN in vcf use annfield record parser
		var _,ann_check=variant.Header.Infos["ANN"]
		if ann_check{
			parse_vcf_record_ANN(variant,encoder)
		}else {
			parse_vcf_record(variant,encoder)
		}

	}
}

func parse_vcf_record(variant *vcfgo.Variant,encoder *json.Encoder)  {
	var common_fields = make(map[string]interface{})
	common_fields["CHROM"]=variant.Chromosome
	common_fields["POS"]=variant.Pos
	common_fields["REF"]=variant.Reference
	common_fields["QUAL"]=variant.Quality
	common_fields["ID"]=variant.Id_
	common_fields["FILTER"]=strings.Split(variant.Filter,";")
	var info_key string
	for _,info_key=range variant.Info().Keys(){
		res,_:=variant.Info().Get(info_key)
		common_fields["INFO_"+info_key]=res
	}
	var v,key reflect.Value
	var t reflect.Type
	var alt string
	var index int
	for index,alt=range variant.Alt(){
		var alt_fields = make(map[string]interface{})
		alt_fields["ALT"]=alt
		if len(variant.Samples)<1{
			unpack(alt_fields,common_fields)
			encoder.Encode(alt_fields)
			continue
		}
		for _,sample:=range variant.Samples{
			var found = false
			var gt int
			for _,gt=range sample.GT{
				if index+1==gt{
					found=true
				}
			}
			if !found {
				continue
			}
			var sample_fields = make(map[string]interface{})
			sample_fields["sample"]=variant.Header.SampleNames[index]
			v=reflect.ValueOf(sample).Elem()
			t=v.Type()
			for i := 0; i < v.NumField(); i++ {
				if(t.Field(i).Name=="Fields"){
					//var field_v reflect.Value
					//var field_t reflect.Type
					var check bool
					for _,key=range v.Field(i).MapKeys(){
						_,check= sample_fields[key.String()]
						if !check{
							//fmt.Println(key.Interface())
							sample_fields[key.String()]=v.Field(i).MapIndex(key).Interface()
						}
					}


				}else{
					sample_fields[t.Field(i).Name]=v.Field(i).Interface()
				}
			}
			unpack(sample_fields,common_fields,alt_fields)
			encoder.Encode(sample_fields)
		}
	}
}

func parse_vcf_record_ANN(variant *vcfgo.Variant,encoder *json.Encoder)  {
	//assign fields common to the row
	var common_fields = make(map[string]interface{})
	common_fields["CHROM"]=variant.Chromosome
	common_fields["POS"]=variant.Pos
	common_fields["REF"]=variant.Reference
	common_fields["QUAL"]=variant.Quality
	common_fields["ID"]=variant.Id_
	common_fields["FILTER"]=strings.Split(variant.Filter,";")
	var info_key string
	var anns []map[string]interface{}
	var ann_strings []string
	var res interface{}
	var s reflect.Value
	//parse the info fields
	for _,info_key=range variant.Info().Keys(){
		if info_key=="ANN"{
			res,_=variant.Info().Get(info_key)
			s=reflect.ValueOf(res)
			ann_strings=make([]string, s.Len())
			for i := 0; i < s.Len(); i++{
				ann_strings[i]=s.Index(i).String()
			}
			if s.Type().Name()=="string" {
				anns=annfield([]string{s.String()})
			}else{
				anns=annfield(ann_strings)
			}
		}else{
			anns=[]map[string]interface{}{}
			res,_=variant.Info().Get(info_key)
			common_fields["INFO_"+info_key]=res
		}
	}
	var v,key reflect.Value
	var t reflect.Type
	var alt string
	var index int
	//loop over alts
	for index,alt=range variant.Alt(){
		var alt_fields = make(map[string]interface{})
		alt_fields["ALT"]=alt
		if len(variant.Samples)<1{
			unpack(alt_fields,common_fields)
			encoder.Encode(alt_fields)
			continue
		}
		//loop over samples
		for _,sample:=range variant.Samples{
			var found = false
			var gt int
			for _,gt=range sample.GT{
				if index+1==gt{
					found=true
				}
			}
			if !found {
				continue
			}
			var sample_fields = make(map[string]interface{})
			sample_fields["sample"]=variant.Header.SampleNames[index]
			v=reflect.ValueOf(sample).Elem()
			t=v.Type()
			for i := 0; i < v.NumField(); i++ {
				if(t.Field(i).Name=="Fields"){
					//var field_v reflect.Value
					//var field_t reflect.Type
					var check bool
					for _,key=range v.Field(i).MapKeys(){
						_,check= sample_fields[key.String()]
						if !check{
							//fmt.Println(key.Interface())
							sample_fields[key.String()]=v.Field(i).MapIndex(key).Interface()
						}
					}


				}else{
					sample_fields[t.Field(i).Name]=v.Field(i).Interface()
				}
			}
			var ann map[string]interface{}
			//for each sample loop over anns and link
			for _,ann=range anns{
				if (reflect.ValueOf(ann["ANN_allele"]).String()==alt){
					unpack(ann,sample_fields,common_fields,alt_fields)
					encoder.Encode(ann)
				}
			}
		}
	}
}

func main() {
	flag.Parse()
	var filename=flag.Arg(0)
	vcf_transform(filename)
}
