package main

import "testing"
import (
	"reflect"
	"os"
	"github.com/brentp/vcfgo"
	"encoding/json"
	"bytes"
	"bufio"
)

func TestAnnfieldMultiEffect(t *testing.T){
	val:=annfield([]string{
		"G|missense_variant&splice_region_variant|"+
			"MODERATE|PIK3CD|PIK3CD|transcript|NM_005026.3|"+
				"protein_coding|4/24|c.368A>G|p.Lys123Arg|576/5411|368/3135|123/1044||"})
	correct:=[]map[string]interface{}{
		{
			"ANN_allele":"G",
			"ANN_effect":"missense_variant",
			"ANN_impact":"MODERATE",
			"ANN_gene_name":"PIK3CD",
			"ANN_gene_id":"PIK3CD",
			"ANN_feature_type":"transcript",
			"ANN_feature_id":"NM_005026.3",
			"ANN_transcript_biotype":"protein_coding",
			"ANN_rank_total":"4/24",
			"ANN_hgvs_c":"c.368A>G",
			"ANN_hgvs_p":"p.Lys123Arg",
			"ANN_cdna_position":"576/5411",
			"ANN_cds_position":"368/3135",
			"ANN_protein_position":"123/1044",
			//"ANN_distance_to_feature":"",
			//"ANN_errors_warnings_info":"",
		},
		{
			"ANN_allele":"G",
			"ANN_effect":"splice_region_variant",
			"ANN_impact":"MODERATE",
			"ANN_gene_name":"PIK3CD",
			"ANN_gene_id":"PIK3CD",
			"ANN_feature_type":"transcript",
			"ANN_feature_id":"NM_005026.3",
			"ANN_transcript_biotype":"protein_coding",
			"ANN_rank_total":"4/24",
			"ANN_hgvs_c":"c.368A>G",
			"ANN_hgvs_p":"p.Lys123Arg",
			"ANN_cdna_position":"576/5411",
			"ANN_cds_position":"368/3135",
			"ANN_protein_position":"123/1044",
			//"ANN_distance_to_feature":"",
			//"ANN_errors_warnings_info":"",
		},
		}
	if !reflect.DeepEqual(val,correct){
		t.Log(val)
		t.Log(correct)
		t.Error("Annfield Multieffect Test Failed")
	}
}

func TestAnnfieldMultiAnnotation(t *testing.T){
	val:=annfield([]string{
		"G|missense_variant&splice_region_variant|"+
			"MODERATE|PIK3CD|PIK3CD|transcript|NM_005026.3|"+
			"protein_coding|4/24|c.368A>G|p.Lys123Arg|576/5411|368/3135|123/1044||",
		"G|missense_variant&splice_region_variant|"+
			"MODERATE|PIK3CD|PIK3CD|transcript|NM_005026.3|"+
			"protein_coding|4/24|c.368A>G|p.Lys123Arg|576/5411|368/3135|123/1044||"})
	correct:=[]map[string]interface{}{
		{
			"ANN_allele":"G",
			"ANN_effect":"missense_variant",
			"ANN_impact":"MODERATE",
			"ANN_gene_name":"PIK3CD",
			"ANN_gene_id":"PIK3CD",
			"ANN_feature_type":"transcript",
			"ANN_feature_id":"NM_005026.3",
			"ANN_transcript_biotype":"protein_coding",
			"ANN_rank_total":"4/24",
			"ANN_hgvs_c":"c.368A>G",
			"ANN_hgvs_p":"p.Lys123Arg",
			"ANN_cdna_position":"576/5411",
			"ANN_cds_position":"368/3135",
			"ANN_protein_position":"123/1044",
			//"ANN_distance_to_feature":"",
			//"ANN_errors_warnings_info":"",
		},
		{
			"ANN_allele":"G",
			"ANN_effect":"splice_region_variant",
			"ANN_impact":"MODERATE",
			"ANN_gene_name":"PIK3CD",
			"ANN_gene_id":"PIK3CD",
			"ANN_feature_type":"transcript",
			"ANN_feature_id":"NM_005026.3",
			"ANN_transcript_biotype":"protein_coding",
			"ANN_rank_total":"4/24",
			"ANN_hgvs_c":"c.368A>G",
			"ANN_hgvs_p":"p.Lys123Arg",
			"ANN_cdna_position":"576/5411",
			"ANN_cds_position":"368/3135",
			"ANN_protein_position":"123/1044",
			//"ANN_distance_to_feature":"",
			//"ANN_errors_warnings_info":"",
		},
		{
			"ANN_allele":"G",
			"ANN_effect":"missense_variant",
			"ANN_impact":"MODERATE",
			"ANN_gene_name":"PIK3CD",
			"ANN_gene_id":"PIK3CD",
			"ANN_feature_type":"transcript",
			"ANN_feature_id":"NM_005026.3",
			"ANN_transcript_biotype":"protein_coding",
			"ANN_rank_total":"4/24",
			"ANN_hgvs_c":"c.368A>G",
			"ANN_hgvs_p":"p.Lys123Arg",
			"ANN_cdna_position":"576/5411",
			"ANN_cds_position":"368/3135",
			"ANN_protein_position":"123/1044",
			//"ANN_distance_to_feature":"",
			//"ANN_errors_warnings_info":"",
		},
		{
			"ANN_allele":"G",
			"ANN_effect":"splice_region_variant",
			"ANN_impact":"MODERATE",
			"ANN_gene_name":"PIK3CD",
			"ANN_gene_id":"PIK3CD",
			"ANN_feature_type":"transcript",
			"ANN_feature_id":"NM_005026.3",
			"ANN_transcript_biotype":"protein_coding",
			"ANN_rank_total":"4/24",
			"ANN_hgvs_c":"c.368A>G",
			"ANN_hgvs_p":"p.Lys123Arg",
			"ANN_cdna_position":"576/5411",
			"ANN_cds_position":"368/3135",
			"ANN_protein_position":"123/1044",
			//"ANN_distance_to_feature":"",
			//"ANN_errors_warnings_info":"",
		},
	}
	if !reflect.DeepEqual(val,correct){
		t.Log(val)
		t.Log(correct)
		t.Error("Annfield MultiAnnotation Test Failed")
	}
}

func Test_parse_vcf_field_array_string(t *testing.T){
	val:=parse_vcf_field_array("COSM898720,COSM1683737","String")
	correct:=[]interface{}{"COSM898720","COSM1683737"}
	if(!reflect.DeepEqual(val,correct)){
		t.Log(val)
		t.Log(correct)
		t.Error("Parse vcf field array string Test Failed")
	}
}
func Test_parse_vcf_field_array_int(t *testing.T){
	val:=parse_vcf_field_array("1,2","Integer")
	correct:=[]interface{}{1,2}
	if(!reflect.DeepEqual(val,correct)){
		t.Log(val)
		t.Log(correct)
		t.Error("Parse vcf field array int Test Failed")
	}
}
func Test_parse_vcf_field_array_float(t *testing.T){
	val:=parse_vcf_field_array("1.0,2.1","Float")
	correct:=[]interface{}{1.0,2.1}
	if(!reflect.DeepEqual(val,correct)){
		t.Log(val)
		t.Log(correct)
		t.Error("Parse vcf field array float Test Failed")
	}
}
func Test_parse_vcf_record(t *testing.T){
	f, err := os.Open("test-h.vcf")
	if err != nil {
		t.Error("Cannot load test vcf")
	}
	vr, err := vcfgo.NewReader(f, false)
	if err != nil {
		t.Error("Cannot create vcf reader")
	}
	var out bytes.Buffer
	out.Grow(4096)
	var o = bufio.NewWriter(&out)
	var encoder=json.NewEncoder(o)
	var v=vr.Read()
	t.Log(v.String())
	parse_vcf_record(v,encoder,false)
	o.Flush()
	//var res interface{}
	//json.NewDecoder(out).Decode(&res)
	t.Log(out.String())
}

func Test_parse_vcf_record2(t *testing.T){
	f, err := os.Open("index.vcf")
	if err != nil {
		t.Error("Cannot load test vcf")
	}
	vr, err := vcfgo.NewReader(f, false)
	if err != nil {
		t.Error("Cannot create vcf reader")
	}
	var out bytes.Buffer
	out.Grow(4096)
	var o = bufio.NewWriter(&out)
	var encoder=json.NewEncoder(o)
	var v=vr.Read()
	t.Log(v.String())
	parse_vcf_record(v,encoder,false)
	o.Flush()
	//var res interface{}
	//json.NewDecoder(out).Decode(&res)
	t.Log(out.String())
}