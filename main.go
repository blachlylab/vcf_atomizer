package main

import (
	"github.com/brentp/vcfgo"
	"os"
	//"fmt"
	//"flag"
	"io"
	//"log"
	"fmt"
	"reflect"
)

func main() {
	//flag.Parse()
	//files := flag.Args()
	f, err := os.Open("/home/cgregory/vcfgo/examples/test.auto_dom.no_parents.vcf")

	r := io.Reader(f)
	vr, err := vcfgo.NewReader(r, false)
	if err != nil {
		panic(err)
	}
	var variant *vcfgo.Variant
	for {
		variant = vr.Read()
		if variant==nil{
			break
		}
		fmt.Println(variant.Chromosome)
		var info_key string
		for _,info_key=range variant.Info().Keys(){
			fmt.Print(info_key+" ")
			res,_:=variant.Info().Get(info_key)
			fmt.Print(reflect.TypeOf(res))
			fmt.Print(" ")
			fmt.Println(res)
		}
	}

}
