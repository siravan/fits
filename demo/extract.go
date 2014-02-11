// Copyright 2014 Shahriar Iravanian (siravan@svtsim.com).  All rights reserved.
// Use of this source code is governed by a MIT license that can be found in the LICENSE file.
//
// extract is a test application for the fits (or go-fits) package.
// It shows different usage of the fits package.
// extract accepts one input in command line, which can be a file name or a URL pointing to a FITS file
// extract reads the file and for each HDU writes the header key/value pairs to stdout
// In addition, it extracts the data for each HDU and writes it in the relevant format
// For images, the output format is 16-bit grayscale PNG
// For tables, extract generates a text file containing the data
// For one-dimensional images, it generates a simple one-column text file
//
package main

import (
	"bytes"	
	"fmt"
	"image"
	"image/color"
	"image/png"
	"io/ioutil"
	"log"
	"net/http"
	"net/url"
	"os"
	"path"
	"strings"    
    "fits"
)

func main() {
	var units []*fits.Unit
	var name string

	if len(os.Args) == 1 {
		fmt.Println("usage: extract filename|url")
		os.Exit(1)
	}

	if strings.HasPrefix(os.Args[1], "http://") { // called as "extract url"
		url, err := url.Parse(os.Args[1])
		if err != nil {
			log.Fatal(err)
		}
		name = path.Base(url.Path)
		res, err := http.Get(os.Args[1])
		if err != nil {
			log.Fatal(err)
		}
		buf, _ := ioutil.ReadAll(res.Body) // we download the whole FITS file first and then pass a buffered Reader to fits.Open
		res.Body.Close()
		units, err = fits.Open(bytes.NewReader(buf))
		if err != nil {
			log.Fatal(err)
		}
	} else { // called as "extract filename"
		name = path.Base(os.Args[1])
		reader, err := os.Open(os.Args[1])
		if err != nil {
			log.Fatal(err)
		}
		defer reader.Close()
		units, err = fits.Open(reader)
		if err != nil {
			log.Fatal(err)
		}
	}

	j := strings.LastIndex(name, ".") // name is set to the name of the FITS file without its extension
	if j != -1 {
		name = string(name[:j])
	}

	for i, h := range units { // for each HDU, extract the content
		fmt.Printf("******************** Header %d ********************\n", i)

		for key, value := range h.Keys { // First, write all key/value pairs
			fmt.Println(key, value)
		}

		out := fmt.Sprintf("%s_%d", name, i)

		if h.HasImage() { // Image type HDU (SIMPLE or XTENSION=IMAGE)
			writeImage(h, out)
		} else if h.HasTable() { // Table type HDU (XTENSION=TABLE or XTENSION=BINTABLE)
			if len(h.Naxis) == 1 { // One-dimensional table, write as an array
				writeArray(h, out)
			} else { // Two- or more dimensional, write as one or more PNG images
				writeTable(h, out)
			}
		} else {
			fmt.Println("Unsupported Header Data Unit")
		}
	}
}

// writeArray writes a one-dimensional image as a text file
func writeArray(h *fits.Unit, name string) {
	g, _ := os.Create(name + ".dat")
	defer g.Close()

	for i := 0; i < h.Naxis[0]; i++ {
		fmt.Fprintln(g, h.FloatAt(i))
	}
}

// writeImage generates PNG file(s) for image type HDUs
// If h contains a two-dimensional image, a single 16-bit normalized PNG file is generated
// For higher dimensional units, NAXIS3xNAXIS4x... different images (each NAXIS1xNAXIS2 in size) are generated
// For example, if h.Naxis=[512, 512, 2, 3] and name is "test_0", the following images are written:
//
//      test_0-0.0.png  contains pixels [0,0,0,0] to [511,511,0,0]
//      test_0-1.0.png  contains pixels [0,0,1,0] to [511,511,1,0]
//      test_0-0.1.png  contains pixels [0,0,0,1] to [511,511,0,1]
//      test_0-1.1.png  contains pixels [0,0,1,1] to [511,511,1,1]
//      test_0-0.2.png  contains pixels [0,0,0,2] to [511,511,0,2]
//      test_0-1.2.png  contains pixels [0,0,1,2] to [511,511,1,2]
//
func writeImage(h *fits.Unit, name string) {
	n := len(h.Naxis)
	maxis := make([]int, n)
	img := image.NewGray16(image.Rect(0, 0, h.Naxis[0], h.Naxis[1]))
	prod := 1
	for k := 2; k < n; k++ {
		prod *= h.Naxis[k]
	}
	min, max := h.Stats()

	for i := 0; i < prod; i++ {
		l := i
		s := name
		for k := 2; k < n; k++ {
			maxis[k] = l % h.Naxis[k]
			l = l / h.Naxis[k]
			s += fmt.Sprintf("-%d", maxis[k])
		}

		for x := 0; x < h.Naxis[0]; x++ {
			for y := 0; y < h.Naxis[1]; y++ {
				maxis[0] = x
				maxis[1] = y
				if !h.Blank(maxis...) {
					v := uint16((h.FloatAt(maxis...) - min) / (max - min) * 65535) // normalizes based on min and max in the whole image cube
					img.SetGray16(x, h.Naxis[1]-y, color.Gray16{v})
				} else {
					img.SetGray16(x, h.Naxis[1]-y, color.Gray16{0}) // blank pixel
				}
			}
		}

		g, _ := os.Create(s + ".png")
		defer g.Close()
		png.Encode(g, img)
	}
}

// writeTable generates a text file containing the table data of h 
// It processes both text (XTENSION=TABLE) and binary (XTENSION=BINTABLE) tables
func writeTable(h *fits.Unit, name string) {
	g, _ := os.Create(name + ".tab")
	defer g.Close()
	ncols := h.Keys["TFIELDS"].(int)

	label := "" // label is the list of field names/labels 
	for col := 0; col < ncols; col++ {
		ttype := h.Keys[fits.Nth("TTYPE", col+1)].(string)
		w := len(h.Format(col, 0)) // the label for each field is resized based on the size of the data on the first row of data 
		label += fmt.Sprintf("%-*.*s", w, w, ttype)
	}
	fmt.Fprintln(g, label)

	for row := 0; row < h.Naxis[1]; row++ {
		s := ""
		for col := 0; col < ncols; col++ {
			s += h.Format(col, row)
		}
		fmt.Fprintln(g, s)
	}
}
