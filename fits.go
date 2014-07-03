// Copyright 2014 Shahriar Iravanian (siravan@svtsim.com).  All rights reserved.
// Use of this source code is governed by a MIT license that can be found in the LICENSE file.
//
// Package fits reads and processes FITS files. It is written in pure golang and is not a wrapper around another library or a direct translation of
// another library to golang. The main purpose is to provide a native golang solution to reading FITS file and to assess the suitability of golang for
// scientific and numerical applications.
//
// FITS is a common format for astronomical image and data.
// This package is based on version 3.0 of the FITS standard:
//  Pence W.D., Chiappetti L., Page C. G., Shaw R. A., Stobie E. Definition of the Flexible Image Transport System (FITS), version 3.0. A&A 524, A42 (2010)
//  http://www.aanda.org/articles/aa/abs/2010/16/aa15362-10/aa15362-10.html
//
// The following features are supported in the current version:
//      1. Images with all six different data format (byte, int16, int32, int64, float32, and float64)
//      2. Text and binary tables with atomic and fixed-size array elements
//
// The following features are not yet implemented:
//      1. Automatic application of BSCALE/BZERO
//      2. Random group structure
//      3. Variable length arrays in binary tables
//      4. World coordinate system
//
// Also note that currently this package provides only read capability and does not write/generate a FITS file.
//
// The basic usage of the package is by calling Open function. It accepts a reader that should provide a valid FITS file.
// The output is a []*fits.Unit, where Unit represents a Header/Data Unit (i.e. a header with the corresponding data).
// Unit provides a set of variables and functions to access the HDU data.
//
// Let 'test.fits' be a FITS file with two HDU. The first one is of type SIMPLE and contains a single two-dimensional image with the following parameters:
//
//      BITPIX  =  -32
//      NAXIS   =  2
//      NAXIS1  =  512
//      NAXIS2  =  256
//
// The second HDU contains a binary table (XTENSION=BINTABLE):
//
//      BITPIX  =  8
//      NAXIS   =  2
//      NAXIS1  =  100
//      NASIX2  =  5
//      TFIELDS =  10
//      TFORM1  =  E
//      TTYPE   =  FLUX
//      TDISP1  =  F10.4
//
// To read this file, we first call
//
//      units := fits.Open("test.fits")
//
// Now, units[0] points to the first HDU. We can access the header keys by using units.Keys map.
// For example, units[0].Keys["BITPIX"].(int) returns -32. Note that Keys stores interface{} and appropriate type-assertion needs to be done.
// Unit.Naxis returns a slice of integers ([]int) containing all NAXIS data. For example, units[0].Naxis is equal to [512, 256].
// We can access the image data points by using one of the three accessor functions: Unit.At, Unit.IntAt and Unit.FloatAt.
// Each function accepts NAXIS integer arguments and returns the pixel value at that location.
// Unit.At returns an interface{} and needs to be type-asserted before use. Unit.IntAt and Unit.FloatAt return int64 and float64, respectively.
//
// For table data, we use two other accessor functions: Field and Format.
// Field accepts one argument, col, that define a field. It can be 0-based int or a string.
// For example, units[1].Field(0) and units[1].Field("FLUX") both points to the same column.
// The return value of Field is another function, which is the actual accessor function and accepts one int argument representing a row.
// For example, units[1].Field("Flux")(1) returns the value of column "FLUX" in the second row of the table as interface{}.
// The following code populates a slice of float with the value of the FLUX column:
//
//      fn := units[1].Field("FLUX")
//      x := make([]float32, units[1].Naxis[1])     // note Naxis[0]=NAXIS1=length of a row, Naxis[1]=NAXIS2=number of rows
//      for row := range x {
//          x[row] = fn(row).(float32)
//      }
//
// Format function on the hand accepts two arguments, col (same as Field) and row and return a string formatted according to TDISP for the field.
// For example, if units[1].Field("Flux")(1) is equal to 987.654321, then units[1].Format("Flux", 1) returns "987.6543".
//
package fits

import (
	"bytes"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
	"sync"
)

// FieldFunc are the type of accessor functions returned by Unit.Field()
// FieldFunc is used to access the value of cells in a text or binary table (XTENSION=TABLE or XTENSION=BINTABLE)
type FieldFunc func(row int) interface{}

// Unit stored the header and data of a single HDU (Header Data Unit) as defined by FITS standard
// Data points to a flat array holding the HDU data
// Its type is []byte for tables and is determined by BITPIX for images:
//
//      BITPIX  Data
//      8       []byte
//      16      []int16
//      32      []int32
//      64      []int64
//      -32     []float32
//      -64     []float64
//
type Unit struct {
	Keys  map[string]interface{}
	Naxis []int // len(Naxis) is equal to the value of NASIX in the header
	// Naxis[k] is equal to NAXIS{k+1} in the header
	Data   interface{}
	list   []FieldFunc          // A slice to help with access to FieldFunc based on index
	fields map[string]FieldFunc // A map of FieldFunc (field-name => accessor-function)
	// field-name is based on TTYPE{k} keys in the header
	class string                     // class holds the type of the Header (SIMPLE, IMAGE, TABLE and BINTABLE)
	blank int                        // The value of BLANK key in the header
	At    func(a ...int) interface{} // Accessor function that returns the value of a pixel based on its coordinates
	// a... represents NAXIS integers corresponding to NAXIS1, NAXIS2,...
	// The return result type is interface{}. The concrete type is determined by BITPIX
	IntAt   func(a ...int) int64   // A helper accessor function that returns the pixel value as int64
	FloatAt func(a ...int) float64 // A helper accessor function that returns the pixel value as float64
	Blank   func(a ...int) bool    // returns true if pixel type is integral and the pixel pointed by a... is equal to blank,
	// or the pixel type is float and its value is NaN
}

// Reader is a buffered Reader implementation that works based on the FITS block structure (each 2880 bytes long)
type Reader struct {
	buf    []byte
	elem   []byte
	left   int
	right  int
	reader io.Reader
	eof    bool
}

// Field returns a FieldFunc corresponding to col
// If col is int, the col'th field is returned (note: col is 0 based, so col=1 means TFORM2)
// If col a string, the field with TDISP equal to col is returned
// Fields are held in a map (Unit.fields) based on their name (TDISP).
// In addition, for each field, an entry with key "#name" is added to Unit.fields to facilitate the search for TDISP based on the name
//
// Note: this function returns an accessor function, that needs to be called to obtain the actual cell value
// For example, assume h is a table. One of its column is named "ID" of type "J" (int32)
// To obtain the value of the cell located at the intersection of the third row (row=2) and column "ID", we write
//
//  fn := h.Field("ID")
//  val := fn(2).(int32)
//
func (h *Unit) Field(col interface{}) FieldFunc {
	var x FieldFunc
	var ok bool

	switch col.(type) {
	case int:
		n := col.(int)
		if n >= 0 && n < len(h.list) {
			return h.list[col.(int)]
		}
	case string:
		x, ok = h.fields[col.(string)]
		if ok {
			return x
		}
	}
	return func(int) interface{} {
		return nil
	}
}

// Format returns a formatted string based on the given col and row and TDISP of the col
// col can be an int or a string (same as Field)
// The return value is a string, which is obtained by
//      1. Finding the FieldFunc based on col
//      2. Running the FieldFunc by passing row as an argument
//      3. Applying format to the result
//
func (h *Unit) Format(col interface{}, row int) string {
	var fn FieldFunc
	var disp interface{}

	switch col.(type) {
	case int:
		n := col.(int)
		if n >= 0 && n < len(h.list) {
			fn = h.list[col.(int)]
			disp, _ = h.Keys[Nth("TDISP", n+1)]
		}
	case string:
		name := col.(string)
		fn, _ = h.fields[name]
		n := h.Keys["#"+name]
		disp, _ = h.Keys[Nth("TDISP", n.(int))]
	}

	if fn == nil {
		return ""
	}

	format := "%v" // default format

	w := 14
	if disp != nil {
		var code rune
		m := -1
		d := disp.(string)

		// accounts for ENw.d and ESw.d formats
		if len(d) > 1 && (d[1] == 'N' || d[1] == 'S') {
			d = string(d[0]) + string(d[2:]) // removes the second character from the format string
			// The standard allows to disregard this secondary format characters
		}

		fmt.Sscanf(d, "%c%d.%d", &code, &w, &m)

		switch code {
		case 'A':
			format = fmt.Sprintf("%%%d.%ds", w, w) // Aw -> %ws
		case 'I':
			format = fmt.Sprintf("%%%dd", w) // Iw -> %wd
		case 'B':
			format = fmt.Sprintf("%%%db", w) // Bw -> %wb, binary
		case 'O':
			format = fmt.Sprintf("%%%do", w) // Ow -> %wo, octal
		case 'Z':
			format = fmt.Sprintf("%%%dX", w) // Zw -> %wX, hexadecimal
		case 'F', 'D':
			if m != -1 {
				format = fmt.Sprintf("%%%d.%df", w, m) // Fw.d -> %w.df
			} else {
				format = fmt.Sprintf("%%%df", w) // Fw -> %wf
			}
		case 'E':
			if m != -1 {
				format = fmt.Sprintf("%%%d.%de", w, m) // Fw.d -> %w.df
			} else {
				format = fmt.Sprintf("%%%de", w) // Ew -> %we
			}
		case 'G':
			if m != -1 {
				format = fmt.Sprintf("%%%d.%dg", w, m) // Fw.d -> %w.df
			} else {
				format = fmt.Sprintf("%%%dg", w) // Gw -> %wg
			}
		}
	}

	return fmt.Sprintf(format, fn(row))
}

// HasImage returns true is the Unit is either SIMPLE or IMAGE and has the data for an actual image
func (h *Unit) HasImage() bool {
	return (h.class == "SIMPLE" || h.class == "IMAGE") && len(h.Naxis) > 0 && h.Naxis[0] > 0
}

// HasImage returns true is the Unit is either TABLE or BINTABLE and has the data for an actual table
func (h *Unit) HasTable() bool {
	return (h.class == "TABLE" || h.class == "BINTABLE")
}

// Bitpix is a helper function the simply returns BITPIX value in the header
func (h *Unit) Bitpix() int {
	return h.Keys["BITPIX"].(int)
}

// Stats returns the minimum and maximum values in the image data
func (h *Unit) Stats() (min float64, max float64) {
	prod := 1
	for _, x := range h.Naxis {
		prod *= x
	}
	if prod == 1 {
		return
	}

	min = math.MaxFloat64
	max = -math.MaxFloat64

	switch h.Bitpix() {
	case 8:
		for i := 0; i < prod; i++ {
			x := int(h.Data.([]byte)[i])
			if x != h.blank && float64(x) < min {
				min = float64(x)
			}
			if x != h.blank && float64(x) > max {
				max = float64(x)
			}
		}
	case 16:
		for i := 0; i < prod; i++ {
			x := int(h.Data.([]int16)[i])
			if x != h.blank && float64(x) < min {
				min = float64(x)
			}
			if x != h.blank && float64(x) > max {
				max = float64(x)
			}
		}
	case 32:
		for i := 0; i < prod; i++ {
			x := int(h.Data.([]int32)[i])
			if x != h.blank && float64(x) < min {
				min = float64(x)
			}
			if x != h.blank && float64(x) > max {
				max = float64(x)
			}
		}
	case 64:
		for i := 0; i < prod; i++ {
			x := int(h.Data.([]int64)[i])
			if x != h.blank && float64(x) < min {
				min = float64(x)
			}
			if x != h.blank && float64(x) > max {
				max = float64(x)
			}
		}
	case -32:
		for i := 0; i < prod; i++ {
			x := float64(h.Data.([]float32)[i])
			if !math.IsNaN(x) && x < min {
				min = x
			}
			if !math.IsNaN(x) && x > max {
				max = x
			}
		}
	case -64:
		for i := 0; i < prod; i++ {
			x := h.Data.([]float64)[i]
			if !math.IsNaN(x) && x < min {
				min = x
			}
			if !math.IsNaN(x) && x > max {
				max = x
			}
		}
	}
	return
}

// Open processes a FITS file provided as an io.Reader and returns a list of HDUs in the FITS file
// It is the main entry point of the fits package
func Open(reader io.Reader) (fits []*Unit, err error) {
	b := NewReader(reader)
	fits = make([]*Unit, 0, 5)
done:
	for !b.IsEOF() {
		h, err := b.NewHeader()
		if err != nil {
			err = nil // EOF, not an error?
			break
		}
		fits = append(fits, h)
		if _, ok := h.Keys["SIMPLE"]; ok {
			err = h.verifyPrimary()
			if err != nil {
				break
			}
			h.class = "SIMPLE"
			if len(h.Naxis) > 0 {
				if h.Naxis[0] == 0 { // Random Group Headers are not supported and are not processed further
					break done
				}
				err = h.loadData(b) // Imaging data
				if err != nil {
					break
				}
			}
		} else if xten, ok := h.Keys["XTENSION"].(string); ok {
			err = h.verifyExtension()
			if err != nil {
				break
			}
			h.class = xten
			switch xten {
			case "IMAGE":
				if len(h.Naxis) > 0 {
					err = h.loadData(b)
					if err != nil {
						break
					}
				}
			case "TABLE":
				err = h.loadTable(b, false)
				if err != nil {
					break
				}
			case "BINTABLE":
				err = h.loadTable(b, true)
				if err != nil {
					break
				}
			}
		} else {
			// unknown header
			break
		}
	}
	return fits, err
}

// index is a helper function the returns the index of the pixel pointed by a... in a flat Data array
func (h *Unit) index(a ...int) int {
	var index int
	for i := len(h.Naxis) - 1; i >= 0; i-- {
		index = index*h.Naxis[i] + a[i]
	}
	return index
}

// loadData processes the image type data sections
// It allocates Data, populates it, and sets the appropriate pixel accessor functions
func (h *Unit) loadData(b *Reader) error {
	var i int

	if len(h.Naxis) == 0 {
		h.Data = make([]int, 0)
		h.IntAt = func(a ...int) int64 {
			return 0
		}
		h.FloatAt = func(a ...int) float64 {
			return 0
		}
		return nil
	}

	prod := 1
	for _, x := range h.Naxis {
		prod *= x
	}

	bitpix := h.Keys["BITPIX"].(int)

	switch bitpix {
	case 8:
		data := make([]byte, prod) // Data type is determined based on bitpix
		h.Data = data
		h.At = func(a ...int) interface{} { // The accessor functions look similar, but note that data is redefined and has a different type for each case
			// Templates (generics) would have helped with cutting back on redundant code!
			return data[h.index(a...)]
		}
		h.IntAt = func(a ...int) int64 {
			return int64(data[h.index(a...)])
		}
		h.FloatAt = func(a ...int) float64 {
			return float64(data[h.index(a...)])
		}
		for i = 0; i < prod; i++ {
			data[i] = b.ReadByte()
		}
	case 16:
		data := make([]int16, prod)
		h.Data = data
		h.At = func(a ...int) interface{} {
			return data[h.index(a...)]
		}
		h.IntAt = func(a ...int) int64 {
			return int64(data[h.index(a...)])
		}
		h.FloatAt = func(a ...int) float64 {
			return float64(data[h.index(a...)])
		}
		for i = 0; i < prod; i++ {
			data[i] = b.ReadInt16()
		}
	case 32:
		data := make([]int32, prod)
		h.Data = data
		h.At = func(a ...int) interface{} {
			return data[h.index(a...)]
		}
		h.IntAt = func(a ...int) int64 {
			return int64(data[h.index(a...)])
		}
		h.FloatAt = func(a ...int) float64 {
			return float64(data[h.index(a...)])
		}
		for i = 0; i < prod; i++ {
			data[i] = b.ReadInt32()
		}
	case 64:
		data := make([]int64, prod)
		h.Data = data
		h.At = func(a ...int) interface{} {
			return data[h.index(a...)]
		}
		h.IntAt = func(a ...int) int64 {
			return int64(data[h.index(a...)])
		}
		h.FloatAt = func(a ...int) float64 {
			return float64(data[h.index(a...)])
		}
		for i = 0; i < prod; i++ {
			data[i] = b.ReadInt64()
		}
	case -32:
		data := make([]float32, prod)
		h.Data = data
		h.At = func(a ...int) interface{} {
			return data[h.index(a...)]
		}
		h.IntAt = func(a ...int) int64 {
			return int64(data[h.index(a...)])
		}
		h.FloatAt = func(a ...int) float64 {
			return float64(data[h.index(a...)])
		}
		for i = 0; i < prod; i++ {
			data[i] = b.ReadFloat32()
		}
	case -64:
		data := make([]float64, prod)
		h.Data = data
		h.At = func(a ...int) interface{} {
			return data[h.index(a...)]
		}
		h.IntAt = func(a ...int) int64 {
			return int64(data[h.index(a...)])
		}
		h.FloatAt = func(a ...int) float64 {
			return float64(data[h.index(a...)])
		}
		for i = 0; i < prod; i++ {
			data[i] = b.ReadFloat64()
		}
	}

	blank, ok := h.Keys["BLANK"]
	switch {
	case ok && bitpix > 0: // Integer pixel type with defined BLANK
		h.blank = blank.(int)
		h.Blank = func(a ...int) bool {
			return h.IntAt(a...) == int64(h.blank)
		}
	case bitpix < 0: // Float pixel type
		h.Blank = func(a ...int) bool {
			return math.IsNaN(h.FloatAt(a...))
		}
	default: // Integer pixel type with undefined BLANK
		h.Blank = func(a ...int) bool {
			return false
		}
	}

	return nil
}

// accessorBin generates the accessor function for a field in a binary table (XTENSION=BINTABLE)
// loadTable function processes TFORM for each field
// For binary tables, TFORM is like rT, where r is the repeat and T is the type code
// With the exception of code='A' (string-type), the accessor functions are different for repeat=1 (returns an atomic value) vs repeat>1 (returns a fixed array)
// Note, variable arrays (type P and Q) and packed bits (type X) are not supported in the current version
// col is the byte index of the value of the field from the beginning of each record
func (h *Unit) accessorBin(code byte, repeat int, col *int) (fn func(int) interface{}, disp string) {
	c := *col
	l := 0
	var f func() interface{} // f holds a helper function that returns the field data assuming that b is set correctly

	// we use a fits.Reader to access data values in the binary table
	b := new(Reader)
	b.buf = h.Data.([]byte)
	b.elem = make([]byte, 8)
	b.right = len(b.buf)

	switch code {
	case 'A':
		f = func() interface{} { // For T='A', the result is always a string, even if repeat is equal to 1
			return b.ReadString(repeat)
		}
		l = 1
		disp = fmt.Sprintf("A%d", repeat)
	case 'B':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadByte()
			}
		} else {
			f = func() interface{} {
				p := make([]uint8, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadByte()
				}
				return p
			}
		}
		l = 1
		disp = "I3" // disp is the default display formatting string to be used if the corresponding TDISP is missing
	case 'L':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadBool()
			}
		} else {
			f = func() interface{} {
				p := make([]bool, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadBool()
				}
				return p
			}
		}
		l = 1
		disp = "B1"
	case 'I':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadInt16()
			}
		} else {
			f = func() interface{} {
				p := make([]int16, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadInt16()
				}
				return p
			}
		}
		l = 2
		disp = "I6"
	case 'J':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadInt32()
			}
		} else {
			f = func() interface{} {
				p := make([]int32, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadInt32()
				}
				return p
			}
		}
		l = 4
		disp = "I11"
	case 'K':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadInt64()
			}
		} else {
			f = func() interface{} {
				p := make([]int64, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadInt64()
				}
				return p
			}
		}
		l = 8
		disp = "I20"
	case 'D':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadFloat64()
			}
		} else {
			f = func() interface{} {
				p := make([]float64, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadFloat64()
				}
				return p
			}
		}
		l = 8
		disp = "F14.7"
	case 'E':
		if repeat == 1 {
			f = func() interface{} {
				return b.ReadFloat32()
			}
		} else {
			f = func() interface{} {
				p := make([]float32, repeat)
				for i := 0; i < repeat; i++ {
					p[i] = b.ReadFloat32()
				}
				return p
			}
		}

		l = 4
		disp = "F14.7"
	case 'M':
		if repeat == 1 {
			f = func() interface{} {
				x := b.ReadFloat64()
				y := b.ReadFloat64()
				return complex(x, y)
			}
		} else {
			f = func() interface{} {
				p := make([]complex128, repeat)
				for i := 0; i < repeat; i++ {
					x := b.ReadFloat64()
					y := b.ReadFloat64()
					p[i] = complex(x, y)
				}
				return p
			}
		}
		l = 16
		disp = "F14.7"
	case 'C':
		if repeat == 1 {
			f = func() interface{} {
				x := b.ReadFloat32()
				y := b.ReadFloat32()
				return complex(x, y)
			}
		} else {
			f = func() interface{} {
				p := make([]complex64, repeat)
				for i := 0; i < repeat; i++ {
					x := b.ReadFloat32()
					y := b.ReadFloat32()
					p[i] = complex(x, y)
				}
				return p
			}
		}
		l = 8
		disp = "F14.7"
	case 'X', 'P', 'Q':
		panic("Binary table forms X, P and Q are not supported")
	}

	*col += l * repeat

	// fn is the actual FieldFunc
	// it sets b.left based on the record size and row and calls f to extract the field value
	fn = func(row int) interface{} {
		var m sync.Mutex
		m.Lock()                          // Lock is needed because each FieldFunc closes over a fits.Reader and b.left is modified
		if row < 0 || row >= h.Naxis[1] { // invalid row number (note Naxis[1] is NAXIS2 in the header equal to the number of rows)
			return nil
		}
		b.left = row*h.Naxis[0] + c
		x := f()
		m.Unlock()
		return x
	}

	return fn, disp
}

// accessorText generates the accessor function for a field in a text table (XTENSION=TABLE)
// loadTable function processes TFORM for each field
// For text tables, TFORM is like Tw or Tw.d (T=code and w=repeat)
func (h *Unit) accessorText(code byte, repeat int, col *int) (fn func(int) interface{}, disp string) {
	c := *col - 1
	var f func() interface{}
	b := new(Reader) // note that b.elem does not need to be set because we only use b.ReadString
	b.buf = h.Data.([]byte)
	b.right = len(b.buf)

	switch code {
	case 'A':
		f = func() interface{} {
			return b.ReadString(repeat)
		}
		disp = fmt.Sprintf("A%d", repeat)
	case 'I':
		f = func() interface{} {
			s := b.ReadString(repeat)
			s = strings.TrimSpace(s)
			n, _ := strconv.ParseInt(s, 10, 32)
			return int(n)
		}
		disp = fmt.Sprintf("I%d", repeat)
	case 'D', 'E', 'F':
		f = func() interface{} {
			s := b.ReadString(repeat)
			s = strings.TrimSpace(s)
			s = strings.Replace(s, "D", "E", 1)
			x, _ := strconv.ParseFloat(s, 64)
			return x
		}
		disp = "F14.7"
	default:
		panic("Unsupported TFORM in an Ascii table")
	}

	// same as fn function in accessorBin
	fn = func(row int) interface{} {
		var m sync.Mutex
		m.Lock()
		if row < 0 || row >= h.Naxis[1] {
			return nil
		}
		b.left = row*h.Naxis[0] + c
		x := f()
		m.Unlock()
		return x
	}

	return fn, disp
}

// verifyPrimary verifies a primary (SIMPLE) header for correctness and the presence of mandatory keys
func (h *Unit) verifyPrimary() error {
	_, ok := h.Keys["SIMPLE"]
	if !ok {
		return fmt.Errorf("No SIMPLE in the primary header")
	}
	n, ok := h.Keys["BITPIX"].(int)
	if !ok {
		return fmt.Errorf("No BITPIX in the primary header")
	}
	if n != 8 && n != 16 && n != 32 && n != 64 && n != -32 && n != -64 {
		return fmt.Errorf("Invalid BITPIX value")
	}
	n, ok = h.Keys["NAXIS"].(int)
	if !ok {
		return fmt.Errorf("No NAXIS in the primary header")
	}
	for i := 1; i <= n; i++ {
		s := Nth("NAXIS", i)
		_, ok := h.Keys[s].(int)
		if !ok {
			return fmt.Errorf("No %v in the primary header", s)
		}
	}
	return nil
}

// verifyExtension verifies a secondary (XTENSION) header for correctness and the presence of mandatory keys
func (h *Unit) verifyExtension() error {
	xten, ok := h.Keys["XTENSION"].(string)
	if !ok {
		return fmt.Errorf("No XTENSION in the extended header")
	}
	n, ok := h.Keys["BITPIX"].(int)
	if !ok {
		return fmt.Errorf("No BITPIX in the extended header")
	}
	if n != 8 && n != 16 && n != 32 && n != 64 && n != -32 && n != -64 {
		return fmt.Errorf("Invalid BITPIX value")
	}
	naxis, ok := h.Keys["NAXIS"].(int)
	if !ok {
		return fmt.Errorf("No NAXIS in the extended header")
	}
	for i := 1; i <= naxis; i++ {
		s := Nth("NAXIS", i)
		_, ok := h.Keys[s].(int)
		if !ok {
			return fmt.Errorf("No %v in the extended header", s)
		}
	}
	pcount, ok := h.Keys["PCOUNT"].(int)
	if !ok {
		return fmt.Errorf("No PCOUNT in the extended header")
	}
	_, ok = h.Keys["GCOUNT"].(int)
	if !ok {
		return fmt.Errorf("No GCOUNT in the extended header")
	}
	switch xten {
	case "IMAGE":
		if pcount != 0 {
			return fmt.Errorf("PCOUNT should be 0 in IMAGE header")
		}
	case "TABLE", "BINTABLE":
		if n != 8 {
			return fmt.Errorf("BITPIX should be 8 in TABLE/BINTABLE headers")
		}
		if naxis != 2 {
			return fmt.Errorf("NAXIS should be 2 in TABLE/BINTABLE headers")
		}
	}
	return nil
}

// loadTable processes a table (text or binary) data section
// it allocates and reads data
// for each field, it calls accessorBin or accessorText to obtain the corresponding accessor function and adds it to fields
func (h *Unit) loadTable(b *Reader, binary bool) error {
	tfields := h.Keys["TFIELDS"].(int) // # of fields
	h.list = make([]FieldFunc, tfields)
	h.fields = make(map[string]FieldFunc, tfields)

	data := make([]byte, h.Naxis[0]*h.Naxis[1])
	b.Read(data)
	h.Data = data

	var col int
	for i := 0; i < tfields; i++ {
		var fn FieldFunc
		var j int
		var disp string
		form := h.Keys[Nth("TFORM", i+1)].(string)

		if binary { // BINTABLE
			j = strings.IndexAny(form, "ABCDEIJKLMPQX")
			if j == -1 {
				return fmt.Errorf("TFROM has invalid format (binary)")
			}
			repeat := 1
			if j > 0 {
				r, _ := strconv.ParseInt(form[:j], 10, 32)
				repeat = int(r)
			}
			if repeat > 0 {
				fn, disp = h.accessorBin(form[j], repeat, &col)
			} else {
				continue
			}
		} else { // TABLE
			j = strings.Index(form, ".")
			if j == -1 {
				j = len(form)
			}
			r, _ := strconv.ParseInt(form[1:j], 10, 32)
			col = h.Keys[Nth("TBCOL", i+1)].(int)
			fn, disp = h.accessorText(form[0], int(r), &col)
		}

		h.list[i] = fn
		name, ok := h.Keys[Nth("TTYPE", i+1)]
		if ok {
			h.fields[name.(string)] = fn
			h.Keys["#"+name.(string)] = i + 1 // is used to find the index of a field if only its name is given
		} else {
			h.Keys[Nth("TTYPE", i+1)] = Nth("COL", i+1) // default name given to fields without a corresponding TTYPE
		}

		_, ok = h.Keys[Nth("TDISP", i+1)]
		if !ok {
			h.Keys[Nth("TDISP", i+1)] = disp // if TDISP is missing, the default disp is added to the header as a TDISP
		}
	}

	return nil
}

// NewReader generates a new fits.Reader that wraps the given reader
// 2880 is the standard FITS file block size
func NewReader(reader io.Reader) *Reader {
	p := new(Reader)
	p.buf = make([]byte, 2880)
	p.elem = make([]byte, 8)
	p.reader = reader
	return p
}

// Read populates p while taking care of the FITS file block structure
func (b *Reader) Read(p []byte) (n int, err error) {
	m := len(p)
	for {
		k := copy(p[n:], b.buf[b.left:b.right])
		n += k
		b.left += k
		if n == m {
			return n, nil
		}
		b.right, err = b.reader.Read(b.buf)
		b.left = 0
		if err != nil {
			if err == io.EOF {
				b.eof = true
				return n, nil
			}
			panic("Error in Reader.Read")
		}
	}
	return 0, fmt.Errorf("unreachable!")
}

// IsEOF returns if b is finished
func (b *Reader) IsEOF() bool {
	return b.eof
}

// NextPage skips the rest of the current 2880-byte block and reads the next block
func (b *Reader) NextPage() (buf []byte, err error) {
	b.right, err = b.reader.Read(b.buf)
	b.left = b.right
	return b.buf, err
}

func (b *Reader) ReadByte() byte {
	b.Read(b.elem[0:1])
	return b.elem[0]
}

func (b *Reader) ReadBool() bool {
	b.Read(b.elem[0:1])
	return b.elem[0] != 0
}

func (b *Reader) ReadString(n int) string {
	p := make([]byte, n)
	b.Read(p)
	return string(p)
}

// ReadInt16 reads an int16 encoded in big-endian binary
// Note that the FITS standard supports only big-endian binaries
func (b *Reader) ReadInt16() int16 {
	b.Read(b.elem[0:2]) // we need to copy into an elem buf instead of pointing directly to b.buf because
	// the target value may straddle a block boundary
	x := uint16(b.elem[1]) | uint16(b.elem[0])<<8
	return int16(x)
}

// ReadInt16 reads an int32 encoded in big-endian binary
func (b *Reader) ReadInt32() int32 {
	b.Read(b.elem[0:4])
	x := uint32(b.elem[3]) | uint32(b.elem[2])<<8 | uint32(b.elem[1])<<16 | uint32(b.elem[0])<<24
	return int32(x)
}

// ReadInt16 reads an int64 encoded in big-endian binary
func (b *Reader) ReadInt64() int64 {
	b.Read(b.elem[0:8])
	x := uint64(b.elem[7]) | uint64(b.elem[6])<<8 | uint64(b.elem[5])<<16 | uint64(b.elem[4])<<24 |
		uint64(b.elem[3])<<32 | uint64(b.elem[2])<<40 | uint64(b.elem[1])<<48 | uint64(b.elem[0])<<56
	return int64(x)
}

// ReadInt16 reads a float32 encoded in big-endian binary
func (b *Reader) ReadFloat32() float32 {
	b.Read(b.elem[0:4])
	x := uint32(b.elem[3]) | uint32(b.elem[2])<<8 | uint32(b.elem[1])<<16 | uint32(b.elem[0])<<24
	return math.Float32frombits(x)
}

// ReadInt16 reads a float64 encoded in big-endian binary
func (b *Reader) ReadFloat64() float64 {
	b.Read(b.elem[0:8])
	x := uint64(b.elem[7]) | uint64(b.elem[6])<<8 | uint64(b.elem[5])<<16 | uint64(b.elem[4])<<24 |
		uint64(b.elem[3])<<32 | uint64(b.elem[2])<<40 | uint64(b.elem[1])<<48 | uint64(b.elem[0])<<56
	return math.Float64frombits(x)
}

// Nth returns a string resulted from concatenation of prefix and n in string form
// it is a stateless helper function
func Nth(prefix string, n int) string {
	return fmt.Sprintf("%s%d", prefix, n)
}

// processString is utilized by NewHeader to process string-type values in the header
// it uses a 3-state machine to process double single quotes
func processString(s string) (string, error) {
	var buf bytes.Buffer

	state := 0
	for _, char := range s {
		quote := (char == '\'')
		switch state {
		case 0:
			if !quote {
				return "", fmt.Errorf("String does not start with a quote")
			}
			state = 1
		case 1:
			if quote {
				state = 2
			} else {
				buf.WriteRune(char)
				state = 1
			}
		case 2:
			if quote {
				buf.WriteRune(char)
				state = 1
			} else {
				return strings.TrimRight(buf.String(), " "), nil
			}
		}
	}
	return "", fmt.Errorf("String ends prematurely")
}

// NewHeader reads and processes the next header from the a the reader stream
// its main function is to populate Keys and setups Naxis
func (b *Reader) NewHeader() (h *Unit, err error) {
	Keys := make(map[string]interface{}, 50)
	h = &Unit{Keys: Keys}

	for {
		buf, err := b.NextPage()
		if err != nil {
			fmt.Println(err)
			return h, err
		}

	_lines:
		for i := 0; i < 36; i++ { // each FITS header block is comprised of up to 36 80-byte lines
			s := string(buf[i*80 : (i+1)*80])
			key := strings.TrimSpace(s[:8])
			if s[8:10] != "= " { // note that the standard is strict regarding the position of the '=' sign
				Keys[key] = nil
				continue
			}

			s = strings.TrimSpace(s[10:])

			if s == "" {
				Keys[key] = nil
				continue
			}

			first := rune(s[0])

			if first == '\'' {
				s, err := processString(s) // processes string type values
				if err == nil {
					Keys[key] = s
				}
				continue _lines
			}

			j := strings.Index(s, "/")
			if j != -1 {
				s = s[:j]
			}

			value := strings.TrimSpace(s)

			if value == "" { // we repeat this to take into account for empty values that have comments
				// we could not remove comments before processString because / is valid in a string value
				Keys[key] = nil
				continue
			}

			if (first >= '0' && first <= '9') || first == '+' || first == '-' {
				if strings.ContainsAny(value, ".DE") {
					value = strings.Replace(value, "D", "E", 1) // converts D type floats to E type
					x, _ := strconv.ParseFloat(value, 64)
					Keys[key] = x
				} else {
					x, _ := strconv.ParseInt(value, 10, 32)
					Keys[key] = int(x)
				}
			} else if first == 'T' {
				Keys[key] = true
			} else if first == 'F' {
				Keys[key] = false
			} else if first == '(' {
				var x, y float64
				fmt.Sscanf(value, "(%f,%f)", &x, &y)
				Keys[key] = complex(x, y)
			}
		}
		_, ends := Keys["END"]
		if ends {
			item, ok := Keys["NAXIS"]
			if ok {
				n := item.(int)
				h.Naxis = make([]int, n)
				for i := 0; i < n; i++ {
					h.Naxis[i] = Keys[Nth("NAXIS", i+1)].(int)
				}
			}

			break
		}
	}
	return h, nil
}
