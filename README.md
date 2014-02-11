Copyright 2014 Shahriar Iravanian (siravan@svtsim.com). All rights reserved. Use of this source code is governed by a MIT license that can be found in the LICENSE file.

Package fits reads and processes FITS files. It is written in pure golang and is not a wrapper around another library or a direct translation of another library to golang. The main purpose is to provide a native golang solution to reading FITS file and to assess the suitability of golang for scientific and numerical applications.

FITS is a common format for astronomical image and data. This package is based on version 3.0 of the FITS standard:

Pence W.D., Chiappetti L., Page C. G., Shaw R. A., Stobie E. Definition of the Flexible Image Transport System (FITS), version 3.0. A&A 524, A42 (2010)
http://www.aanda.org/articles/aa/abs/2010/16/aa15362-10/aa15362-10.html

The following features are supported in the current version:

1. Images with all six different data format (byte, int16, int32, int64, float32, and float64)
2. Text and binary tables with atomic and fixed-size array elements

The following features are not yet implemented:

1. Automatic application of BSCALE/BZERO
2. Random group structure
3. Variable length arrays in binary tables
4. World coordinate system

Also note that currently this package provides only read capability and does not write/generate a FITS file.

The basic usage of the package is by calling Open function. It accepts a reader that should provide a valid FITS file. The output is a []*fits.Unit, where Unit represents a Header/Data Unit (i.e. a header with the corresponding data). Unit provides a set of variables and functions to access the HDU data.

Let 'test.fits' be a FITS file with two HDU. The first one is of type SIMPLE and contains a single two-dimensional image with the following parameters:

BITPIX  =  -32
NAXIS   =  2
NAXIS1  =  512
NAXIS2  =  256

The second HDU contains a binary table (XTENSION=BINTABLE):

BITPIX  =  8
NAXIS   =  2
NAXIS1  =  100
NASIX2  =  5
TFIELDS =  10
TFORM1  =  E
TTYPE   =  FLUX
TDISP1  =  F10.4

To read this file, we first call

units := fits.Open("test.fits")

Now, units[0] points to the first HDU. We can access the header keys by using units.Keys map. For example, units[0].Keys["BITPIX"].(int) returns -32. Note that Keys stores interface{} and appropriate type-assertion needs to be done. Unit.Naxis returns a slice of integers ([]int) containing all NAXIS data. For example, units[0].Naxis is equal to [512, 256]. We can access the image data points by using one of the three accessor functions: Unit.At, Unit.IntAt and Unit.FloatAt. Each function accepts NAXIS integer arguments and returns the pixel value at that location. Unit.At returns an interface{} and needs to be type-asserted before use. Unit.IntAt and Unit.FloatAt return int64 and float64, respectively.

For table data, we use two other accessor functions: Field and Format. Field accepts one argument, col, that define a field. It can be 0-based int or a string. For example, units[1].Field(0) and units[1].Field("FLUX") both points to the same column. The return value of Field is another function, which is the actual accessor function and accepts one int argument representing a row. For example, units[1].Field("Flux")(1) returns the value of column "FLUX" in the second row of the table as interface{}. The following code populates a slice of float with the value of the FLUX column:

fn := units[1].Field("FLUX")
x := make([]float32, units[1].Naxis[1])     // note Naxis[0]=NAXIS1=length of a row, Naxis[1]=NAXIS2=number of rows
for row := range x {
    x[row] = fn(row).(float32)
}

Format function on the hand accepts two arguments, col (same as Field) and row and return a string formatted according to TDISP for the field. For example, if units[1].Field("Flux")(1) is equal to 987.654321, then units[1].Format("Flux", 1) returns "987.6543".

You can find the documentation for fits.go in http://godoc.org/github.com/siravan/fits.

demo/extract.go is a test program which show cases typical use of the package.
