Copyright 2014 Shahriar Iravanian (siravan@svtsim.com).  All rights reserved.
Use of this source code is governed by a MIT license that can be found in the LICENSE file.

<p>
Package fits reads and processes FITS files. It is written in pure golang and is not a wrapper around another library or a direct translation of
another library to golang. The main purpose is to provide a native golang solution to reading FITS file and to assess the suitability of golang for
scientific and numerical applications.
</p>
<p>
FITS is a common open-source format for storage and transmission of astronomical images and data.
This package is based on <a href="http://www.aanda.org/articles/aa/abs/2010/16/aa15362-10/aa15362-10.html">version 3.0 of the FITS standard</a>.
</p>
<p>
The following features are supported in the current version:
</p>
<pre>1. Images with all six different data format (byte, int16, int32, int64, float32, and float64)
2. Text and binary tables with atomic and fixed-size array elements
</pre>
<p>
The following features are not yet implemented:
</p>
<pre>1. Automatic application of BSCALE/BZERO
2. Random group structure
3. Variable length arrays in binary tables
4. World coordinate system
</pre>
<p>
Also note that currently this package provides only read capability and does not write/generate a FITS file.
</p>
<p>
The basic usage of the package is by calling Open function. It accepts a reader that should provide a valid FITS file.
The output is a []*fits.Unit, where Unit represents a Header/Data Unit (i.e. a header with the corresponding data).
Unit provides a set of variables and functions to access the HDU data.
</p>
<p>
Let &#39;test.fits&#39; be a FITS file with two HDU. The first one is of type SIMPLE and contains a single two-dimensional image with the following parameters:
</p>
<pre>BITPIX  =  -32
NAXIS   =  2
NAXIS1  =  512
NAXIS2  =  256
</pre>
<p>
The second HDU contains a binary table (XTENSION=BINTABLE):
</p>
<pre>BITPIX  =  8
NAXIS   =  2
NAXIS1  =  100
NASIX2  =  5
TFIELDS =  10
TFORM1  =  E
TTYPE   =  FLUX
TDISP1  =  F10.4
</pre>
<p>
To read this file, we first call
</p>
<pre>units := fits.Open(&#34;test.fits&#34;)
</pre>
<p>
Now, units[0] points to the first HDU. We can access the header keys by using units.Keys map.
For example, units[0].Keys[&#34;BITPIX&#34;].(int) returns -32. Note that Keys stores interface{} and appropriate type-assertion needs to be done.
Unit.Naxis returns a slice of integers ([]int) containing all NAXIS data. For example, units[0].Naxis is equal to [512, 256].
We can access the image data points by using one of the three accessor functions: Unit.At, Unit.IntAt and Unit.FloatAt.
Each function accepts NAXIS integer arguments and returns the pixel value at that location.
Unit.At returns an interface{} and needs to be type-asserted before use. Unit.IntAt and Unit.FloatAt return int64 and float64, respectively.
</p>
<p>
For table data, we use two other accessor functions: Field and Format.
Field accepts one argument, col, that define a field. It can be 0-based int or a string.
For example, units[1].Field(0) and units[1].Field(&#34;FLUX&#34;) both points to the same column.
The return value of Field is another function, which is the actual accessor function and accepts one int argument representing a row.
For example, units[1].Field(&#34;Flux&#34;)(1) returns the value of column &#34;FLUX&#34; in the second row of the table as interface{}.
The following code populates a slice of float with the value of the FLUX column (note that Naxis[0] is the number of bytes in a row and Naxis[1] is the number of rows):
</p>
<pre>fn := units[1].Field(&#34;FLUX&#34;)
x := make([]float32, units[1].Naxis[1])     
for row := range x {
    x[row] = fn(row).(float32)
}
</pre>
<p>
Format function on the hand accepts two arguments, col (same as Field) and row and return a string formatted according to TDISP for the field.
For example, if units[1].Field(&#34;Flux&#34;)(1) is equal to 987.654321, then units[1].Format(&#34;Flux&#34;, 1) returns &#34;987.6543&#34;.</p>

