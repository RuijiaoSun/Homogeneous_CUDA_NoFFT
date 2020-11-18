Order of "output.lyr" parameters.

The wavenumber in free space: 		 	k0	double			8B
The direction of propogation:   	 	d	double*2		16B
The refractive index in the first layer: 	n[0]	complex<double>		16B

for i in LAYERS:
	z positions[i]:			 	z[i]	double			8B * LAYERS
	z-component of propogation directions:	sz[i]	complex<double>		16B * LAYERS
	Transmission:				Ptx[i]	complex<double>		16B * LAYERS
	Transmission:				Pty[i]	complex<double>		16B * LAYERS
	Transmission:				Ptz[i]	complex<double>		16B * LAYERS
	Reflection:				Prx[i]	complex<double>		16B * LAYERS
	Reflection:				Pry[i]	complex<double>		16B * LAYERS
	Reflection:				Prz[i]	complex<double>		16B * LAYERS
	

All parameters we need will be:
	15 * LAYERS + 5
