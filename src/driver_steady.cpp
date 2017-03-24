
	// set inviscid flux scheme
	if(invflux == "VANLEER")
		inviflux = new VanLeerFlux(NVARS, m->gndim(), g);
	else if(invflux == "ROE")
	{
		inviflux = new RoeFlux(NVARS, m->gndim(), g);
		std::cout << "SpatialBase: Using Roe fluxes." << std::endl;
	}
	else if(invflux == "HLLC")
	{
		inviflux = new HLLCFlux(NVARS, m->gndim(), g);
		std::cout << "SpatialBase: Using HLLC fluxes." << std::endl;
	}
	else
		std::cout << "SpatialBase: ! Flux scheme not available!" << std::endl;

	// set reconstruction scheme
	std::cout << "SpatialBase: Reconstruction scheme is " << reconst << std::endl;
	if(reconst == "GREENGAUSS")
	{
		rec = new GreenGaussReconstruction();
		//rec->setup(m, &u, &ug, &dudx, &dudy, &rc, &rcg);
	}
	else
	{
		rec = new WeightedLeastSquaresReconstruction();
	}
	if(order == 1) std::cout << "SpatialBase: No reconstruction" << std::endl;

	// set limiter
	if(limiter == "NONE")
	{
		lim = new NoLimiter(m, &u, &ug, &dudx, &dudy, &rcg, &rc, gr, &uleft, &uright);
		std::cout << "SpatialBase: No limiter will be used." << std::endl;
	}
	else if(limiter == "WENO")
	{
		lim = new WENOLimiter(m, &u, &ug, &dudx, &dudy, &rcg, &rc, gr, &uleft, &uright);
		std::cout << "SpatialBase: WENO limiter selected.\n";
	}
