#include "mesh2dh.hpp"

namespace acfd {


/** Reads Professor Luo's mesh file, which I call the 'domn' format.
   NOTE: Make sure nnofa is mentioned after ndim and ntype in the mesh file. ntype makes no sense for us now.
   Can only be used for linear mesh for now.
   MODIFIED FOR READING NON-HYBRID MESH FOR NOW.
*/
void UMesh2dh::readDomn(std::string mfile)
{
	std::ifstream infile(mfile);

	int nnode2, nfael2;
	
	// Do file handling here to populate npoin and nelem
	std::cout << "UMesh2dh: Reading domn mesh file...\n";
	char ch = '\0'; int dum = 0; double dummy;

	infile >> dum;
	infile >> ch;
	for(int i = 0; i < 4; i++)		//skip 4 lines
		do
			ch = infile.get();
		while(ch != '\n');
	infile >> ndim;
	infile >> nnode2;
	infile >> nfael2;
	infile >> maxnnofa;
	infile >> ch;			//get the newline
	do
		ch = infile.get();
	while(ch != '\n');
	infile >> nelem; infile >> npoin; infile >> nface;
	infile >> dummy; 				// get time
	ch = infile.get();			// clear newline

	nnode.resize(nelem,-1);
	nfael.resize(nelem,-1);

	nbtag = 2;
	ndtag = 2;

	//std::cout << "\nUTriMesh: Allocating coords..";
	coords.setup(npoin, ndim);
	// temporary array to hold connectivity matrix
	amat::Matrix<acfd_int > elms(nelem,nnode2);
	//std::cout << "UTriMesh: Allocating bface...\n";
	bface.setup(nface, maxnnofa + nbtag);
	
	//std::cout << "UTriMesh: Allocation done.";

	do
		ch = infile.get();
	while(ch != '\n');

	//now populate inpoel
	for(int i = 0; i < nelem; i++)
	{
		//infile >> nnode[i];
		infile >> dum;
		nnode[i] = nnode2;
		//nfael[i] = nnode[i];		// NOTE: assuming linear element
		nfael[i] = nfael2;

		for(int j = 0; j < nnode[i]; j++)
			infile >> elms(i,j);

		do
			ch = infile.get();
		while(ch != '\n');
	}
	std::cout << "UMesh2dh: Populated inpoel.\n";

	maxnnode = 3;
	for(int i = 0; i < nelem; i++)
		if(nnode[i] > maxnnode)
			maxnnode = nnode[i];
	
	inpoel.setup(nelem, maxnnode);

	for(int i = 0; i < nelem; i++)
		for(int j = 0; j < nnode[i]; j++)
			inpoel(i,j) = elms(i,j);
	
	//Correct inpoel:
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnode[i]; j++)
			inpoel(i,j)--;
	}

	ch = infile.get();
	do
		ch = infile.get();
	while(ch != '\n');

	// populate coords
	for(int i = 0; i < npoin; i++)
	{
		infile >> dum;
		for(int j = 0; j < ndim; j++)
			infile >> coords(i,j);
	}
	std::cout << "UMesh2dh: Populated coords.\n";

	// skip initial conditions
	ch = infile.get();
	for(int i = 0; i < npoin+2; i++)
	{
		do
			ch = infile.get();
		while(ch != '\n');
	}
	
	// populate bface
	for(int i = 0; i < nface; i++)
	{
		infile >> dum;
		for(int j = 0; j < ndim + nbtag; j++)
		{
			infile >> bface(i,j);
		}
		if (i==nface-1) break;
		do
			ch = infile.get();
		while(ch!='\n');
	}
	std::cout << "UMesh2dh: Populated bface. Done reading mesh.\n";
	//correct first 2 columns of bface
	for(int i = 0; i < nface; i++)
		for(int j = 0; j < 2; j++)
			bface(i,j)--;

	infile.close();

	vol_regions.setup(nelem, ndtag);
	vol_regions.zeros();
	
	std::cout << "UMesh2dh: Number of elements: " << nelem << ", number of points: " << npoin << ", max number of nodes per element: " << maxnnode << std::endl;
	std::cout << "Number of boundary faces: " << nface << ", Number of dimensions: " << ndim << std::endl;
	
	// set flag_bpoin
	flag_bpoin.setup(npoin,1);
	flag_bpoin.zeros();
	for(int i = 0; i < nface; i++)
		for(int j = 0; j < maxnnofa; j++)
			flag_bpoin(bface(i,j)) = 1;
}

/// Reads mesh from Gmsh 2 format file
void UMesh2dh::readGmsh2(std::string mfile, int dimensions)
{
	std::cout << "UMesh2d: readGmsh2(): Reading mesh file...\n";
	int dum; double dummy; std::string dums; char ch;
	ndim = dimensions;

	std::ifstream infile(mfile);
	for(int i = 0; i < 4; i++)		//skip 4 lines
		do
			ch = infile.get();
		while(ch != '\n');

	infile >> npoin;
	std::cout << "UMesh2d: readGmsh2(): No. of points = " << npoin << std::endl;
	coords.setup(npoin,ndim);

	// read coords of points
	for(int i = 0; i < npoin; i++)
	{
		infile >> dum;
		for(int j = 0; j < ndim; j++)
			infile >> coords(i,j);
		if(ndim < 3) infile >> dummy;
	}
	infile >> dums;		// get 'endnodes'
	infile >> dums;		// get 'elements'

	int width_elms = 25;
	int nelm, elmtype, nbtags, ntags;
	/// elmtype is the standard element type in the Gmsh 2 mesh format - of either faces or elements
	ndtag = 0; nbtag = 0;
	infile >> nelm;
	amat::Matrix<acfd_int> elms(nelm,width_elms);
	nface = 0; nelem = 0;
	std::vector<int> nnodes(nelm,0);
	std::vector<int> nnofas(nelm,0);
	std::vector<int> nfaels(nelm,0);
	std::vector<int> nintnodes(nelm,0);
	maxnnofa = 2;
	
	std::cout << "UMesh2d: readGmsh2(): Total number of elms is " << nelm << std::endl;

	for(int i = 0; i < nelm; i++)
	{
		infile >> dum;
		infile >> elmtype;

		/** elmtype is different for all faces and for all elements.
		 * nnodes contains number of nodes per face for the faces and number of nodes per element for elements.
		 */

		switch(elmtype)
		{
			case(1): // linear edge
				nnodes[i] = 2;
				if (nnodes[i] > maxnnofa) maxnnofa = nnodes[i];
				infile >> nbtags;
				if(nbtags > nbtag) nbtag = nbtags;
				for(int j = 0; j < nbtags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nface++;
				break;
			case(8): // quadratic edge
				nnodes[i] = 3;
				if (nnodes[i] > maxnnofa) maxnnofa = nnodes[i];
				infile >> nbtags;
				if(nbtags > nbtag) nbtag = nbtags;
				for(int j = 0; j < nbtags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nface++;
				break;
			case(2): // linear triangles
				nnodes[i] = 3;
				nfaels[i] = 3;
				nintnodes[i] = 0;
				if (maxnnofa < 2) maxnnofa = 2;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(3):	// linear quads
				nnodes[i] = 4;
				nfaels[i] = 4;
				nintnodes[i] = 0;
				if(maxnnofa < 2) maxnnofa = 2;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(9):	// quadratic triangles
				nnodes[i] = 6;
				nfaels[i] = 3;
				nintnodes[i] = 0;
				if(maxnnofa < 3) maxnnofa = 3;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(16):	// quadratic quad (8 nodes)
				nnodes[i] = 8;
				nfaels[i] = 4;
				nintnodes[i] = 0;
				if(maxnnofa < 3) maxnnofa = 3;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(10):	// quadratic quad (9 nodes)
				nnodes[i] = 9;
				nfaels[i] = 4;
				nintnodes[i] = 1;
				if(maxnnofa < 3) maxnnofa = 3;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			default:
				std::cout << "! UMesh2d: readGmsh2(): Element type not recognized. Setting as linear triangle." << std::endl;
				nnodes[i] = 3;
				nfaels[i] = 3;
				nintnodes[i] = 0;
				if(maxnnofa < 2) maxnnofa = 2;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnodes[i]);		// get tags
				for(int j = 0; j < nnodes[i]; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
		}
	}
	std::cout << "UMesh2d: readGmsh2(): Done reading elms" << std::endl;
	/*for(int i = 0; i < nelm; i++)
		std::cout << nnodes[i] << " " << nfaels[i] << std::endl;*/

	nnode.resize(nelem);
	nfael.resize(nelem);
	nintnodel.resize(nelem);
	nnobfa.resize(nface);

	// calculate max nnode and nfael
	maxnnode = nnodes[nface]; maxnfael = nfaels[nface];
	for(int i = 0; i < nelm; i++)
	{
		if(nnodes[i] > maxnnode)
			maxnnode = nnodes[i];
		if(nfaels[i] > maxnfael)
			maxnfael = nfaels[i];
	}

	if(nface > 0)
		bface.setup(nface, maxnnofa+nbtag);
	else std::cout << "UMesh2d: readGmsh2(): NOTE: There is no boundary data!" << std::endl;

	inpoel.setup(nelem, maxnnode);
	vol_regions.setup(nelem, ndtag);

	std::cout << "UMesh2dh: readGmsh2(): Done. No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << 
		",\n max number of nodes per element: " << maxnnode << ", max number of nodes per face: " << maxnnofa << ", max number of faces per element: " << maxnfael << std::endl;

	// write into inpoel and bface
	// the first nface rows to be read are boundary faces
	for(int i = 0; i < nface; i++)
	{
		nnobfa[i] = nnodes[i];
		for(int j = 0; j < nnobfa[i]; j++)
			bface(i,j) = elms(i,j)-1;			// -1 to correct for the fact that our numbering starts from zero
		for(int j = nnobfa[i]; j < nnobfa[i]+nbtag; j++)
			bface(i,j) = elms(i,j);
	}
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnodes[i+nface]; j++)
			inpoel(i,j) = elms(i+nface,j)-1;
		for(int j = 0; j < ndtag; j++)
			vol_regions(i,j) = elms(i+nface,j+nnodes[i+nface]);
		nnode[i] = nnodes[i+nface];
		nfael[i] = nfaels[i+nface];
		nintnodel[i] = nintnodes[i+nface];
	}
	infile.close();

	/*for(int i = 0; i < nelem; i++) 
	{	
		if (nnode == 3) nmtens[i] = 1;
		else if(nnode == 4) nmtens[i] = 4;
	}*/
	
	// set flag_bpoin
	flag_bpoin.setup(npoin,1);
	flag_bpoin.zeros();
	for(int i = 0; i < nface; i++)
		for(int j = 0; j < nnobfa[i]; j++)
			flag_bpoin(bface(i,j)) = 1;
}

void UMesh2dh::printmeshstats()
{
	std::cout << "UMesh2dh: No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << 
		",\n   max number of nodes per element: " << maxnnode << ", max number of nodes per face: " << maxnnofa << ", max number of faces per element: " << maxnfael << std::endl;
}

void UMesh2dh::writeGmsh2(std::string mfile)
{
	std::cout << "UMesh2dh: writeGmsh2(): writing mesh to file " << mfile << std::endl;
	// decide element type first, based on nfael/nnode and nnofa
	int elm_type = 2;

	std::ofstream outf(mfile);
	outf << std::setprecision(MESHDATA_DOUBLE_PRECISION);
	//std::cout << "nodes\n";
	outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
	outf << "$Nodes\n" << npoin << '\n';
	for(int ip = 0; ip < npoin; ip++)
	{
		outf << ip+1;
		for(int j = 0; j < ndim; j++)
			outf << " " << coords(ip,j);
		for(int j = 3-ndim; j > 0; j--)
			outf << " " << 0.0;
		outf << '\n';
	}
	outf << "$EndNodes\n";

	//std::cout << "boundary faces\n";
	outf << "$Elements\n" << nelem+nface << '\n';
	// boundary faces first
	for(int iface = 0; iface < nface; iface++)
	{
		int face_type = 1;
		if(nnobfa[iface] == 3)
			face_type = 8;

		outf << iface+1 << " " << face_type << " " << nbtag;
		for(int i = nnobfa[iface]; i < nnobfa[iface]+nbtag; i++)
			outf << " " << bface(iface,i);			// write tags
		for(int i = 0; i < nnobfa[iface]; i++)
			outf << " " << bface(iface,i)+1;		// write nodes
		outf << '\n';
	}
	//std::cout << "elements\n";
	for(int iel = 0; iel < nelem; iel++)
	{
		if(nnode[iel] == 3)
			elm_type = 2;
		else if(nnode[iel] == 4)
			elm_type = 3;
		else if(nnode[iel] == 6)
			elm_type = 9;
		else if(nnode[iel] == 8)
			elm_type = 16;
		else if(nnode[iel]==9)
			elm_type = 10;
		outf << nface+iel+1 << " " << elm_type << " " << ndtag;
		for(int i = 0; i < ndtag; i++)
			outf << " " << vol_regions(iel,i);
		for(int i = 0; i < nnode[iel]; i++)
			outf << " " << inpoel(iel,i)+1;
		outf << '\n';
	}
	outf << "$EndElements\n";

	outf.close();
}

/// \todo: TODO: There is an issue with psup for some boundary nodes belonging to elements of different types. Correct this.
void UMesh2dh::compute_topological()
{

	std::cout << "UMesh2dh: compute_topological(): Calculating and storing topological information...\n";
	/// 1. Elements surrounding points. Note that we only consider the vertices, not high-order points for this.
	//std::cout << "UMesh2d: compute_topological(): Elements surrounding points\n";
	esup_p.setup(npoin+1,1);
	esup_p.zeros();

	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nfael[i]; j++)
		{
			esup_p(inpoel(i,j)+1) += 1;		// inpoel(i,j) + 1 : the + 1 is there because an extra point in this element affects the starting point of the next.
		}
	}
	// Now make the members of esup_p cumulative
	for(int i = 1; i < npoin+1; i++)
		esup_p(i) += esup_p(i-1);
	// Now populate esup
	esup.setup(esup_p(npoin),1);
	esup.zeros();
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nfael[i]; j++)
		{
			int ipoin = inpoel(i,j);
			esup(esup_p(ipoin)) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
			esup_p(ipoin) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
		}
	}
	//But now esup_p holds increased values - each member increased by the number elements surrounding the corresponding point.
	// So correct this.
	for(int i = npoin; i >= 1; i--)
		esup_p(i,0) = esup_p(i-1,0);
	esup_p(0,0) = 0;
	// Elements surrounding points is now done.

	/// 2. Points surrounding points
	std::cout << "UMesh2dh: compute_topological(): Points surrounding points\n";
	psup_p.setup(npoin+1,1);
	psup_p.zeros();
	psup_p(0,0) = 0;
	amat::Matrix<int > lpoin(npoin,1);  // The ith member indicates the global point number of which the ith point is a surrounding point
	for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize this std::vector to -1
	int istor = 0;

	// first pass: calculate storage needed for psup
	for(int ip = 0; ip < npoin; ip++)
	{
		lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
		// Loop over elements surrounding this point
		for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
		{
			int ielem = esup(ie,0);		// element number

			// find local node number of ip in ielem
			int inode;
			for(int jnode = 0; jnode < nfael[ielem]; jnode++)
				if(inpoel(ielem,jnode) == ip) inode = jnode;

			std::vector<bool> nbd(nfael[ielem]);		// contains true if that local node number is connected to a particular local node.
			for(int j = 0; j < nfael[ielem]; j++)
				nbd[j] = false;

			if(nfael[ielem] == 3)
				for(int i = 0; i < nbd.size(); i++)
					nbd[i] = true;
			else if(nfael[ielem] == 4)
				for(int jnode = 0; jnode < nfael[ielem]; jnode++)
				{
					if(jnode == /*perm(0,nfael[ielem]-1,inode,1)*/(inode+1)%nfael[ielem] || jnode == (inode+nfael[ielem]-1)%nfael[ielem] /*perm(0,nfael[ielem]-1, inode, -1)*/)
						nbd[jnode] = true;
				}

			//loop over nodes of the element
			for(int inode = 0; inode < nfael[ielem]; inode++)
			{
				//Get global index of this node
				int jpoin = inpoel(ielem, inode);
				if(lpoin(jpoin,0) != ip && nbd[inode])		// test if this point as already been counted as a surrounding point of ip, and whether it's connected to ip.
				{
					istor++;
					lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
				}
			}
		}
		psup_p(ip+1,0) = istor;
	}

	psup.setup(istor,1);
	//std::cout << "+++ " << istor << std::endl;

	//second pass: populate psup
	istor = 0;
	for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize lpoin to -1
	for(int ip = 0; ip < npoin; ip++)
	{
		lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
		// Loop over elements surrounding this point
		for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
		{
			int ielem = esup(ie,0);		// element number

			// find local node number of ip in ielem
			int inode;
			for(int jnode = 0; jnode < nfael[ielem]; jnode++)
				if(inpoel(ielem,jnode) == ip) inode = jnode;

			std::vector<bool> nbd(nfael[ielem]);		// nbd[j] contains true if ip is connected to local node number j of ielem.
			for(int j = 0; j < nfael[ielem]; j++)
				nbd[j] = false;

			if(nfael[ielem] == 3)
				for(int i = 0; i < nbd.size(); i++)
					nbd[i] = true;
			else if(nfael[ielem] == 4)
				for(int jnode = 0; jnode < nfael[ielem]; jnode++)
				{
					if(jnode == /*perm(0,nfael[ielem]-1,inode,1)*/(inode+1)%nfael[ielem] || jnode == (inode+nfael[ielem]-1)%nfael[ielem] /*perm(0,nfael[ielem]-1, inode, -1)*/)
						nbd[jnode] = true;
				}

			//loop over nodes of the element
			for(int inode = 0; inode < nfael[ielem]; inode++)
			{
				//Get global index of this node
				int jpoin = inpoel(ielem, inode);
				if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
				{
					psup(istor,0) = jpoin;
					istor++;
					lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
				}
			}
		}
	}
	//Points surrounding points is now done.

	/// 3. Elements surrounding elements
	std::cout << "UMesh2dh: compute_topological(): Elements surrounding elements...\n";

	//amat::Matrix<int> lpoin(npoin,1);
	esuel.setup(nelem, maxnfael);
	for(int ii = 0; ii < nelem; ii++)
		for(int jj = 0; jj < maxnfael; jj++)
			esuel(ii,jj) = -1;
	//lpofa.mprint();
	const int lonnofa = 2;
	amat::Matrix<int > lhelp(lonnofa,1);
	lhelp.zeros();
	lpoin.zeros();

	for(int ielem = 0; ielem < nelem; ielem++)
	{
		// first get lpofa for this element
		amat::Matrix<int > lpofai(nfael[ielem], lonnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in [0:nnofa], i in [0:nfael])
		amat::Matrix<int > lpofaj;							// to be initialized for each jelem
		for(int i = 0; i < nfael[ielem]; i++)
		{
			for(int j = 0; j < lonnofa; j++)
			{
				lpofai(i,j) = (i+j) % nfael[ielem];		// fine as long as operands of % are not negative
			}
		}

		for(int ifael = 0; ifael < nfael[ielem]; ifael++)
		{
			for(int i = 0; i < lonnofa; i++)
			{
				lhelp(i,0) = inpoel(ielem, lpofai(ifael,i));	// lhelp stores global node nos. of current face of current element
				lpoin(lhelp(i,0)) = 1;
			}
			int ipoin = lhelp(0);
			for(int istor = esup_p(ipoin); istor < esup_p(ipoin+1); istor++)
			{
				int jelem = esup(istor);

				if(jelem != ielem)
				{
					// setup lpofa for jelem
					lpofaj.setup(nfael[jelem],lonnofa);
					for(int i = 0; i < nfael[jelem]; i++)
						for(int j = 0; j < lonnofa; j++)
							lpofaj(i,j) = (i+j)%nfael[jelem];

					for(int jfael = 0; jfael < nfael[jelem]; jfael++)
					{
						//Assume that no. of nodes in face ifael is same as that in face jfael
						int icoun = 0;
						for(int jnofa = 0; jnofa < lonnofa; jnofa++)
						{
							int jpoin = inpoel(jelem, lpofaj(jfael,jnofa));
							if(lpoin(jpoin)==1) icoun++;
						}
						if(icoun == lonnofa)		// lonnofa is 2
						{
							esuel(ielem,ifael) = jelem;
							esuel(jelem,jfael) = ielem;
						}
					}
				}
			}
			for(int i = 0; i < lonnofa; i++)
				lpoin(lhelp(i)) = 0;
		}
	}

	/** - Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac.
	 * Note that we assume that elements sharing a face have the same order, 
	 * which effectively means that all elements are assumed to have the same geometric p order.
	 * The orientation of the face is such that the element with smaller index is always to the left of the face, 
	 * while the element with greater index is always to the right of the face.
	 * - Computes nnofa for each face which is initially the same for all faces.
	 * - Also computes element-face connectivity array elemface in the same loop which computes intfac.
	 *
	 * \note After the following portion, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.
	 */

	std::cout << "UMesh2dh: compute_topological(): Computing intfac..." << std::endl;
	nbface = naface = 0;
	// first run: calculate nbface
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael[ie]; in++)
		{
			int je = esuel(ie,in);
			if(je == -1)
			{
				nbface++;
			}
		}
	}
	std::cout << "UMesh2dh: compute_topological(): Number of boundary faces = " << nbface << std::endl;
	// calculate number of internal faces
	naface = nbface;
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael[ie]; in++)
		{
			int je = esuel(ie,in);
			if(je > ie && je < nelem) naface++;
		}
	}
	std::cout << "UMesh2dh: compute_topological(): Number of all faces = " << naface << std::endl;

	// remember nnofa?
	nnofa.resize(naface);

	//allocate intfac and elemface
	intfac.setup(naface,maxnnofa+2);
	elemface.setup(nelem,maxnfael);

	//reset face totals
	nbface = naface = 0;

	int in1, je, jnode, jlocnode;

	//second run: populate intfac
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael[ie]; in++)
		{
			je = esuel(ie,in);
			if(je == -1)
			{
				esuel(ie,in) = nelem+nbface;
				intfac(nbface,0) = ie;
				intfac(nbface,1) = nelem+nbface;

				in1 = (in+1)%nfael[ie];
				intfac(nbface,2) = inpoel(ie,in);
				intfac(nbface,3) = inpoel(ie,in1);

				// high-order nodes
				int nhighperface = (nnode[ie]-nfael[ie]-nintnodel[ie])/nfael[ie];
				for(int i = 0; i < nhighperface; i++)
					intfac(nbface,i+4) = inpoel(ie, nfael[ie] + in*nhighperface + i);

				// finally set nnofa
				nnofa[nbface] = nhighperface + 2;

				elemface(ie,in) = nbface;

				nbface++;
			}
		}
	}
	naface = nbface;
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael[ie]; in++)
		{
			je = esuel(ie,in);
			if(je > ie && je < nelem)
			{
				in1 = (in+1)%nfael[ie];
				intfac(naface,0) = ie;
				intfac(naface,1) = je;
				intfac(naface,2) = inpoel(ie,in);
				intfac(naface,3) = inpoel(ie,in1);

				// high-order nodes
				int nhighperface = (nnode[ie]-nfael[ie]-nintnodel[ie])/nfael[ie];
				for(int i = 0; i < nhighperface; i++)
					intfac(naface,i+4) = inpoel(ie, nfael[ie] + in*nhighperface + i);

				// finally set nnofa
				nnofa[naface] = nhighperface + 2;
				
				elemface(ie,in) = naface;
				for(jnode = 0; jnode < nfael[je]; jnode++)
					if(inpoel(ie,in1) == inpoel(je,jnode))
						elemface(je,jnode) = naface;

				naface++;
			}
		}
	}

	/// Finally, calculates bpoints.

	//first get number of bpoints
	nbpoin = 0;
	amat::Matrix<int > isbpflag(npoin,1);
	isbpflag.zeros();
	for(int i = 0; i < nface; i++)
	{
		for(int j = 0; j < nnobfa[i]; j++)
			isbpflag(bface(i,j)) = 1;
	}
	for(int i = 0; i < npoin; i++)
		if(isbpflag(i)==1) nbpoin++;

	std::cout << "UMesh2dh: compute_topological(): Number of boundary points = " << nbpoin << std::endl;
	std::cout << "UMesh2dh: compute_topological(): Done." << std::endl;
}

void UMesh2dh::compute_boundary_maps()
{
	const int lonnofa = 2;
	
	// iterate over bfaces and find corresponding intfac face for each bface
	bifmap.setup(nbface,1);
	ifbmap.setup(nbface,1);

	std::vector<int> fpo(lonnofa);

	for(int ibface = 0; ibface < nface; ibface++)
	{
		for(int i = 0; i < lonnofa; i++)
			fpo[i] = bface(ibface,i);

		int inface = -1;

		// iterate over intfacs - check if bface ibface matches intfac iface, for each iface
		for(int iface = 0; iface < nbface; iface++)
		{
			bool final1 = true;

			std::vector<bool> inter(lonnofa);
			for(int b = 0; b < lonnofa; b++)
				inter[b] = false;						// initially set all bools to false

			for(int j = 0; j < lonnofa; j++)
			{
				for(int k = 0; k < lonnofa; k++)
					if(fpo[j] == intfac(iface, 2+k)) {
						inter[j] = true;			// if jth node of ibface has a node of iface, it belongs to iface; set the corresp. boolean to true
						break;
					}
			}

			/*for(int i = 0; i < lonnofa; i++)
				std::cout << inter[i];
			std::cout << std::endl;*/

			for(int b = 0; b < lonnofa; b++)
				if(inter[b] == false) final1 = false;						// if any node of ibface failed to find a node of iface, ibface is not the same as iface

			if(final1 == true) inface = iface;
		}

		if(inface != -1) {
			bifmap(inface) = ibface;
			ifbmap(ibface) = inface;
		}
		else {
			std::cout << "UMesh2d: compute_boundary_maps(): ! intfac face corresponding to " << ibface << "th bface not found!!" << std::endl;
			continue;
		}
	}
	isBoundaryMaps = true;
}

void UMesh2dh::writeBoundaryMapsToFile(std::string mapfile)
{
	if(isBoundaryMaps == false) {
		std::cout << "UMesh2d: writeBoundaryMapsToFile(): ! Boundary maps not available!" << std::endl;
		return;
	}
	std::ofstream ofile(mapfile);
	ofile << nbface << '\n'<< "bifmap\n";
	for(int i = 0; i < nbface; i++)
		ofile << bifmap.get(i) << ' ';
	ofile << '\n';
	ofile << "ifbmap\n";
	for(int i = 0; i < nbface; i++)
		ofile << ifbmap.get(i) << ' ';
	ofile << '\n';
	ofile.close();
}

void UMesh2dh::readBoundaryMapsFromFile(std::string mapfile)
{
	std::ifstream ofile(mapfile);
	std::string dum; int sz;
	ofile >> sz >>  dum;
	std::cout << "UMesh2d: readBoundaryMapsFromFile(): Number of boundary faces in file = " << sz << std::endl;
	bifmap.setup(sz,1);
	ifbmap.setup(sz,1);

	for(int i = 0; i < nbface; i++)
		ofile >> bifmap(i);

	ofile >> dum;
	for(int i = 0; i < nbface; i++)
		ofile >> ifbmap(i);

	ofile.close();
	isBoundaryMaps = true;
}

void UMesh2dh::compute_intfacbtags()
{
	/// Populate intfacbtags with boundary markers of corresponding bfaces

	intfacbtags.setup(nbpoin,nbtag);

	if(isBoundaryMaps == false)
	{
		std::cout << "UMesh2d: compute_intfacbtags(): ! Boundary maps are not available!" << std::endl;
		return;
	}

	for(int ibface = 0; ibface < nface; ibface++)
	{
		for(int j = 0; j < nbtag; j++)
			intfacbtags(ifbmap(ibface),j) = bface(ibface,nnobfa[ibface]+j);
	}
}


} // end namespace
