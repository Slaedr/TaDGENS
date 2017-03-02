/** @file mesh2dh.hpp
 * @brief Contains a class to handle 2D hybrid meshes containing triangles and quadrangles
 * @author Aditya Kashi
 */

#ifndef __AMESH2DH_H
#define __AMESH2DH_H

#ifndef __ACONSTANTS_H
#include "aconstants.hpp"
#endif

#ifndef __ARRAY2D_H
#include "aarray2d.hpp"
#endif

namespace acfd {

/// General hybrid unstructured mesh class supporting triangular and quadrangular elements
class UMesh2dh
{
private:
	int ndim;								///< Number of space dimensions of mesh - can only be 2, though
	acfd_int npoin;							///< Number of nodes
	acfd_int nelem;							///< Number of elements
	acfd_int nface;							///< Number of boundary faces
	std::vector<int> nnode;					///< number of nodes to an element, for each element
	std::vector<int> nintnodel;				///< Number of internal (not lying on a face) nodes in an element
	int maxnnode;							///< Maximum number of nodes per element for any element
	int maxnnofa;							///< Maximum number of nodes per face for any face
	std::vector<int> nfael;					///< number of faces to an element (equal to number of edges to an element in 2D) for each element
	int maxnfael;							///< Maximum number of faces per element for any element
	std::vector<int> nnofa;					///< number of nodes in a face
	std::vector<int> nnobfa;				///< number of nodes in boundary faces
	acfd_int naface;						///< total number of (internal and boundary) faces
	acfd_int nbface;						///< number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	acfd_int nbpoin;						///< number of boundary points
	int nbtag;								///< number of tags for each boundary face
	int ndtag;								///< number of tags for each element
	amat::Array2d<acfd_real> coords;		///< Specifies coordinates of each node
	amat::Array2d<acfd_int> inpoel;			///< Interconnectivity matrix: lists node numbers of nodes in each element
	amat::Array2d<acfd_int> bface;			///< Boundary face data: lists nodes belonging to a boundary face and contains boudnary markers
	amat::Array2d<int> vol_regions;			///< to hold volume region markers, if any
	amat::Array2d<int> flag_bpoin;			///< Holds 1 or 0 for each point depending on whether or not that point is a boundary point

	/// List of indices of [esup](@ref esup) corresponding to vertices (vertices = "low order" nodes only)
	amat::Array2d<acfd_int> esup_p;
	/// List of elements surrounding vertices
	/** Integers pointing to particular vertices' element lists are stored in [esup_p](@ref esup_p). */
	amat::Array2d<acfd_int> esup;

	/// Lists of indices of psup corresponding to vertices
	amat::Array2d<acfd_int> psup_p;
	/// List of vertices surrounding vertices
	/** Integers pointing to particular vertices' vertex lists are stored in [psup_p](@ref psup_p)
	 */
	amat::Array2d<acfd_int> psup;
	
	/// Elements surrounding elements
	amat::Array2d<acfd_int> esuel;
	/// Face data structure - contains info about elements and nodes associated with a face
	amat::Array2d<acfd_int> intfac;
	/// Holds boundary tags (markers) corresponding to intfac
	amat::Array2d<int> intfacbtags;
	/// Holds face numbers of faces making up an element
	amat::Array2d<acfd_int> elemface;

	amat::Array2d<acfd_int> bifmap;				///< relates boundary faces in intfac with bface, ie, bifmap(intfac no.) = bface no.
	amat::Array2d<acfd_int> ifbmap;				///< relates boundary faces in bface with intfac, ie, ifbmap(bface no.) = intfac no.
	bool isBoundaryMaps;			///< Specifies whether bface-intfac maps have been created

public:
		
	/* Functions to get mesh data. */

	double gcoords(int pointno, int dim) const
	{
		return coords.get(pointno,dim);
	}
	int ginpoel(int elemno, int locnode) const
	{
		return inpoel.get(elemno, locnode);
	}
	int gbface(int faceno, int val) const
	{
		return bface.get(faceno, val);
	}

	amat::Array2d<double >* getcoords()
	{ return &coords; }

	acfd_int gesup(acfd_int i) const { return esup(i); }
	acfd_int gesup_p(acfd_int i) const { return esup_p(i); }
	acfd_int gpsup(acfd_int i) const { return psup(i); }
	acfd_int gpsup_p(acfd_int i) const { return psup_p(i); }
	acfd_int gesuel(acfd_int ielem, int jnode) const { return esuel(ielem, jnode); }
	acfd_int gelemface(acfd_int ielem, int inode) const { return elemface(ielem,inode); }
	acfd_int gintfac(acfd_int face, int i) const { return intfac(face,i); }
	int gintfacbtags(acfd_int face, int i) const { return intfacbtags(face,i); }
	acfd_int gbifmap(acfd_int intfacno) const { return bifmap(intfacno); }
	acfd_int gifbmap(acfd_int bfaceno) const { return ifbmap(bfaceno); }
	int gflag_bpoin(const acfd_int pointno) const { return flag_bpoin(pointno); }

	acfd_int gnpoin() const { return npoin; }
	acfd_int gnelem() const { return nelem; }
	acfd_int gnface() const { return nface; }
	acfd_int gnbface() const { return nbface; }
	int gnnode(acfd_int ielem) const { return nnode[ielem]; }
	acfd_int gnaface() const {return naface; }
	int gnfael(acfd_int ielem) const { return nfael[ielem]; }
	int gnnofa(acfd_int iface) const { return nnofa[iface]; }
	int gnnobfa(acfd_int ibface) const { return nnobfa[ibface]; }
	int gmaxnnofa() const { return maxnnofa; }
	int gnbtag() const{ return nbtag; }
	int gndtag() const { return ndtag; }
	acfd_int gnbpoin() const { return nbpoin; }

	/* Functions to set some mesh data structures. */
	/// set coordinates of a certain point; 'set' counterpart of the 'get' function [gcoords](@ref gcoords).
	void scoords(const acfd_int pointno, const int dim, const acfd_real value)
	{
		coords(pointno,dim) = value;
	}

	void setcoords(amat::Array2d<double >* c)
	{ coords = *c; }

	void setinpoel(amat::Array2d<int >* inp)
	{ inpoel = *inp; }

	void setbface(amat::Array2d<int >* bf)
	{ bface = *bf; }

	void modify_bface_marker(int iface, int pos, int number)
	{ bface(iface, pos) = number; }

	/** \brief Reads Professor Luo's mesh file, which I call the 'domn' format.
	 * 
	 * \note NOTE: Make sure nfael and nnofa are mentioned after ndim and nnode in the mesh file.
	*/
	void readDomn(std::string mfile);

	/// Reads mesh from Gmsh 2 format file
	void readGmsh2(std::string mfile, int dimensions);
	
	/// Stores (in array bpointsb) for each boundary point: the associated global point number and the two bfaces associated with it.
	void compute_boundary_points();

	void printmeshstats();
	void writeGmsh2(std::string mfile);

	/** Computes data structures for 
	 * elements surrounding point (esup), 
	 * points surrounding point (psup), 
	 * elements surrounding elements (esuel), 
	 * elements surrounding faces along with points in faces (intfac),
	 * element-face connectivity array elemface (for each facet of each element, it stores the intfac face number)
	 * a list of boundary points with correspong global point numbers and containing boundary faces (according to intfac) (bpoints).
	 * \note
	 * - Use only after setup()
	 * - Currently only works for linear mesh
	 */
	void compute_topological();

	/// Iterates over bfaces and finds the corresponding intfac face for each bface
	/** Stores this data in the boundary label maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap).
	 */
	void compute_boundary_maps();
	
	/// Writes the boundary point maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap) to a file
	void writeBoundaryMapsToFile(std::string mapfile);
	/// Reads the boundary point maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap) from a file
	void readBoundaryMapsFromFile(std::string mapfile);
	
	/// Populate [intfacbtags](@ref intfacbtags) with boundary markers of corresponding bfaces
	void compute_intfacbtags();
};


} // end namespace
#endif
