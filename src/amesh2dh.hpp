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
	int g_degree;							///< Degree of (polynomial) geometric mapping from reference element to physical element
	a_int npoin;							///< Number of nodes
	a_int nelem;							///< Number of elements
	a_int nface;							///< Number of boundary faces
	std::vector<int> nnode;					///< number of nodes to an element, for each element
	std::vector<int> nintnodel;				///< Number of internal (not lying on a face) nodes in an element
	int maxnnode;							///< Maximum number of nodes per element for any element
	int maxnnofa;							///< Maximum number of nodes per face for any face
	std::vector<int> nfael;					///< number of faces to an element (equal to number of edges to an element in 2D) for each element
	int maxnfael;							///< Maximum number of faces per element for any element
	std::vector<int> nnofa;					///< number of nodes in a face
	std::vector<int> nnobfa;				///< number of nodes in boundary faces
	a_int naface;							///< total number of (internal and boundary) faces
	a_int nbface;							///< number of boundary faces as calculated, as opposed to nface which is read from file
	a_int nbpoin;							///< number of boundary points
	int nbtag;								///< number of tags for each boundary face
	int ndtag;								///< number of tags for each element
	amat::Array2d<a_real> coords;			///< Specifies coordinates of each node
	amat::Array2d<a_int> inpoel;			///< Interconnectivity matrix: lists node numbers of nodes in each element
	amat::Array2d<a_int> bface;				///< Boundary face data: lists nodes belonging to a boundary face and contains boudnary markers
	amat::Array2d<int> vol_regions;			///< to hold volume region markers, if any
	amat::Array2d<int> flag_bpoin;			///< Holds 1 or 0 for each point depending on whether or not that point is a boundary point

	/// List of indices of [esup](@ref esup) corresponding to vertices (vertices = "low order" nodes only)
	amat::Array2d<a_int> esup_p;
	/// List of elements surrounding vertices
	/** Integers pointing to particular vertices' element lists are stored in [esup_p](@ref esup_p). */
	amat::Array2d<a_int> esup;

	/// Lists of indices of psup corresponding to vertices
	amat::Array2d<a_int> psup_p;
	/// List of vertices surrounding vertices
	/** Integers pointing to particular vertices' vertex lists are stored in [psup_p](@ref psup_p)
	 */
	amat::Array2d<a_int> psup;
	
	/// Elements surrounding elements
	amat::Array2d<a_int> esuel;
	
	/// Face data structure - contains info about elements and nodes associated with a face
	/** For each this contains: left element index, right element index, 
	 * "starting" node index, "ending" node index and inner node indices, in that order.
	 */
	amat::Array2d<a_int> intfac;

	/// Holds boundary tags (markers) corresponding to intfac
	amat::Array2d<int> intfacbtags;

	/// Local face numbers, in the left element and right element, for each intfac face
	/** facelocalnum(iface,0) holds local face number of face iface as seen by the left element,
	 * facelocalnum(iface,1) holds the local face number of iface as seen by the right element.
	 */
	amat::Array2d<a_int> facelocalnum;

	/// Holds face numbers of faces making up an element
	amat::Array2d<a_int> elemface;

	amat::Array2d<a_int> bifmap;				///< relates boundary faces in intfac with bface, ie, bifmap(intfac no.) = bface no.
	amat::Array2d<a_int> ifbmap;				///< relates boundary faces in bface with intfac, ie, ifbmap(bface no.) = intfac no.
	bool isBoundaryMaps;						///< Specifies whether bface-intfac maps have been created

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

	a_int gesup(a_int i) const { return esup(i); }
	a_int gesup_p(a_int i) const { return esup_p(i); }
	a_int gpsup(a_int i) const { return psup(i); }
	a_int gpsup_p(a_int i) const { return psup_p(i); }
	a_int gesuel(a_int ielem, int jnode) const { return esuel(ielem, jnode); }
	a_int gelemface(a_int ielem, int inode) const { return elemface(ielem,inode); }
	a_int gintfac(a_int face, int i) const { return intfac(face,i); }
	int gintfacbtags(a_int face, int i) const { return intfacbtags(face,i); }
	int gfacelocalnum(a_int face, int leftright) const { return facelocalnum(face, leftright); }
	a_int gbifmap(a_int intfacno) const { return bifmap(intfacno); }
	a_int gifbmap(a_int bfaceno) const { return ifbmap(bfaceno); }
	int gflag_bpoin(const a_int pointno) const { return flag_bpoin(pointno); }
	int degree() const { return g_degree; }

	a_int gnpoin() const { return npoin; }
	a_int gnelem() const { return nelem; }
	a_int gnface() const { return nface; }
	a_int gnbface() const { return nbface; }
	int gnnode(a_int ielem) const { return nnode[ielem]; }
	a_int gnaface() const {return naface; }
	int gnfael(a_int ielem) const { return nfael[ielem]; }
	int gnnofa(a_int iface) const { return nnofa[iface]; }
	int gnnobfa(a_int ibface) const { return nnobfa[ibface]; }
	int gmaxnnofa() const { return maxnnofa; }
	int gnbtag() const{ return nbtag; }
	int gndtag() const { return ndtag; }
	a_int gnbpoin() const { return nbpoin; }

	/* Functions to set some mesh data structures. */
	/// set coordinates of a certain point; 'set' counterpart of the 'get' function [gcoords](@ref gcoords).
	void scoords(const a_int pointno, const int dim, const a_real value)
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
	/** Also sets a [flag](@ref flag_bpoin) for each point according to whether or not it is a boundary point.
	 */
	void readGmsh2(std::string mfile, int dimensions);
	
	/// Stores (in array bpointsb) for each boundary point: the associated global point number and the two bfaces associated with it.
	void compute_boundary_points();

	void printmeshstats();
	void writeGmsh2(std::string mfile);

	/** Computes data structures for 
	 * elements surrounding point (esup), 
	 * points surrounding point (psup), 
	 * elements surrounding elements (esuel), 
	 * elements surrounding faces along with points in faces, [see](@ref intfac),
	 * element-face connectivity array elemface (for each facet of each element, it stores the intfac face number)
	 * local face numbers in left and right elements for each face, ie, the [face local number](@ref facelocalnum)
	 * a list of boundary points with correspong global point numbers and containing boundary faces (according to intfac) [see](@ref bpoints).
	 * \note
	 * - Use only after setup()
	 * - Currently only works for linear mesh
	 */
	void compute_topological();

	/// Iterates over bfaces and finds the corresponding intfac face for each bface
	/** Stores this data in the boundary label maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap).
	 * Also stores boundary markers in [intfacbtags](@ref intfacbtags).
	 */
	void compute_boundary_maps();
	
	/// Writes the boundary point maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap) to a file
	void writeBoundaryMapsToFile(std::string mapfile);
	/// Reads the boundary point maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap) from a file
	void readBoundaryMapsFromFile(std::string mapfile);
};


} // end namespace
#endif
