/** @file mesh2dh.hpp
 * @brief Contains a class to handle 2D hybrid meshes containing triangles and quadrangles
 * @author Aditya Kashi
 */

#ifndef AMESH2DH_H
#define AMESH2DH_H

#include <utility>
#include "aconstants.hpp"
#include "utilities/aarray2d.hpp"

namespace acfd {

/// Index of something w.r.t. the element it is associated with
typedef int EIndex;
/// Index of something (usually a node) w.r.t. the face it is associated with
typedef int FIndex;

/// General hybrid unstructured mesh class supporting triangular and quadrangular elements
class UMesh2dh
{
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
	a_real gedgelengthsquared(const a_int iface) const { return els[iface]; }
	a_real gelemdiam(const a_int iel) const { return eldiam[iel]; }

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

	/// Reads mesh from Gmsh 2 format file
	/** Also sets a [flag](@ref flag_bpoin) for each point according to whether or not it is a boundary point.
	 */
	void readGmsh2(std::string mfile, int dimensions);
	
	/// Stores (in array bpointsb) for each boundary point: the associated global point number
	///  and the two bfaces associated with it.
	void compute_boundary_points();

	void printmeshstats();
	void writeGmsh2(std::string mfile);

	/// Checks whether a bface and the corresponding element face have the same orientation
	void correctBoundaryFaceOrientation();

	/** Computes data structures for 
	 * elements surrounding point (esup),
	 * elements surrounding elements (esuel),
	 * elements surrounding faces along with points in faces, [see](@ref intfac),
	 * element-face connectivity array elemface
	 *   (for each facet of each element, it stores the intfac face number)
	 * local face numbers in left and right elements for each face, ie,
	 *   the [face local number](@ref facelocalnum)
	 * Also checks for consistency of orientation of boundary faces - see correctBoundaryFaceOrientation
	 */
	void compute_topological();

	/// Iterates over bfaces and finds the corresponding intfac face for each bface
	/** Stores this data in the boundary label maps [ifbmap](@ref ifbmap) and [bifmap](@ref bifmap).
	 * Also stores boundary markers in [intfacbtags](@ref intfacbtags).
	 */
	void compute_boundary_maps();

	/** Computes edge lengths and element diameters.
	 * Needs to be called once before any call to meshSizeParameter
	 */
	void compute_edge_elem_sizes();

	/// Computes the "mesh size" h
	/** Call only after compute_topological() has been called.
	 */
	a_real meshSizeParameter() const;

	/// Get the index of a vertex w.r.t. an element (ie., get the vertex's "EIndex") from
	///  its index in a face of that element
	/** \param ielem Element index in the subdomain
	 * \param faceEIndex Index of a face in the element w.r.t. the element (ie., the face's EIndex)
	 * \param nodeFIndex Index of a node in the face above w.r.t. the face (ie., the node's FIndex)
	 *
	 * The current implementation works only in 2D, and the only for vertices ("low-order nodes").
	 */
	EIndex getNodeEIndex(const a_int ielem, const EIndex iface, const FIndex inode) const
	{
		static_assert(NDIM==2, "Only 2D meshes are currently supported!");
		assert(ielem < nelem); assert(ielem >= 0);
		return (iface + inode) % nnode[ielem];
	}

	/// Returns the EIndex of a face in a certain element
	/** Returns negative if the face is not present in that element.
	 * \warning If iface is a physical boundary face, this function will always work. But if iface is
	 *  an interior face defined according to \ref intfac, obviously intfac must be available.
	 * \param[in] phyboundary If true, iface is interpreted as a bface index between 0 and nbface.
	 *   If false, iface is interpreted as an intfac index.
	 */
	EIndex getFaceEIndex(const bool phyboundary, const a_int iface, const a_int elem) const;

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
	std::vector<a_real> els;				///< Holds squared (approx for high-order) length of each edge of the mesh
	std::vector<a_real> eldiam;				///< Diameter of each element

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

	/// Holds intfac face numbers of faces making up an element
	amat::Array2d<a_int> elemface;

	amat::Array2d<a_int> bifmap;				///< relates boundary faces in intfac with bface, ie, bifmap(intfac no.) = bface no.
	amat::Array2d<a_int> ifbmap;				///< relates boundary faces in bface with intfac, ie, ifbmap(bface no.) = intfac no.
	bool isBoundaryMaps;						///< Specifies whether bface-intfac maps have been created

	/// Compute lists of elements surrounding points \ref esup
	void compute_elementsSurroundingPoints();

	/// Compute lists of elements surrounding elements
	void compute_elementsSurroundingElements();

	/// Computes \ref intfac and \ref elemface
	/** - Computes, for each face,
	 *    > the elements on either side,
	 *    > the starting node and the ending node of the face.
	 * This is stored in \ref intfac. Elements sourrounding elements (\ref esuel) and elements surrounding
	 * points (\ref esup) are assumed to be available.
	 * Esuel is modified to hold (\ref nelem + face no.) for each ghost cell, instead of -1 as before.
	 */
	void compute_face_structure();

	std::vector<std::pair<a_int,int>> compute_phyBFaceNeighboringElements() const;

	/// Currently unused, but supposed to compute lists of points surrounding each point
	void compute_pointsSurroundingPoints();
};

UMesh2dh prepare_mesh(const std::string meshfile);


} // end namespace
#endif
