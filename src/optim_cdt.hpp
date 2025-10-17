#include <iostream>
#include <vector>
#include <numeric>
#include "cdt.h"
#include "tetmesh.h"

class bnd_edge {
	public:
	uint32_t e0, e1;
	uint32_t count;
	bool is_link;
	bool visited;

	bnd_edge(uint32_t v0, uint32_t v1) : e0(v0), e1(v1), count(1), is_link(false), visited(false) { if(e0 > e1) std::swap(e0,e1); };

	bool same_ep(const bnd_edge& b) const { return (e0 == b.e0 && e1 == b.e1);}
	bool has_ep(uint32_t ep) const { return (ep == e0 || ep == e1); }
	uint32_t opposite_ep(uint32_t v) const { assert(v==e0 || v==e1); if(v==e0) return e1; else return e0; }
	uint32_t common_ep(const bnd_edge& b) const { if(b.has_ep(e0)) return e0; if(b.has_ep(e1)) return e1; else return UINT32_MAX; }
	uint32_t non_common_ep(const bnd_edge& b) const { if(b.has_ep(e0)) return e1; if(b.has_ep(e1)) return e0; else return UINT32_MAX; }
	bool isUnique() const { return count == 1; }

    // Static functions to be used as predicates in std algorithms
    static inline bool isNotUniquePtr(const bnd_edge& b) { return !b.isUnique(); }

	bool operator<(const bnd_edge& b) const { return (e0 < b.e0 || (e0 == b.e0 && e1 < b.e1) ); }
};

// Given a triangulated surface with boundaries, whose triangle are defined as triple
// of indices throgh surf_tri, fills be with edges composing the boundary.
void get_chamPLC_bnd_edges(std::vector<uint32_t>& surf_tri, std::vector<bnd_edge>& be){
	// Make half-edges
    for(uint32_t i=0; i<surf_tri.size()/3; i++){
		const uint32_t* v = surf_tri.data() + i*3;
		be.insert(be.end(), {bnd_edge(v[0],v[1]), bnd_edge(v[1],v[2]), bnd_edge(v[2],v[0])} );
    }
	// Keep half-edges that occour only one time, i.e. those that belong to the boundary.
    // Forall i, be[i] endpoints e0 and e1 are s.t. e0 < e1 by construction, 
    // now we sort them by lexicographic endpoint order.
	sort( be.begin(), be.end() ); 
    for(size_t i=0; i<be.size()-1; i++) if( be[i].same_ep(be[i+1]) ){
        be[i+1].count += be[i].count; // correct number of occurences
        be[i].count = 0; // its occorrence is > 1, so it will be deleted later
    }
    // Remove all edges that have not exactly 1 occurence.
    be.erase( std::remove_if(be.begin(), be.end(), bnd_edge::isNotUniquePtr), be.end());
}

void get_vrt_chain_from_edge_chain(const std::vector<bnd_edge*>& eChain, std::vector<uint32_t>& vChain){
	size_t n_vrts = eChain.size()+1; // the number of vertices is the number of edges +1
	vChain.resize(n_vrts, UINT32_MAX);
	if(n_vrts == 2) { vChain[0] = eChain[0]->e0; vChain[1] = eChain[0]->e1; return; }
	vChain[0] = eChain[0]->non_common_ep( *(eChain[1]) );
	for(size_t i=1; i<n_vrts; i++) vChain[i] = eChain[i-1]->opposite_ep( vChain[i-1] );
}

// This DS supports reconstruction around chamfered edges:
// - r0 and r1 are the endpoints of a chamfered input edge, when it has been chamfered out
//   it generated a number of new edges E_i equal to the number of input triangles
//	 incident at <r0,r1> ;
// - v are the inices of the vertices that the Delauny refinement algorithm has placed
// 	 on a certain segment E_i;
// - comm_v are the indices of the vertices that are placed on <r0,r1> to triangulate 
//	 the "strip" delimited by <r0,r1> and E_i (note that comm_v must be the same for
//	 each of the E_i)
class bnd_vrt_chain {
	public:
	std::vector<uint32_t> v;
	std::vector<uint32_t> comm_v;
	uint32_t r0, r1;

	bnd_vrt_chain(){};
	
	void add_comm_v(size_t first_i, size_t last_i){ 
		comm_v.resize( last_i - first_i );
		std::iota(comm_v.begin(), comm_v.end(), (uint32_t)first_i);
	}
	void set_ref_vrts(uint32_t _r0, uint32_t _r1){ r0=_r0; r1=_r1; }

	uint32_t e0() const { return v.front(); }
	uint32_t e1() const { return v.back(); }
	size_t size_v() const { return v.size(); }
	size_t size_comm_v() const { return comm_v.size(); }
	void reverse_v() { std::reverse(v.begin(),v.end()); std::swap(r0,r1); }
	void reverse_comm_v() { std::reverse(comm_v.begin(),comm_v.end()); }

	bool same_ref_vrts(const bnd_vrt_chain& c) const { 
		return (r0 == c.r0 && r1 == c.r1) || (r1 == c.r0 && r0 == c.r1);}
	bool operator<(const bnd_vrt_chain& c) const { return (r0 < c.r0 || ( r0 == c.r0 && r1 < c.r1 )); }
};

uint32_t howMany_newVrts_onInputEdge(const bnd_vrt_chain& c1, const bnd_vrt_chain& c2, const std::vector<pointType*>& vertices) {
	const pointType* p1 = vertices[c1.e0()];
	const pointType* p2 = vertices[c2.e0()];
	const pointType* U = vertices[c1.r0];
	const pointType* V = vertices[c1.r1];
	const double dist_p1_UV = sqrt( pointSqDistanceFromLine(p1, U, V) );
	const double dist_p2_UV = sqrt( pointSqDistanceFromLine(p2, U, V) );
	const double won1 = dist_p2_UV / (dist_p2_UV + dist_p1_UV);
	return (uint32_t)floor( (c1.size_v() * won1 + c2.size_v() * (1.0-won1) ) );
}


// Mark all explicit vertices within distance 'mind' from 'p'
double markChamferExplicitNeighbors(Tetrahedrization& mesh, const pointType *p, double mind) {
	// Search the tet containing 'p'
	Tetrahedron *t0 = mesh.searchTet(p);
	if (t0 == NULL) ip_error("markChamferExplicitNeighbors: could not find 'p' in  'mesh'\n");

	// Set 'mind' to the distance of the closest implicit vertex to p
	const pointType* tp;
	double d0, d1, d2, d3;
	Tetrahedra cavity;
	Tetrahedron* s, * t;
	cavity.push_back(t0); t0->mark<7>();
	for (size_t i = 0; i < cavity.size(); i++) {
		t = cavity[i];
		tp = t->v0()->getPoint(); d0 = sqEuclideanDistance(tp, p); if (tp->isExplicit3D() && d0 < mind) { t->v0()->mark<0>(); }
		tp = t->v1()->getPoint(); d1 = sqEuclideanDistance(tp, p); if (tp->isExplicit3D() && d1 < mind) { t->v1()->mark<0>(); } 
		tp = t->v2()->getPoint(); d2 = sqEuclideanDistance(tp, p); if (tp->isExplicit3D() && d2 < mind) { t->v2()->mark<0>(); } 
		tp = t->v3()->getPoint(); d3 = sqEuclideanDistance(tp, p); if (tp->isExplicit3D() && d3 < mind) { t->v3()->mark<0>(); }

		s = t->t0(); if (s != NULL && !s->isMarked<7>() && (d1 <= mind || d2 <= mind || d3 <= mind)) { cavity.push_back(s); s->mark<7>(); }
		s = t->t1(); if (s != NULL && !s->isMarked<7>() && (d0 <= mind || d2 <= mind || d3 <= mind)) { cavity.push_back(s); s->mark<7>(); }
		s = t->t2(); if (s != NULL && !s->isMarked<7>() && (d1 <= mind || d0 <= mind || d3 <= mind)) { cavity.push_back(s); s->mark<7>(); }
		s = t->t3(); if (s != NULL && !s->isMarked<7>() && (d1 <= mind || d2 <= mind || d0 <= mind)) { cavity.push_back(s); s->mark<7>(); }
	}
	for (Tetrahedron* t : cavity) { t->unmark<7>(); }

	return mind;
}

// Marks all explicit vertices within distance 'mind' from segment p1-p2.
void markChamferExplicitNeighbors(Tetrahedrization& mesh, const pointType* p1, const pointType* p2, double mind) {
	// Search the tet containing 'p1'
	Tetrahedron* t0 = mesh.searchTet(p1);
	if (t0 == NULL) ip_error("markChamferExplicitNeighbors: could not find 'p1' in 'mesh'\n");

	// Set 'mind'
	const pointType* tp;
	double d0, d1, d2, d3;
	Tetrahedra cavity;
	Tetrahedron* s, * t;
	cavity.push_back(t0); t0->mark<7>();
	for (size_t i = 0; i < cavity.size(); i++) {
		t = cavity[i];
		tp = t->v0()->getPoint(); d0 = vector3d(tp).sq_dist_segment(p1, p2); if (tp->isExplicit3D() && d0 < mind) { t->v0()->mark<0>(); }
		tp = t->v1()->getPoint(); d1 = vector3d(tp).sq_dist_segment(p1, p2); if (tp->isExplicit3D() && d1 < mind) { t->v1()->mark<0>(); }
		tp = t->v2()->getPoint(); d2 = vector3d(tp).sq_dist_segment(p1, p2); if (tp->isExplicit3D() && d2 < mind) { t->v2()->mark<0>(); }
		tp = t->v3()->getPoint(); d3 = vector3d(tp).sq_dist_segment(p1, p2); if (tp->isExplicit3D() && d3 < mind) { t->v3()->mark<0>(); }

		s = t->t0(); if (s != NULL && !s->isMarked<7>() && (d1 <= mind || d2 <= mind || d3 <= mind)) { cavity.push_back(s); s->mark<7>(); }
		s = t->t1(); if (s != NULL && !s->isMarked<7>() && (d0 <= mind || d2 <= mind || d3 <= mind)) { cavity.push_back(s); s->mark<7>(); }
		s = t->t2(); if (s != NULL && !s->isMarked<7>() && (d1 <= mind || d0 <= mind || d3 <= mind)) { cavity.push_back(s); s->mark<7>(); }
		s = t->t3(); if (s != NULL && !s->isMarked<7>() && (d1 <= mind || d2 <= mind || d0 <= mind)) { cavity.push_back(s); s->mark<7>(); }
	}
	for (Tetrahedron* t : cavity) { t->unmark<7>(); }
}

// just for DEBUG purposes 
bool check_overlaps(const std::vector<pointType *> vrts, const std::vector<uint32_t>& tri){
	bool passed = true;
	std::vector<bnd_edge> e;
	const uint32_t* t;
	for(size_t i=0; i<tri.size()/3; i++){
		t = tri.data() + 3*i;
		e.push_back( bnd_edge(*t     , *(t+1)) ); e.back().count = (uint32_t)i;
		e.push_back( bnd_edge(*(t+1) , *(t+2)) ); e.back().count = (uint32_t)i;
		e.push_back( bnd_edge(*(t+2) , *t	 ) ); e.back().count = (uint32_t)i;
	}

	uint32_t u0, u1, v, w;
	std::vector<uint32_t> otri;
	for(size_t i=0; i<e.size()-1; i++) for(size_t j=i+1; j<e.size(); j++){
		if( e[i].same_ep(e[j]) ){
			u0 = e[i].e0; u1 = e[i].e1;
			t = tri.data() + 3*e[i].count;
			v = *t; if(v==u0 || v==u1){ v = *(t+1); if(v==u0 || v==u1) v = *(t+2); }
			t = tri.data() + 3*e[j].count;
			w = *t; if(w==u0 || w==u1){ w = *(t+1); if(w==u0 || w==u1) w = *(t+2); }
			assert(v!=u0 && v!=u1 && w!=u0 && w!=u1 && v!=w);
			int o3d = pointType::orient3D(*vrts[u0],*vrts[u1],*vrts[v],*vrts[w]);
			if(o3d == 0 && isAcuteDihedral_exact( vrts[u0], vrts[u1], vrts[v], vrts[w])){
				std::cout<<"ERROR overlapping founded\n";
				t = tri.data() + 3*e[i].count;
				otri.insert(otri.end(),{*t,*(t+1),*(t+2)});
				std::cout<<"tri["<<i<<"] = <"<<*t<<" "<<*(t+1)<<" "<<*(t+2)<<">\n";
				t = tri.data() + 3*e[j].count;
				otri.insert(otri.end(),{*t,*(t+1),*(t+2)});
				std::cout<<"tri["<<i<<"] = <"<<*t<<" "<<*(t+1)<<" "<<*(t+2)<<">\n";
				passed = false;
			}
		}
	}

	if(passed) std::cout<<"overlap check passed\n";
	else{
		FILE* fp = fopen("overlaps.off", "w");
		fprintf(fp, "OFF\n%u %u 0\n", (uint32_t)vrts.size(), (uint32_t) otri.size() / 3);

		for(uint32_t i=0; i<vrts.size(); i++) {
			double x,y,z;
			vrts[i]->getApproxXYZCoordinates(x,y,z);
			fprintf(fp, "%f %f %f\n", x, y, z);
		}

		for(size_t i=0; i<otri.size()/3; i++) {
			std::vector<uint32_t> v;
			v.assign(otri.begin() + 3*i, otri.begin() +3*i +3);
			fprintf(fp, "3 %u %u %u \n", v[0], v[1], v[2]);
		}
			
		fclose(fp);
	}

	return passed;
}

bool get_vrts_and_tris_for_cdt(Tetrahedrization& mesh, std::vector<uint32_t>& ref_vrts, 
							   std::vector<genericPoint*>& cdt_vrts, std::vector<uint32_t>& cdt_tris,
							   const char* filename, bool verbose){

	// DEBUG start
	// mesh.checkAllFaces(); // DEBUG
    // mesh.checkConnectivity(); // DEBUG
	// std::cout<<"\n[get_vrts_and_tris_for_cdt] INPUT:\n";
	// std::cout<<"\nref_vrts size: "<<ref_vrts.size()<<"\n";
	// for(size_t i=0; i<ref_vrts.size(); i++) std::cout<<"ref_vrt["<<i<<"] = "<<ref_vrts[i]<<"\n";
	// std::cout<<"\n";
	// std::cout<<"(cdt_vrts.size() = "<<cdt_vrts.size()<<", "<<"cdt_tris.size() = "<<cdt_tris.size()<<"\n";
	// DEBUG end
	
	// Exctract triangles of the optimized mesh that are part of the chamfered surface (constr_tri)
	std::vector<uint32_t> constr_tri;
	mesh.export_DelTris_asTriVrtsInds(constr_tri, true); // true -> do not include bounding box triangles
	size_t n_constr_tri = constr_tri.size()/3;

	// DEBUG start
	// std::cout<<"\nDR mesh DEL TRIS: "<<n_constr_tri<<" (size = "<<constr_tri.size()<<") \n";
	// for(size_t i=0; i<n_constr_tri; i++){ 
	// 	std::cout<<"constr_tri["<<i<<"] = < "
	// 			<<constr_tri[3*i]<<" "
	// 			<<constr_tri[3*i+1]<<" "
	// 			<<constr_tri[3*i+2]<<" >\n";
	// }
	// std::cout<<"\n";
	// mesh.checkAllFaces(); // DEBUG
    // mesh.checkConnectivity(); // DEBUG
	// mesh.savePLCFaces("DR_PLCfaces.off");
	// DEBUG end

	// constr_tri triangles exactly conform the chamfered version of the input PLC.
	// The chamfered surface is a "reduction" of the input PLC.
	// The objetive of this function is to add new trinagles to constr_tri in order 
	// to exactcly conform the input PLC. 
	// To this end we have to exctract boundary edges of the surface defined by 
	// constr_tri which is a refinement of the boundary edges of the chamfered surface.
	// Each one of these boundary edges will be connected to an appropriate vertex,
	// in order to create triangles conforming the complementar region of constr_tri
	// wrt the input PLC (which will be refered as co-chamfered region).
	// ref_vrts is a vector of lenght equal to the number of input PLC vertices +
	// the number of vertices introduced by the chamfering: 
	// it associate to each input PLC vertex UINT32_MAX, while to each other
	// vertex V associate the index of the input vertex U that while removed by the
	// chamfering produced V.

	// Collect the boundary edges of constr_tri.
	std::vector<bnd_edge> be; get_chamPLC_bnd_edges(constr_tri, be);

	// DEBUG start
	// std::cout<<"\nDR mesh BND EDGES: "<<be.size()<<"\n";
	// for(size_t i=0; i<be.size(); i++){ 
	// 	std::cout<<"bnd_edge["<<i<<"] = < "
	// 			<<be[i].e0<<" "
	// 			<<be[i].e1<<" >\n";
	// 	//std::cout<<*(mesh.vrts()[be[i].e0]->getPoint())<<"\n";
	// 	//std::cout<<*(mesh.vrts()[be[i].e1]->getPoint())<<"\n\n";
	// }
	// std::cout<<"\n";
	// mesh.checkAllFaces(); // DEBUG
    // mesh.checkConnectivity(); // DEBUG
	// DEBUG end

	// Make the vertex-edge relation wrt boundary edges
	size_t n_optimMesh_vrts = mesh.num_vertices();
	std::vector< std::vector<uint32_t> > vbe_rel( n_optimMesh_vrts );
	for(size_t be_i=0; be_i<be.size(); be_i++){
		vbe_rel[ be[be_i].e0 ].push_back((uint32_t)be_i);
		vbe_rel[ be[be_i].e1 ].push_back((uint32_t)be_i);
	}

	// To connect a boundary edge be_i (of be) to a suitable vertex in order to
	// create a valid triangle, we first collect all its connected boundary edges
	// in to a 'chain' convering the whole edge <V0,V1> of the chamfered surface 
	// and secondly distinguish two cases:
	// - (1) the 'chain' <V0,V1> was introduced while chamfering the input PLC vertex U,
	//	 in such case the endpoints of all the be_i forming <V0,V1> have to
	//	 be connected to U; ref_vrts will be enlarged to store this information (U)
	//	 for each of the vertices of the chain and trinagles will be created later on.
	// - (2) the 'chain' <V0,V1> was introduced while chamfering the input edge <U0,U1>,
	//	 in such case we follow a more complicated procedure as to form 'decent' trinagles
	//	 a certain number of new vertices on the input edge <U0,U1> have to be introduced.
	//	 This refinement of <U0,U1> have to be the same for each chain generated by the
	//	 chamfering of <U0,U1> (two for manifold cases, or more for non-manifold ones).
	//	 We use the 'strip' data strucuture to associate each chain to the input
	//	 edge <U0,U1> and later uniquely refine <U0,U1> and finally trianguale each strip.
	size_t n_cham_vrts = ref_vrts.size();
	ref_vrts.resize(n_optimMesh_vrts, UINT32_MAX);

	std::vector< bnd_edge* > chain;
	std::vector<bnd_vrt_chain> half_strip_bnd;
	uint32_t v0,v1,v2, r0,r1;
	uint32_t curr_e;
	for(size_t vi=0; vi<n_cham_vrts; vi++) if(ref_vrts[vi]!=UINT32_MAX){
		// Consider each boundary edge incident at vi (at least 2)
		for(uint32_t bei : vbe_rel[vi]) if(!be[bei].visited) {

			// Step 1 : form the chain.
			// Each chamfered-PLC-boundary-edge E = <b0,b1>
			// has ref_vrt[b0]!=UINT32_MAX and ref_vrt[b1]!=UINT32_MAX.
			// During Delaunay refinement new vertices may be inserted on E,
			// each one of them has exactly two incident edges. 
			// Starting from b0 = v0 we want to reach b1 collecting all
			// sub-edges of E in a unique chain.

			v0 = (uint32_t)vi; assert( ref_vrts[v0] != UINT32_MAX ); // starting vertex
			v1 = be[bei].opposite_ep((uint32_t)vi);  assert(vi!=v1);
			curr_e = bei;
			chain.push_back( &(be[curr_e]) ); be[curr_e].visited = true;
			while( ref_vrts[v1] == UINT32_MAX ){
				v0 = v1;
				assert(vbe_rel[v0].size() == 2);
				if( curr_e == vbe_rel[v0][0] ) curr_e = vbe_rel[v0][1];
				else 						   curr_e = vbe_rel[v0][0];
				v1 = be[curr_e].opposite_ep(v0);
				chain.push_back( &(be[curr_e]) ); be[curr_e].visited = true;
			}

			// Step 2 : distinguish between cases (1) and (2)
			// vi is the "firts" vertex of the chain, while v1 is the "last".
			if(ref_vrts[vi] == ref_vrts[v1]) {
				assert(vi!=v1);
				// Having the same reference vertex U it is the case (1)
				for(bnd_edge* b : chain) ref_vrts[b->e0] = ref_vrts[b->e1] = ref_vrts[vi]; 
			}
			else{
				// Having different reference vertices U0,U1 it is the case (2)
				// we have to initialize the 'strip' data structure
				for(bnd_edge* b : chain) b->is_link = true;
				half_strip_bnd.push_back( bnd_vrt_chain() );
				get_vrt_chain_from_edge_chain( chain, half_strip_bnd.back().v );
				r0 = ref_vrts[ half_strip_bnd.back().e0() ];
				r1 = ref_vrts[ half_strip_bnd.back().e1() ];
				assert(r0!=r1);
				if(r0 > r1){ 
					half_strip_bnd.back().reverse_v(); // to lexicographic sort half_strip_bnd later on.
					std::swap(r0, r1);
				}
				half_strip_bnd.back().set_ref_vrts(r0, r1);
			}
			chain.clear();
		}
	}

	// Identify opposite half-strips, i.e. those that have same couple of ref_vrts.
	// To later triangulate we need to insert new vertices on the shared input edge.
	
	// Creting a new vertices vector to add new vertices.
	std::vector<pointType*> vertices(n_optimMesh_vrts);
	for(size_t i=0; i<n_optimMesh_vrts; i++) vertices[i] = mesh.vrts()[i]->getPoint();

	std::sort(half_strip_bnd.begin(), half_strip_bnd.end());

	// assert(half_strip_bnd.size()%2 == 0);

	for(size_t i=0; i<half_strip_bnd.size(); i+=2){
		bnd_vrt_chain& c1 = half_strip_bnd[i];
		bnd_vrt_chain& c2 = half_strip_bnd[i+1];

		assert(c1.same_ref_vrts(c2));
		
		r0 = c1.r0; r1 = c1.r1;
		if( r0 == c2.r1 ) c2.reverse_v();
		
		uint32_t n_pts = howMany_newVrts_onInputEdge(c1,c2, vertices); assert(n_pts > 1);
		const size_t init_vert_size = vertices.size();
		const pointType* P = vertices[r0];
		const pointType* Q = vertices[r1];
		assert(vertices[r0]->isExplicit3D() && vertices[r1]->isExplicit3D());
		const vector3d OP(P), OQ(Q);
		double l = sqrt( OP.dist_sq(OQ) );
		double d0 = sqrt( OP.dist_sq( vector3d( vertices[c1.e0()] ) ) );
		double d1 = sqrt( OQ.dist_sq( vector3d( vertices[c1.e1()] ) ) );
		double t = d0 / l; assert(0<t && t<1);
		vertices.push_back( new implicitPoint3D_LNC(P->toExplicit3D(), Q->toExplicit3D(), t) );
		if(n_pts > 2){
			double h = (l - (d0+d1)) / (n_pts-1);
			double dt = h/l;
			for(size_t j=1; j<n_pts-1; j++){
				t += dt; assert(0<t && t<1);
				vertices.push_back( new implicitPoint3D_LNC(P->toExplicit3D(), Q->toExplicit3D(), t) );
			}
		}
		t = 1.0 - d1 / l; assert(0<t && t<1);
		vertices.push_back( new implicitPoint3D_LNC(P->toExplicit3D(), Q->toExplicit3D(), t) );
		half_strip_bnd[i].add_comm_v(init_vert_size, vertices.size());
		half_strip_bnd[i+1].add_comm_v(init_vert_size, vertices.size());
	}

	// build the comlementary triangles of the chamfered optimized surface wrt input
    std::vector<uint32_t> compl_tri;
    for(const bnd_edge& b : be) if(!b.is_link) {
		if(vertices[b.e0]->isExplicit3D() || vertices[b.e1]->isExplicit3D() ) continue;
		assert(ref_vrts[b.e0] == ref_vrts[b.e1]);
		compl_tri.insert( compl_tri.end(), {b.e0, b.e1, ref_vrts[b.e0]} );
		// orientation maybe have to be changed
    }

	uint32_t n_new_vrts = (uint32_t)vertices.size();
	// uint32_t n_tot_vrts = (uint32_t)n_optimMesh_vrts + n_new_vrts;
	for(const bnd_vrt_chain& c : half_strip_bnd){
		v0 = c.e0();
		v1 = c.comm_v.front();
		v2 = ref_vrts[v0];
		compl_tri.insert( compl_tri.end(), {v0, v1, v2} );

		size_t i = 0; // indexing vertices of c.v
		size_t j = 0; // indexing vertices of c.comm_v
		size_t i_max = c.v.size()-1;
		size_t j_max = c.comm_v.size()-1;

		v0 = c.e0();
		v1 = c.comm_v[0];
		double dv0u1, dv1u0;
		uint32_t u0, u1;
		while (i < i_max || j < j_max){

			if(i<i_max){
				u0 = c.v[i+1];
				dv1u0 = vector3d( vertices[u0] ).dist_sq(vector3d( vertices[v1] ));
			}
			else{
				u0 = UINT32_MAX;
				dv1u0 = DBL_MAX;
			}
			
			if(j<j_max){
				u1 = c.comm_v[j+1];
				dv0u1 = vector3d( vertices[v0] ).dist_sq(vector3d( vertices[u1] ));
			}
			else{
				u1 = UINT32_MAX;
				dv0u1 = DBL_MAX;
			}
			
			
			if(dv0u1 < dv1u0){
				assert(u1!=UINT32_MAX);
				compl_tri.insert( compl_tri.end(), {v0, v1, u1} );
				v1 = u1;
				j++;
			}
			else{
				assert(u0!=UINT32_MAX);
				compl_tri.insert( compl_tri.end(), {v0, u0, v1} );
				v0 = u0;
				i++;
			}
		}

		v0 = c.e1();
		v1 = c.comm_v.back();
		v2 = ref_vrts[v0];
		compl_tri.insert( compl_tri.end(), {v0, v1, v2} );
	}

	constr_tri.insert(constr_tri.end(), compl_tri.begin(), compl_tri.end());
	assert( (constr_tri.size() % 3) == 0 );

	assert( check_overlaps(vertices, constr_tri) );

	// Remove Delaunay refined mesh vertices too close to CDT vertices
	uint32_t n_input_exp_vrts = 0;
	double d;
	while(vertices[n_input_exp_vrts]->isExplicit3D()) n_input_exp_vrts++;
	std::vector<double> ch_dist(n_input_exp_vrts, 0.0);
	for(size_t i=n_input_exp_vrts; i<n_cham_vrts; i++) {
		assert( ref_vrts[i] != UINT32_MAX );
		d = vector3d( vertices[ ref_vrts[i] ] ).dist_sq(vector3d( vertices[i] ));
		ch_dist[ ref_vrts[i] ] = max(ch_dist[ ref_vrts[i] ], d);
	}
	for(TetVertex* v : mesh.vrts()) v->unmark<0>(); // reset markers
	for(size_t i=0; i<n_input_exp_vrts; i++) markChamferExplicitNeighbors(mesh, vertices[i], ch_dist[i]);
	for(const bnd_vrt_chain& strip : half_strip_bnd){ 
		assert(strip.r0 < n_input_exp_vrts && strip.r1 < n_input_exp_vrts );
		d = max(ch_dist[strip.r0], ch_dist[strip.r1]);
		markChamferExplicitNeighbors(mesh, vertices[ strip.r0 ], vertices[ strip.r1 ], d); 
	}
	std::vector<uint32_t> uv(vertices.size(), 1); // Used Vertex 
	assert(vertices.size() == mesh.num_vertices());
	size_t count_removed = 0;
	//for(size_t i = n_cham_vrts; i < n_optimMesh_vrts; i++) if(mesh.vrts()[i]->getPoint()->isExplicit3D()){
	for(size_t i = n_cham_vrts; i < n_optimMesh_vrts; i++) if(mesh.vrts()[i]->isMarked<0>()){ 
		// do not use marked vertices inserted by the optimizer
		uv[i] = UINT32_MAX;
		if(verbose) count_removed++;
	} 
	if(verbose) std::cout<<count_removed<<" vertices removed.\n";
	if(verbose){
		size_t count_marked = 0;
		for(size_t i = 0; i < n_cham_vrts; i++) if(mesh.vrts()[i]->isMarked<0>()){ 
			count_marked++;
		} 
		if(count_marked) std::cout<<count_marked<<" vertices marked but not removed.\n";
	}
	for(TetVertex* v : mesh.vrts()) v->unmark<0>(); // reset markers

	uint32_t idx = 0;
	for (size_t i = 0; i < uv.size(); i++) if (uv[i] != UINT32_MAX) uv[i] = idx++; // now uv stores new indexing
	for (size_t i=0; i<constr_tri.size(); i++){ 
		assert(constr_tri[i]<uv.size());
		constr_tri[i] = uv[ constr_tri[i] ];
	}

	for(uint32_t i=0; i<vertices.size(); i++) if(uv[i]!=UINT32_MAX) {
		cdt_vrts.push_back( vertices[i] );
	}

	cdt_tris.assign(constr_tri.begin(), constr_tri.end());

	
	if(strlen(filename) != 0){
		// Saves the surface on an .off file, it is valid input for a cdt algorithm
		char out_filename[2048]; // following lines remove folder path and extension to file_name
		strcpy(out_filename, filename);
		char* tok = strtok(out_filename, "/");
		while(tok != NULL){ strcpy(out_filename, tok); tok = strtok(NULL, "/");  }
		strcpy(out_filename, strtok(out_filename, "."));
		strcat(out_filename, "_rebuilt.off");
		// ----------------------------------------
		FILE* fp = fopen(out_filename, "w");
		fprintf(fp, "OFF\n%u %u 0\n", (uint32_t)cdt_vrts.size(), (uint32_t) cdt_tris.size() / 3);

		for(uint32_t i=0; i<cdt_vrts.size(); i++) {
			double x,y,z;
			cdt_vrts[i]->getApproxXYZCoordinates(x,y,z);
			fprintf(fp, "%f %f %f\n", x, y, z);
		}

		for(size_t i=0; i<cdt_tris.size()/3; i++) {
			std::vector<uint32_t> v;
			v.assign(cdt_tris.begin() + 3*i, cdt_tris.begin() +3*i +3);
			fprintf(fp, "3 %u %u %u \n", v[0], v[1], v[2]);
		}
			
		fclose(fp);
	}

	// check_overlaps(cdt_vrts, cdt_tris); // DEBUG

	if(verbose) std::cout<<"valid CDT input generated\n";

	return true;
}