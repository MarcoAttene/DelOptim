#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <set>


// ------ //
// POINTS //
// ------ //

// Returns the baricentric coordinates of BPT 'B'.
// NOTE. the implicitPointBPT constructor takes as argument 3 points and 2 
//       doubles: implicitPointBPT(P, Q, R, v, u),
//       to build the point 'B' = P * v + Q * u + R * (1 - v - u)
//       since the argument order is confusing (P came before Q, while u came 
//       after v) I preferet to itroduce the change of variables 
//       'xi' := v, 'eta' := u, which brings to
//       implicitPointBPT(P, Q, R, 'xi', 'eta') and
//       B = P * 'xi' + Q * 'eta' + R * (1 - 'xi' - 'eta')
inline void get_baricentric_coords(double& xi, double& eta, const pointType* B){
    assert(B->isBPT());
    xi = B->toBPT().V();  eta = B->toBPT().U();
}

// ------ //
// ANGLES //
// ------ //

// Let be "tau" the plane for the triangle t = <'q','r','s'>.
// Returns TRUE if both following conditions holds:
// 1) segment <'q','p'> forms an acute angle at 'q' with the triangle t,
// 2) the projection of <'q','p'> on "tau" is inside the sector limited by 
//    half-straight-lines <'q','r'> and <'q','s'> containing the triangle t.
// It is equivalent to check acuteness of both dihedral angles
// - between <'p','q','s'> and <'r','q','s'> and
// - between <'p','q','r'> and <'s','q','r'>.
bool isAcuteAngle(const pointType* p, const pointType* q, 
                    const pointType* r, const pointType* s) {
    return (isAcuteDihedral_exact(s,q,r,p) && isAcuteDihedral_exact(r,q,s,p));
}

// ------------------ //
// IMPLICIT RELATIONS //
// ------------------ //

// Supports only simply-connected faces (each vertex compares at most 1 time on 
// face boundary). Fills 'fv' with the indices of face 'f' vertices.
void PLCc::get_face_vertices(const CHAMface& f, 
                                std::vector<uint32_t>& fv) const {

    const std::vector<uint32_t>& fbe = f.bounding_edges;
    fv.resize( fbe.size() );
    for(size_t i = 0; i < fbe.size()-1; i++) if(fbe[i] != EMPTY_PLACE) {
        const CHAMedge& e1 = edges[ fbe[i] ];
        const CHAMedge& e2 = edges[ fbe[i+1] ];

        assert( e1.has_commonVertex( e2 ) );        
        
        fv[i] = e1.commonVertex( e2 );
    }

    assert( edges[ fbe.back() ].has_commonVertex( edges[ fbe.front() ] ) );
        
    fv.back() = edges[ fbe.back() ].commonVertex( edges[ fbe.front() ] );
}

// ------------------------- //
// OREDERING OF PLCC VECTORS //
// ------------------------- //

void PLCc::swap_edges(uint32_t e1, uint32_t e2){

    // update boundary edges of faces incident at edges[e1]
    for(size_t i=0; i<edges[e1].inc_face.size(); i++)
        faces[ edges[e1].inc_face[i] ].replaceEdge_11(e1, e2);

    // update boundary edges of faces incident at edges[e2]
    for(size_t i=0; i<edges[e2].inc_face.size(); i++)
        faces[edges[e2].inc_face[i]].replaceEdge_11(e2, e1);

    std::swap(edges[e1], edges[e2]); // swap edges elements
    std::swap(mark_edges[e1], mark_edges[e2]);
}

uint32_t PLCc::move_back_isolated_edges(){
    uint32_t last = ( (uint32_t) edges.size() ) -1;
    for(uint32_t i = 0; i < edges.size(); i++) if(edges[i].isIsolated()){ 
        while( edges[last].isIsolated() ) last--;
        if(last > i){ swap_edges(i, last); }
    }
    while( edges[last].isIsolated() ) last--;
    return ++last;
}

// ------------------- //
// PLCc INITIALIZATION //
// ------------------- //

// Returns TRUE if triangle 'u' vertices u[0], u[1], u[2] are an odd permutation 
// triangle 'v' vertices v[0], v[1], v[2].
inline bool are_flipped_triangles(const uint32_t* u, const uint32_t* v){
    return( (v[0]==u[1] && v[1]==u[0] && v[2]==u[2]) ||
            (v[0]==u[2] && v[1]==u[1] && v[2]==u[0]) ||
            (v[0]==u[0] && v[1]==u[2] && v[2]==u[1])    );
}

void PLCc::orient_initial_triface_bnd(CHAMface& f){
    const uint32_t* t = plc.triangle_vertices.data() + (f.triangle*3);
    const CHAMedge& e0 = edges[ f.bounding_edges[0] ];
    const CHAMedge& e1 = edges[ f.bounding_edges[1] ];
    const CHAMedge& e2 = edges[ f.bounding_edges[2] ];
    uint32_t u[] = { e0.commonVertex(e1), 
                     e1.commonVertex(e2), 
                     e2.commonVertex(e0) };
    if( are_flipped_triangles(u, t) ) 
        std::swap(f.bounding_edges[0], f.bounding_edges[2]);
}

// Removes duplicated pre-edges (transorming half-edges into edges) 
void PLCc::mergePreEdges(){
    // Sort half-edges by lexicografic non-descending order wrt half-edge 
    // endpoints. Each half-edge has ep[0] < ep[1].
    std::sort(edges.begin(), edges.end(), CHAMedge::vertexSortFunc);

    for(uint32_t ei=0; ei < edges.size()-1; ){
        CHAMedge& e = edges[ei];

        while( ++ei < edges.size() && edges[ei].coincident(e) ) {
            CHAMedge& ne = edges[ei];
            // each edge initially has only one incident triangle
            e.inc_face.push_back( ne.inc_face[0] );
            ne.inc_face.clear(); // 'ne' is now isolated
        }
    }

    // Remove duplicated edges (no inc_faces) from edges vectcor
    edges.erase( 
        std::remove_if(edges.begin(), edges.end(), CHAMedge::isIsolatedPtr), 
        edges.end() );
}

// Fill the PLCc data structure by using the input triangulation information
void PLCc::initialize(){
    
    // -- Vertices -- 
    vertices.reserve(n_in_vrts);
    for(uint32_t i=0; i<n_in_vrts; i++) {
        const double* x = plc.coordinates.data() + i*3;
        vertices.push_back( new explicitPoint(*x, *(x +1), *(x +2)) ); 
    }
    mark_vrts.resize(vertices.size(),0);
    ref_exp3D_vrt.resize(vertices.size(), UINT32_MAX);

    // -- Edges --
    edges.resize(plc.numTriangles() * 3); // More than necessary, will be merged

    // Fill only EV and partial ET relations
    // also fill vt (input_vertex - incident_input_triangle relation)
    input_vt.resize(vertices.size());
    uint32_t pei, v0,v1,v2;
    const uint32_t* t_ref = plc.triangle_vertices.data();
    for(uint32_t ti=0; ti<plc.numTriangles(); ti++){
        pei = ti*3;
        v0 = *(t_ref + pei   );  
        v1 = *(t_ref + pei +1);  
        v2 = *(t_ref + pei +2);
    
        // pre-edges
        edges[pei  ].fill_preEdge(v0, v1, ti);
        edges[pei+1].fill_preEdge(v1, v2, ti);
        edges[pei+2].fill_preEdge(v2, v0, ti);

        // vt
        input_vt[v0].push_back(ti);
        input_vt[v1].push_back(ti);
        input_vt[v2].push_back(ti);
    }

    // merge duplicated pre-edges and fill edges
    mergePreEdges(); // Now edges contains proper edges.
    mark_edges.resize(edges.size(), 0);

    // -- Faces --
    // Faces are input triangles

    // Find flat edges
    for(CHAMedge& e : edges){
        if( findIF_flat_edge(e) ) e.type = CHAMedge_t::flat; 
        else                      e.type = CHAMedge_t::undet;

        if (e.inc_face.size() & 1){ def_interior = false; }
        if (e.inc_face.size() > 2){ manifold = false; }
    }

    #ifdef PLCC_VERBOSE_DEBUG
    disp_howManyFlatEdges();
    #endif

    // Create a face for each input triangle
    faces.resize( plc.numTriangles() );
    mark_faces.resize( faces.size(), 0);
    for(uint32_t i=0; i<plc.numTriangles(); i++) faces[i].triangle = i;

    // associate to each face (input triangle) its edges
    for(uint32_t ei=0; ei<edges.size(); ei++)
        for(uint32_t& t : edges[ei].inc_face) 
            faces[t].bounding_edges.push_back(ei);

    for(CHAMface& f : faces) orient_initial_triface_bnd(f);

    if(verbose) std::cout<<"chamfer-PLC initialization COMPLETED\n";
}

// -------------------------- //
// SEARCH ACUTE PLCC ELEMENTS //
// -------------------------- //

// Assumes that the PLCedge 'e' is one of the sides of the input-triangle 'ti'.
// Returns the index of the vertex of 'ti' different from 'e' endpoints.
uint32_t PLCc::inTri_opp_vrt(const CHAMedge& e, const uint32_t ti) const {
  const uint32_t* tv = plc.triangle_vertices.data() + (3*ti);

  assert((e.hasVertex(*tv) || e.hasVertex(*(tv+1))) && "assumption violated\n");
  
  uint32_t v = *tv; if( !e.hasVertex(v) ) return v;
  v = *(++tv);      if( !e.hasVertex(v) ) return v;
  v = *(++tv);      
  
  assert( !e.hasVertex(v) ); 
  
  return v;
}

// Assumes that the vertex 'v' is a vertex of the input-triangle 'ti'.
// Returns (through 'u[0]' , 'u[1]') the vertices 'ti' different from 'v'.
void PLCc::inTri_opp_edge(const uint32_t v, const uint32_t ti, 
                                        uint32_t& u1, uint32_t& u2) const{

    const uint32_t* tv = plc.triangle_vertices.data() + (3*ti);

    assert((v == tv[0] || v == tv[1] || v == tv[2]) && "assumption violated\n");

    if(v == tv[0]){ u1 = tv[1];  u2 = tv[2]; return; }
    if(v == tv[1]){ u1 = tv[2];  u2 = tv[0]; return; }
    if(v == tv[2]){ u1 = tv[0];  u2 = tv[1]; return; }

   assert(false); // SHOULD NEVER BEEN REACHED
}

// Returns TRUE if the two incident faces at 'e' forms an acute dihedral angle.
bool PLCc::findIF_acute_edge(const CHAMedge& e) const {

    if(e.inc_face.size() == 1) return false;
    if(e.inc_face.size() > 2) return true; 
    // Delaunay Refinement algorithm (tetmesh.h) does not support non-manyfold 
    // configurations. (Thus we chamfer edges in such cases)

    uint32_t u, v;
    const pointType* e0_pt = vertices[ e.ep[0] ];
    const pointType* e1_pt = vertices[ e.ep[1] ];
    for(size_t i=1; i<e.inc_face.size(); i++){
        for(size_t j=0; j<i; j++){
            u = inTri_opp_vrt(e, faces[ e.inc_face[i] ].triangle);
            v = inTri_opp_vrt(e, faces[ e.inc_face[j] ].triangle);
            if( isAcuteDihedral_exact(e0_pt, e1_pt, 
                                    vertices[u], vertices[v]) ) return true;
        }
    }

    return false;
}

// Returns TRUE if:
// 1) two non-flat edges incident at 'vi' forms an acute angle, OR
// an edge and an input triangle incident at 'vi' forms an acute angle. 
bool PLCc::findIF_acute_vrt(uint32_t vi, 
                            const std::vector<uint32_t>& vv_i, 
                            const std::vector<uint32_t>& vt_i    ) const {
    
    const pointType* vip = vertices[vi];

    // Check each edge-edge angles
    for(uint32_t j = 1; j < vv_i.size(); j++){
        const pointType* vjp = vertices[ vv_i[j] ];
        for(uint32_t k = 0; k < j; k++){
            // is exact and returns TRUE if acute at second argument
            if( isAcuteAngle( vjp, vip, vertices[ vv_i[k] ] ) ) return true;
        }
    }

    // Check angles formed by each edge and triangle sharing exactly a vertex.
    uint32_t t0, t1, t2, ep;
    for(uint32_t j = 0; j < vv_i.size(); j++){
        ep = vv_i[j];
        const pointType* epp = vertices[ vv_i[j] ];
        for(uint32_t k = 0; k < vt_i.size(); k++){
            const uint32_t* tv = plc.triangle_vertices.data() + vt_i[k] * 3 ;
            t0 = tv[0];  t1 = tv[1];  t2 = tv[2];

            // vi is a vertex of vt_i[k], 
            assert(vi == t0 || vi == t1 || vi == t2);
            // if ep is a vertex of vt_i[k] too
            // then <vi,ep> is an edge of vt_i[k]
            if( t0 == ep || t1 == ep || t2 == ep ) continue;
            
            if( t1 == vi ) std::swap(t0, t1); 
            else if( t2 == vi ) std::swap(t0, t2); 
            // now t0 is vi

            if(isAcuteAngle(epp, vip, vertices[t1], vertices[t2])) return true;
        }
    }    

  return false;
}

// Returns the length of the shortest segment <vi, vv_i[j]>
double PLCc::closest_vv_dist(uint32_t vi, 
                                const std::vector<uint32_t>& vv_i ) const {
    double sq_dist = vPtsSqDist(vi, vv_i[0]);
    double tmp_sq_dist;
    for(size_t j=1; j<vv_i.size(); j++){
        tmp_sq_dist = vPtsSqDist(vi, vv_i[j]); 
        if(tmp_sq_dist < sq_dist) sq_dist = tmp_sq_dist;
    }
    return sqrt( sq_dist );
}

// returns the shortest between the length of the segments that connect
// a vertex vi with the opposite sides of its incident triangles
double PLCc::closest_vt_dist(uint32_t vi, 
                                const std::vector<uint32_t>& vt_i ) const {
    uint32_t u1, u2;   inTri_opp_edge(vi, vt_i[0], u1, u2);
    double sq_dist = vSqVrtDistSeg(vi, u1, u2);
    double tmp_sq_dist;
    for(size_t j=0; j<vt_i.size(); j++){
        inTri_opp_edge(vi, vt_i[j], u1, u2);
        tmp_sq_dist = vSqVrtDistSeg(vi,u1,u2);
        if( tmp_sq_dist < sq_dist ) sq_dist = tmp_sq_dist;
    }
    return sqrt( sq_dist );
}

// Search for acute vertices and edges and computes chamfering cut distances.
// For acute vertices 'vrt_ch_dist' is 1/3 of the shortest incident edge length.
void PLCc::search_acute_angles(){

    if(vrt_ch_dist.size() != vertices.size()) 
        vrt_ch_dist.resize(vertices.size(), -1.0); // init

    // Classify not-flat edges between acute and not-acute.

    // Compute " v -> all v' connected with v by a not-flat edge " relation.
    std::vector< std::vector<uint32_t> > vv(vertices.size()); 
    for(CHAMedge& e : edges) if( !e.isFlat() ){
        vv[e.ep[0]].push_back(e.ep[1]);
        vv[e.ep[1]].push_back(e.ep[0]);

        if( findIF_acute_edge(e) ) e.type = CHAMedge_t::acute;
    }

    #ifdef PLCC_VERBOSE_DEBUG
    disp_howManyAcuteEdges();
    #endif

    // Force acutness of acute edge endpoints (necessary for edge chamfering)
    for(const CHAMedge& e : edges) if(e.isAcute()){ 
        uint32_t e0 = e.ep[0], e1 = e.ep[1];
        if(!is_acute_vrt(e0)) 
            vrt_ch_dist[e0] = get_vrt_ch_dist( e0, vv[e0], input_vt[e0] );
        if(!is_acute_vrt(e1)) 
            vrt_ch_dist[e1] = get_vrt_ch_dist( e1, vv[e1], input_vt[e1] );
    }

    // Classify acute vertices
    for(uint32_t vi=0; vi<vertices.size(); vi++) if( !is_acute_vrt(vi) ) {
        if( findIF_acute_vrt(vi, vv[vi], input_vt[vi]) ) 
            vrt_ch_dist[vi] = get_vrt_ch_dist( vi, vv[vi], input_vt[vi] );
    }

    // Chamfering distance must take in to account flat edge lenght,
    // otherwise flat edges may be "cutted out". 
    for(const CHAMedge& e : edges ) if( e.isFlat() ){
        uint32_t e0 = e.ep[0], e1 = e.ep[1];
        double e_len = sqrt( vPtsSqDist(e0, e1) );
        if( vrt_ch_dist[e0] * 3 > e_len ) vrt_ch_dist[e0] = e_len / 3;
        if( vrt_ch_dist[e1] * 3 > e_len ) vrt_ch_dist[e1] = e_len / 3;
    }

    #ifdef PLCC_VERBOSE_DEBUG
    disp_howManyAcuteVertices();
    #endif

    if(verbose) std::cout<<"Acute vertices and edges detection COMPLETED\n";
}

// ---------------------------------------- //
// CHAMFERING POINTS CREATION AND PLACEMENT //
// -------------------------------- ------- //

// Creates a new point 'p' on the segment <'v0','v1'> at distance 'd' from:
// - 'v0' if 'd0' = 0
// - 'v1' if 'd0' = 1
// The point 'p' is EXACTLY on the segment <'v0','v1'>.
// Returns the position of the new vertex in the vertices vector.
uint32_t PLCc::new_vrt_on_segment(uint32_t v0, uint32_t v1, 
                                            double d, uint32_t d0) {
    
    if(d0==1) std::swap(v0, v1);
    double l = sqrt( vPtsSqDist(v0, v1) );
    double t = d / l;

    assert( t > 0 && t < 1 && "cannot place LNC on segment: invalid t\n" );
  
    // new point p is such that p = (1-t) * v0 + t * v1
    pointType* p;
    uint32_t exp3D_i = (t < 0.5) ? v0 : v1; // closest explicit point
    const pointType* p0 = vertices[v0];
    const pointType* p1 = vertices[v1];

    assert(p0->isExplicit3D() && p1->isExplicit3D() && "Implicit endpoint\n");
    
    p = new implicitPoint_LNC( p0->toExplicit3D(), p1->toExplicit3D(), t);
    return add_vertex( p, exp3D_i);
}

bool is_BPT_representable(double xi, double eta) {
    return (xi > 0.0  &&  eta > 0.0  &&  xi + eta < 1.0);
}

// DEBUG
bool disp_non_representable_bpt_info(double xi0, double eta0,
                                     double xi1, double eta1,
                                     double t0, double t1,
                                     double cos_alpha,
                                     double teta, double psi, double psi_lim) {

    std::cout << "WARNING: [get_bpt_coeff] "
              << "Steiner point not representable as BPT.\n"
              << "Computed baricentric coordinates: \n"
              << "xi0 = " << xi0 << ", eta0 = " << eta0 << "\n"
              << "xi1 = " << xi1 << ", eta1 = " << eta1 << "\n"
              << "Input parameters: \n"
              << "t0 = " << t0 << ", t1 = " << t1 << "\n"
              << "Internal parameters: \n"
              << "cos_alpha = " << cos_alpha << "\n"
              << "teta = " << teta 
              << ", psi = " << psi << ", psi_lim = " << psi_lim << "\n";
    return false;
}

// A segment with 'L' (LNC) and 'B' (BPT) endpoints shoul be orthogonal to the 
// side <'R','P'> or <'R','Q'> which 'L' belongs to. 
// By the way, due to numerical rounding of BPT baricentric coordinates 'xi'
// and 'eta', the two angles at 'L' on the half-plane containing 'B' are not 
// exactly pi/2: one is acute and the other is obtuse. 
// If the acute angle is inside the chamfered face it may spoil output acuteness
// while if acute angle is outside it does not matter as it will be increased by
// the amplitude of the external angle of the other face incident at 'L'.
// A trivial solution is to perturbate 'B' in order to satisfy the condition.
bool correct_coeff(double& xi, double& eta, 
                   const pointType* P, const pointType* Q, const pointType* R,
                    const pointType* L, const pointType* E) {
    
    assert( P->isExplicit3D() );
    assert( Q->isExplicit3D() );
    assert( R->isExplicit3D() );
    assert( L->isLNC() );
    assert( L->toLNC().P().toExplicit3D() == R->toExplicit3D() || 
            L->toLNC().Q().toExplicit3D() == R->toExplicit3D()      );
    assert( E->isExplicit3D() );
    assert( E->toExplicit3D() == P->toExplicit3D() || 
            E->toExplicit3D() == Q->toExplicit3D()      );
    assert( L->toLNC().P().toExplicit3D() == E->toExplicit3D() || 
            L->toLNC().Q().toExplicit3D() == E->toExplicit3D()      );

    // vector3d B = P * xi + Q * eta + R * (1 - (xi+eta));
    pointType* B = new implicitPoint3D_BPT(*P,*Q,*R,xi,eta);

    bool success = true;
    bigfloat k = eta/xi;
    // Removed angle R-L-B must be acute so that its supplementar angle
    // (internal to chamfered face) will result non-acute.
    while( isAcuteAngle(E,L,B) ) {

        assert(isAcuteAngle(B,L,E));
        assert(!isAcuteAngle(B,L,R));
        assert(!isAcuteAngle(R,L,B));
        
        // try a perturbation reducing t0 or t1 while they remain positive
        xi = std::nextafter(xi, -DBL_MAX);
        if(xi <= 0) { success=false; break; }
        eta = (xi * k).get_d();
        if(eta <= 0) { success=false; break; }

        delete B;
        B = new implicitPoint3D_BPT(*P,*Q,*R,xi, eta);
    }

    assert(!isAcuteAngle(B,L,E));
    assert(!success || isAcuteAngle(B,L,R) || 
            pointType::dotProductSign3D(*B, *R, *L) == 0);
    assert(!success || isAcuteAngle(R,L,B) || 
            pointType::dotProductSign3D(*R, *B, *L) == 0);

    delete B;
    return success;
}
// interface for above function
bool correct_coeff(double& xi, double& eta, 
                    double t0, double t1,
                    const vector3d& OP, const vector3d& OQ, const vector3d& OR,
                    bool relative_to_Q) {
    
    const pointType* P = new explicitPoint3D(OP.c[0],OP.c[1],OP.c[2]);
    const pointType* Q = new explicitPoint3D(OQ.c[0],OQ.c[1],OQ.c[2]);
    const pointType* R = new explicitPoint3D(OR.c[0],OR.c[1],OR.c[2]);
    const pointType* L = relative_to_Q ? 
            new implicitPoint3D_LNC(R->toExplicit3D(), Q->toExplicit3D(), t1) :
            new implicitPoint3D_LNC(R->toExplicit3D(), P->toExplicit3D(), t0);
    const pointType* E = relative_to_Q ? 
                            new explicitPoint3D(OQ.c[0],OQ.c[1],OQ.c[2]) :
                            new explicitPoint3D(OP.c[0],OP.c[1],OP.c[2]);
    return correct_coeff(xi,eta, P,Q,R, L,E);
}

bool PLCc::get_bpt_coeff(double& xi0, double& eta0, double& xi1, double& eta1, 
                   const vector3d& Ov, const vector3d& Ou0, const vector3d& Ou1, 
                   double t0, double t1, const double zero_toll){

    double ang = getAngle(new explicitPoint3D(Ov.c[0],Ov.c[1],Ov.c[2]),
            new explicitPoint3D(Ou0.c[0],Ou0.c[1],Ou0.c[2]) ,
            new explicitPoint3D(Ou1.c[0],Ou1.c[1],Ou1.c[2]) );

    double cos_alpha = cosOfAngle_at(Ov, Ou0, Ou1); // cos at v

    double one_plus_cos_alpha = 1+cos_alpha;
    double one_min_cos_alpha = 1-cos_alpha;

    if(abs( one_plus_cos_alpha ) < zero_toll ||
       abs( one_min_cos_alpha ) < zero_toll    ){ 
        std::cout<<"WARNING: [get_bpt_coeff] "
                   "internal parameters below zero_toll.\n";
    }

    // Baricentric coordinates of BPTs depend on 3 different parameters:
    // - the angle at V through its cosine "cos_alpha", 
    // - the chamfering distance "d" = dist(V, LNC_i)
    // - a free parameter "k" (which is determined by the equal length request).
    // To have <LNC_0,BPT_0>, <BPT_0,BPT_1>, <BPT_1,LNC_1> of equal length,
    // the following relation between "k" and "d" must hold:
    // k = (d * sqrt(2) * sqrt(1 - cos_alpha))  /  
    //     (1 + sqrt(2) * sqrt(1 + cos_alpha)).

    // To short notation we use 2 derived parameter "teta" and "psi" instead of
    // "cos_aplta" and "d".

    double teta = 1.0 / one_plus_cos_alpha;
    double psi, psi_max;
    
    // Being psi = (k / d) * sqrt( (1 + cos_alpha) / (1 - cos_alpha ) )
    // and substituting the value of "k" that gives equal length:
    psi = 1.0  /  ( 1.0 + 1.0 / sqrt(one_plus_cos_alpha * 2.0) );

    // Geometrically, "k" represents the distance between BPT_i and LNC_i,
    // which cannot be greater that "d" itself, otherwise BPTs may be 
    // placed outside non-acute triangles.
    // When "k" = "d" we get the maximum distance that cannot be exceded.
    psi_max = sqrt( one_plus_cos_alpha / one_min_cos_alpha ); 
    // BTPs have to stay at most at distance d from triangle side, when 
    // psi is forced to psi_max bridge-edges have no more same length:
    // as uinque consequence shorter edges are generated.
    if(psi > psi_max) psi = psi_max; 
    
    // In case of non-acute angles at V, BPTs may stay outside triangle or
    // too close to triangle side opposite to V.
    // if teta * (t0+t1) < 1.0 THEN BTPs SHURELY stay inside the triangle.
    double psi_lim = -1.0;
    if( teta * (t0+t1) >= 1.0 ){ 
        
        assert(cos_alpha <= 0); // This happens for non-acute angles at V

        // BPTs may not stay inside triangle <Ov,Ou0,Ou1>
        double psi_lim_0 = (1.0-t0) / ( (teta-1.0)*t0 + teta*t1 );
        double psi_lim_1 = (1.0-t1) / ( (teta-1.0)*t1 + teta*t0 );
        // The value of psi_lim_0 correspond to a point aligned with LNC_0 and
        // BPT_0 on the triangle side opposite to V. Similarly for psi_lim_1.
        psi_lim = min(psi_lim_0, psi_lim_1) / 3;
        if( psi > psi_lim )psi = psi_lim;
    }
    
    double psi_teta = psi * teta;
    double term = 1.0 + psi * (teta - 1.0);

    xi0 = t0 * term;
    eta0 = t1 * psi_teta;
    xi1 = t0 * psi_teta;
    eta1 = t1 * term;

    if(!correct_coeff(xi0, eta0, t0, t1, Ou0, Ou1, Ov, 0)) {
        report_error_and_exit("cannot get non-acute angle at LNC\n");
    }
    if(!correct_coeff(xi1, eta1, t0, t1, Ou0, Ou1, Ov, 1)) {
        report_error_and_exit("cannot get non-acute angle at LNC\n");
    }

    // DEBUG
    if( !is_BPT_representable(xi0, eta0) || !is_BPT_representable(xi1, eta1) ){
        return disp_non_representable_bpt_info( xi0, eta0, xi1, eta1, t0, t1,
                                                cos_alpha, teta, psi, psi_lim );
    }

    assert(is_BPT_representable(xi0, eta0) && is_BPT_representable(xi1, eta1));

    return true;
}

// Creates 2 new implicit BPT points 'p0' and 'p1' inside the CHAMface whose
// index is 'fi'. New BPTs are placed "near" input vertex whose index is 'vi'
// and to compute their coordinates an LNCs 'u0' and 'u1' must be present on the 
// input triangle side <'vi','ou0'> and <'vi',ou1'> respectivelly.
// The point is EXACTLY inside the input triangle containing 'fi'.
uint32_t PLCc::new_vrts_in_inputTri(const uint32_t fi, const uint32_t vi, 
                                    const uint32_t u0, const uint32_t ou0, 
                                    const uint32_t u1, const uint32_t ou1){

    assert( vertices[u0]->isLNC() && vertices[u1]->isLNC() );
    assert( vertices[ou0]->isExplicit3D() );
    assert( vertices[ou1]->isExplicit3D() );
    assert( vertices[vi]->isExplicit3D() );

    double xi0, eta0; // coefficient of BPT relative to u0
    double xi1, eta1; // coefficient of BPT relative to u1
    const vector3d Ov( vertices[vi] );
    const vector3d Ou0( vertices[ou0] );
    const vector3d Ou1( vertices[ou1] );

    const explicitPoint3D& pv = vertices[vi]->toExplicit3D();
    const explicitPoint3D& pu0 = vertices[ou0]->toExplicit3D();
    const explicitPoint3D& pu1 = vertices[ou1]->toExplicit3D();

    double t0 = vertices[u0]->toLNC().T();
    if(vertices[u0]->toLNC().Q().toExplicit3D() == pv) t0 = 1-t0;
    double t1 = vertices[u1]->toLNC().T();
    if(vertices[u1]->toLNC().Q().toExplicit3D() == pv) t1 = 1-t1;

    if(!get_bpt_coeff(xi0,eta0,xi1,eta1, Ov,Ou0,Ou1,t0,t1)) return INVALID_BPT;

    pointType *p0 = new implicitPoint3D_BPT(pu0, pu1, pv, xi0, eta0);
    pointType *p1 = new implicitPoint3D_BPT(pu0, pu1, pv, xi1, eta1);
    
    add_vertex( p0, vi );
    return add_vertex( p1, vi );
}

// --------------------- //
// ORTHOGONAL CHAMFERING //
// --------------------- //

// splt edge 'ei' at distance 'd' from endpoint edeges[ei].ep[ep_i]
void PLCc::chamfer_edge_ep(const size_t ei, double d, const uint32_t ep_i){

    CHAMedge& e = edges[ei]; // e = <e0,e1>
    const uint32_t e0 = e.ep[0];

    #ifdef PLCC_VERBOSE_DEBUG
    std::cout<<"\nchamfering (d="<<d<<") ep "<<ep_i<<" of "; print_edge(ei);
    #endif

    const uint32_t Pt_i = new_vrt_on_segment(e.oep[0], e.oep[1], d, ep_i);

    assert( vPointInInnerSegment(Pt_i, e.ep[0], e.ep[1]) );
    
    // Update PLCedges
    CHAMedge new_e(e); // copy e
    e.ep[1] = Pt_i; // Update PLCedge e endpoints: <e0,e1> becomes <e0,Pt>
    new_e.ep[0] = Pt_i; // New edge is <Pt_i,e1>
    edges.push_back( new_e ); // <Pt,e1> 
    mark_edges.push_back(0);

    assert( ( e.ep[0]==e.oep[0] || 
                vPointInInnerSegment(e.ep[0], e.oep[0], e.ep[1]) ) &&
            ( e.ep[1] == e.oep[1] || 
                vPointInInnerSegment(e.ep[1], e.ep[0], e.oep[1]) ) );

    // Update type
    if(ep_i==1) edges.back().type = CHAMedge_t::junk;
    else        edges[ei].type = CHAMedge_t::junk;

    #ifdef PLCC_VERBOSE_DEBUG
    print_just_chamfered_edge((uint32_t)ei, ep_i);
    #endif

    // update connectivity (boundaries of incident faces)
    uint32_t new_ei = (uint32_t) edges.size()-1;
    for(uint32_t fi : edges[ei].inc_face){

        faces[fi].make_last((uint32_t)ei); // <e0,Pt> is the last on face boundary.
        if(edges[ faces[fi].bounding_edges.front() ].hasVertex(e0)){
            faces[fi].make_first((uint32_t)ei);
        }
        faces[fi].bounding_edges.push_back( new_ei );

        #ifdef PLCC_VERBOSE_DEBUG
        std::cout<<"\n"; print_face_edges(fi);
        #endif
    }
}

// Splits each edge that has at least one acute endpoint 'v'. A Steiner point
// is inserted at distance 'r' = min('epsilon', 'vrt_ch_dist(v)'/3) from 'v'.
// NOTE. A steiner points is not acute: it has only the two sub-edges as 
//       incident edges which forms two angles of "pi" amplitude.
void PLCc::chamfering_vrts(){
    size_t n_edges_before = edges.size();
    uint32_t v;  double r;
    for(size_t ei=0; ei<n_edges_before; ei++) if(!edges[ei].isIsolated()){

        v = edges[ ei ].ep[ 1 ];
        if( is_acute_vrt(v) ){ 
            r = vrt_ch_dist[v]; if(epsilon < r) r = epsilon;
            chamfer_edge_ep(ei, r, 1); // splits at distance r from endpoint 1
            // now edges[ei] is <e0,new_point>
        }
        
        v = edges[ ei ].ep[ 0 ];
        if( is_acute_vrt(v) ){ 
            r = vrt_ch_dist[v]; if(epsilon < r) r = epsilon;
            chamfer_edge_ep(ei, r, 0); // splits at distance r from endpoint 0
        }
    }
}

// safe chamfering
implicitPoint3D_BPT* move_BPT_toward_LNC(
                                pointType* bpt, pointType* lnc, double min_d){
    
    assert(bpt->isBPT() && lnc->isLNC());

    const implicitPoint3D_BPT* BPT = &(bpt->toBPT());
    const implicitPoint3D_LNC* LNC = &(lnc->toLNC());

    const explicitPoint3D* P = &(BPT->P().toExplicit3D());
    const explicitPoint3D* Q = &(BPT->Q().toExplicit3D());
    const explicitPoint3D* R = &(BPT->R().toExplicit3D());
    double old_xi, old_eta;
    get_baricentric_coords(old_xi, old_eta, bpt);
    // BPT = old_xi * P + old_eta * Q + (1 - old_xi - old_eta) * R

    const explicitPoint3D* A = &(LNC->P().toExplicit3D());
    const explicitPoint3D* B = &(LNC->Q().toExplicit3D());
    double t = LNC->T();
    // LNC = t * A + (1 - t) * B

    // By construction 'R' must be shared by both LNC and BPT, while
    // the other explicit point of the LNC must be equal to 'P' or 'Q'
    assert(*A == *R || *B == *R);
    assert(*A == *P || *A == *Q || *B == *P || *B == *Q);

    // To simplify the algoith we assume that when 't' = 1 we get the explicit
    // point of LNC corresponding to 'R', i.e. 'A' = 'R'. 
    // If this is not the case just reverse 't'.
    if(*B == *R) t = 1.0 - t;

    // Since one LNC explicit points is 'R' there are only two possibilities:
    // 1) the other LNC explicit point is 'P',
    // 2) the other LNC explicit point is 'Q'
    bool common_is_P = (*A == *P || *B == *P);

    assert(!isAcuteAngle((*A==*P) ? A : B, lnc, bpt));

    vector3d OP(P), OQ(Q), OR(R);
    double distsq_PR = OP.dist_sq(OR);
    double distsq_QR = OQ.dist_sq(OR);
    double dotprod_PR_QR = (OP-OR).dot(OQ-OR);

    double a = old_xi, b = t - old_eta; // when common is Q
    if(common_is_P){ a = t - old_xi;  b = old_eta; }

    double k = a*a * distsq_PR + b*b * distsq_QR - 2*a*b * dotprod_PR_QR;
    k = sqrt( abs(k) ); // avoid rounding errors when close to 0

    double c = min_d / k;
    double eta = t - c * b, xi = c * a; // when common is Q
    if(common_is_P){ xi = t - c * a;  eta = c * b; }

    if(!correct_coeff(xi, eta, P, Q, R, lnc, common_is_P ? P : Q)) {
        std::cout<<"[cham.cpp - move_BPT_toward_LNC()] "
                 <<"cannot get non-acute angle at LNC\n";
        exit(1);
    }

    assert( is_BPT_representable(xi, eta) );

    return new implicitPoint3D_BPT(*P, *Q, *R, xi, eta);
}

// If the face has at least an acute edge, all the face edges connecting 
// a BPT to a LNC have to be "normalized" to the same (minimum) length.
// NOTE. This step is necessary only in theory to guarantee the elimination
//       of all acute angles; practically it have no measurable effects 
//       as possible acute angle are close to "pi/". For that reason 
//       'safe_mode' is deactivated by default.
void PLCc::normalize_face(uint32_t fi){
    const std::vector<uint32_t>& fbnd = faces[fi].bounding_edges;
    
    // Compute minimum BPT-LNC endpoint edge length
    double sq_min_d = DBL_MAX, edge_len_sq;
    for(uint32_t bei : fbnd ) if( bei != EMPTY_PLACE ){
        if(is_edge_normalizable(bei)) sq_min_d = min(sq_min_d, eEdgeSqLen(bei)); 
    }

    // Normalize to 'sq_min_d' the length of BPT-LNC edges while keeping their
    // "direction" (line which they belong) untached.
    double min_d = sqrt(sq_min_d);
    uint32_t lnc_i, bpt_i;
    for(size_t bei : fbnd) if( bei != EMPTY_PLACE ){
        lnc_i = edges[bei].ep[0], bpt_i = edges[bei].ep[1];
        if(!vertices[lnc_i]->isLNC()){ 
            std::swap(lnc_i, bpt_i);
            if(!vertices[lnc_i]->isLNC()) continue;
        }
        if(vertices[bpt_i]->isBPT() && eEdgeSqLen(bei) > sq_min_d){ 
            vertices.push_back( 
                move_BPT_toward_LNC(vertices[bpt_i], vertices[lnc_i], min_d) );
            std::swap( vertices[bpt_i], vertices.back() );
            vertices.pop_back();
        }
    }

    assert( check_face(fi, true) );
}

//
void PLCc::chamfering_face(uint32_t fi){

    #ifdef PLCC_VERBOSE_DEBUG
    std::cout<<"\nchamfering "; print_face_edges(fi);
    #endif

    // Part 1:
    // each couple of edges incident at an acute vertex will be replaced by 3
    // new edges, i.e. the size of bounding_edges have to be rised by the
    // number of acute vertices on face boundary.

    // NOTE. since we added LNC vertices on edges incident at an acute vertex,
    //       it is impossible that exists an edge with both acute endpoints.

    // Be shure that first and last bounding edges do not share an acute vertex.
    uint32_t first_on_bnd = faces[fi].bounding_edges.front();
    uint32_t last_on_bnd = faces[fi].bounding_edges.back();
    uint32_t v0 = edges[ first_on_bnd ].commonVertex( edges[ last_on_bnd ] );
    if(is_acute_vrt( v0 )) faces[fi].make_first( last_on_bnd ); // by rotation

    size_t n_acute_vrts = 0;
    for(uint32_t& ei : faces[fi].bounding_edges) 
        if(has_acute_endpoint(ei)) n_acute_vrts++;
    n_acute_vrts = n_acute_vrts / 2; // each vertex have been counted two times. 

    if(n_acute_vrts == 0) return;

    faces[fi].bounding_edges.resize( 
        faces[fi].bounding_edges.size() + n_acute_vrts, EMPTY_PLACE);

    std::vector<uint32_t>& fbnd = faces[fi].bounding_edges;
    bool impossible = false;
    for(size_t i=0; i < fbnd.size()-1; i++) if(fbnd[i] != EMPTY_PLACE) {

        #ifdef PLCC_VERBOSE_DEBUG
        std::cout<<"curr "; print_edge(fbnd[i]);
        #endif

        uint32_t ei = fbnd[i],  nei = fbnd[i+1];
        uint32_t vi = edges[ei].commonVertex( edges[nei] );

        if( !is_acute_vrt(vi) ) continue;
        
        // replace <v1,vi> + <vi,v4> with <v1,v2> + <v2,v3> + <v3,v4>

        #ifdef PLCC_VERBOSE_DEBUG
        std::cout<<"\nvertex "<<vi<<" between\n"; print_edge(ei); 
        std::cout<<"and\n"; print_edge(nei); std::cout<<"is acute.\n\n";
        #endif

        // 1) create two new vertices (BPT points v2 and v3) inside the face
        uint32_t v1, v4, ov1, ov4;
        v1 = edges[ei].oppositeVertex(vi);
        ov1 = edges[ei].oppositeOriginalVertex(vi);
        v4 = edges[nei].oppositeVertex(vi);
        ov4 = edges[nei].oppositeOriginalVertex(vi);
        uint32_t last = new_vrts_in_inputTri(fi, vi, v1, ov1, v4, ov4);
        if(last == INVALID_BPT){ impossible = true; break;}

        // 2) create 3 new edges
        uint32_t v2 = last -1, v3 = last;
        uint32_t num_e = (uint32_t)edges.size()+3;
        edges.resize(num_e);
        edges[num_e-3].init_bridge_edge(v1, v2, vi, fi);
        edges[num_e-2].init_bridge_edge(v2, v3, vi, fi);
        edges[num_e-1].init_bridge_edge(v3, v4, vi, fi);
        mark_edges.resize(num_e, 0);

        #ifdef PLCC_VERBOSE_DEBUG
        std::cout<<"\nnew edges:\n"; 
        print_edge( (uint32_t) edges.size()-3 );
        print_edge( (uint32_t) edges.size()-2 );
        print_edge( (uint32_t) edges.size()-1 );
        #endif

        assert( check_bridge(vi, v1, v2, v3, v4) );

        // 3) update boundary
        edges[ei].removeIncidentFace(fi);
        edges[nei].removeIncidentFace(fi);
        uint32_t pne = (uint32_t) edges.size()-1;
        fbnd[i] = pne-2; fbnd[i+1] = pne-1;
        for(size_t k=fbnd.size()-1; k>i+2; k--) fbnd[k] = fbnd[k-1];
        fbnd[i+2] = pne;

        #ifdef PLCC_VERBOSE_DEBUG
        std::cout<<"\nupdated face bnd "; print_face_edges(fi);
        #endif

        assert( check_face(fi, true) );

        i += 2;
    }

    if(impossible){ 
        if(verbose) get_unchamferableFace_info(fi);
        report_error_and_exit("ERROR: BPT cannot be built\n.");
    }

    // Part 1/2 (optional): face normalization.
    // When 'safe_mode' chamfering with strong theoretical guarantees of not 
    // introducing new acute angles is performed. 
    if(safe_mode && has_acute_edge(fi)) normalize_face(fi);

    // Part 2:
    // each acute edge will be removed the face by creating a new edge on each
    // incident face connecting two appropriate bridge points.

    #ifdef PLCC_VERBOSE_DEBUG
    std::cout<<"\nChamfer acute edges\n\n";
    #endif

    // Compute the distance between BPT facing the same acute edge to 
    // estimate the size of the "hole" created by removing the acute edge.
    // This will contribute to determinate the LFS of the chamfered PLC.
    for(const CHAMedge& e : edges) if( e.isAcute() ) {
        
        assert(!e.isFlat() && !e.isJunk());
        if( e.inc_face.size() < 2 ) continue;
        
        // Collect incident face BPT vertices
        std::vector< vector<uint32_t> > bpt_f(e.inc_face.size());
        for(size_t i = 0; i < e.inc_face.size(); i++) {
            std::vector<uint32_t> fv;
            get_face_vertices(faces[e.inc_face[i]], fv);
            for(uint32_t v : fv) if(vertices[v]->isBPT()) bpt_f[i].push_back(v);
        }

        double min_dist_sq = DBL_MAX;
        for(size_t i = 0; i < e.inc_face.size()-1; i++)
            for(size_t j = i+1; j < e.inc_face.size(); j++) 
                for(uint32_t v : bpt_f[i]) for(uint32_t u : bpt_f[j])
                    min_dist_sq = min(min_dist_sq, vPtsSqDist(u, v));
        
        lfs = min(lfs, sqrt(min_dist_sq));
    }

    // Edge chamfering
    for(size_t i=0; i < fbnd.size(); i++ ) if( fbnd[i] != EMPTY_PLACE ){

        uint32_t ei = fbnd[i];

        #ifdef PLCC_VERBOSE_DEBUG
        std::cout<<"\ncurr "; print_edge(ei);
        #endif

        if( edges[ei].isAcute() ){

            assert( !edges[ei].isFlat() );

            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"is acute\n";
            #endif
            
            const CHAMface& f = faces[fi];
            size_t ni = i, pi = i;
            f.advance_on_pierced_bnd(ni);
            f.reverse_on_pierced_bnd(pi);
            uint32_t vn = edges[ fbnd[ni] ].notCommonVertex( edges[ei] );
            uint32_t vp = edges[ fbnd[pi] ].notCommonVertex( edges[ei] );

            edges.push_back( CHAMedge( vn, vp, fi) );
            mark_edges.push_back(0);

            edges[ fbnd[ni] ].inc_face.clear(); // incident only at 'fi'
            edges[ fbnd[pi] ].inc_face.clear(); // incident only at 'fi'
            edges[ ei ].removeIncidentFace(fi);
            faces[fi].bounding_edges[i] = (uint32_t) edges.size()-1;
            faces[fi].bounding_edges[ni] = EMPTY_PLACE;
            faces[fi].bounding_edges[pi] = EMPTY_PLACE;

            #ifdef PLCC_VERBOSE_DEBUG
            std::cout<<"\nupdated face bnd "; print_face_edges(fi);
            #endif

            assert( check_face(fi, true) );

            i = ni;
        }

    }

    faces[fi].remove_empty_places();

    #ifdef PLCC_VERBOSE_DEBUG
    std::cout<<"after removing..\n"; print_face_edges(fi);
    #endif

    assert( check_face(fi, true) );
}

// --------------- //
// CHAMFERING MAIN //
// --------------- //

void PLCc::chamfering(){

    assert( acute_edges_have_acute_ep() && checkup() );

    chamfering_vrts();

    if(verbose) std::cout<<"Vertex chamfering COMPLETED\n";

    assert( checkEdges_beforeJunkDeletion() );

    for(uint32_t fi=0; fi<faces.size(); fi++) chamfering_face(fi);

    // Delete junk edges (and all isoleted edges) from edges vector
    for(CHAMedge& e : edges) if( e.isJunk() ){ e.inc_face.clear(); }
    cleanUp_edges();

    if(verbose) std::cout<<"Face chamfering COMPLETED\n";

    assert( checkup() );
}

// --------------------------------------- //
// CHAMFERED PLC SIMPLIFICATION (Optional) //
// --------------------------------------- //

#define EMPTY_BRIDGE_PIECE UINT32_MAX

class bridge{
    public:
    uint32_t el, ec, er;    // Bridge edges before simplification
    double shortest_e_sq;   // Shortest bridge edge before simplification
    uint32_t s, sl, sr;     // Bridge edges after simplification

    bridge() :  el(EMPTY_BRIDGE_PIECE), 
                ec(EMPTY_BRIDGE_PIECE), 
                er(EMPTY_BRIDGE_PIECE), 
                shortest_e_sq(DBL_MAX),
                s(EMPTY_BRIDGE_PIECE), 
                sl(EMPTY_BRIDGE_PIECE), 
                sr(EMPTY_BRIDGE_PIECE) {}

    // Edge indexed 'edge_i' wrt PLCc edges vector, of endpoints 'pe0' and 'pe1'
    void set_bridge_edge(uint32_t edge_i,
                        const pointType* pe0, const pointType* pe1 ) {
        if( pe0->isLNC() && pe1->isBPT() ) el = edge_i;
        else if( pe0->isBPT() && pe1->isBPT() ) ec = edge_i;
        else {
            assert(pe0->isBPT() && pe1->isLNC() && 
                    "Invalid edge configuration on chamfered face\n" ); 

            er = edge_i;
        }
        double e_len_sq = vector3d(pe0).dist_sq(vector3d(pe1));
        if(e_len_sq < shortest_e_sq) shortest_e_sq = e_len_sq;
    }

    inline bool is_3bridge() const { 
        return el!=EMPTY_BRIDGE_PIECE && er!=EMPTY_BRIDGE_PIECE; 
    }
    inline bool is_2bridge() const { 
        return (el==EMPTY_BRIDGE_PIECE && er!= EMPTY_BRIDGE_PIECE) || 
                (el!=EMPTY_BRIDGE_PIECE && er== EMPTY_BRIDGE_PIECE); 
    }

    static inline bool bridgeSortFunction(const bridge& b1, const bridge& b2) { 
        if(b1.shortest_e_sq == b2.shortest_e_sq) return b1.ec < b2.ec;
        return (b1.shortest_e_sq < b2.shortest_e_sq); 
    }
};

inline std::ostream& operator<<(std::ostream& os, const bridge& b) {
    os <<"<";
    if(b.el != EMPTY_BRIDGE_PIECE) os << (b.el); else os << "empty";
    os << ", " << b.ec << ", " ;
    if(b.er != EMPTY_BRIDGE_PIECE) os << (b.er); else os << "empty";
    os << "> - short_len = " << b.shortest_e_sq;
    if( b.s != EMPTY_BRIDGE_PIECE) os << " => <" << b.s << ">";
    else if( b.sl != EMPTY_BRIDGE_PIECE) os << " => <" 
                                            << b.sl << ", " << b.sr << ">";
    return os;
}

void print_all_bridges(std::vector<bridge>& all_bridges, 
                        const std::vector<CHAMedge>& edges) {
    std::cout << "\n";
    for(const bridge& b : all_bridges) 
        std::cout << b <<" (ref vrt."<< edges[b.ec].loc_face_bridge_id <<")\n";
    std::cout << "\n";
}

// Edge 'ei' has two incident edges on the boundary of its (unique) incident
// face 'fi'. 'conn_edge' is one of this two incident edge that has exactly two 
// incident faces 'fi' and 'fj'. If 'conn_edge' exists meets 'ei' at vertex V.
// 'cons_edge' is the edge of 'fj', different from 'cons_edge' incident at V.
// NOTE: It is assumed that edges[ei] has a unique incident face.
uint32_t PLCc::get_cons_edge_on_adj_face(uint32_t ei){

    assert(edges[ei].inc_face.size() == 1);

    uint32_t conn_edge;
    uint32_t fi = edges[ei].inc_face[0];
    faces[fi].make_first(ei); // now 'ei' is fbnd[0]
    const std::vector<uint32_t>& fbnd = faces[fi].bounding_edges;
    conn_edge = fbnd[1];
    if(edges[conn_edge].inc_face.size() == 1) conn_edge = fbnd.back();
    if(edges[conn_edge].inc_face.size() != 2) return UINT32_MAX;
    // conn_edge found!
    uint32_t fj = edges[conn_edge].inc_face[0];
    if(fj == fi) fj = edges[conn_edge].inc_face[1];
    faces[fj].make_first(conn_edge);
    uint32_t cons_edge = faces[fj].bounding_edges[1];
    if(!edges[cons_edge].has_commonVertex(edges[ei])) 
        cons_edge = faces[fj].bounding_edges.back();
    
    assert(edges[cons_edge].has_commonVertex(edges[ei]));
    assert(edges[cons_edge].has_commonVertex(edges[conn_edge]));
    
    return cons_edge;
}

PLCc::bridge_simpl_t PLCc::is_3bridge_simplifiable( uint32_t next_bel, 
                                uint32_t bel, uint32_t bec, uint32_t ber, 
                                uint32_t next_ber) const {
    
    bridge_simpl_t simpl = bridge_simpl_type::no;
    const CHAMedge& el = edges[bel];
    const CHAMedge& ec = edges[bec];
    const CHAMedge& er = edges[ber];
    const pointType* p1 = vertices[ el.notCommonVertex(ec) ];
    const pointType* p2 = vertices[ er.notCommonVertex(ec) ];
    if(next_bel != UINT32_MAX) {
        const pointType* p0 = vertices[ edges[next_bel].notCommonVertex(el) ];
        if(!isAcuteAngle(p0, p1, p2)) simpl = bridge_simpl_type::left;
    }
    if(next_ber != UINT32_MAX) {
        const pointType* p3 = vertices[ edges[next_ber].notCommonVertex(er) ];
        if(!isAcuteAngle(p1, p2, p3)){ 
            if(simpl == bridge_simpl_type::no) simpl = bridge_simpl_type::right;
            else simpl = bridge_simpl_type::both;
        }
    }
    return simpl;
}

void PLCc::simplify_3bridge(uint32_t bel, uint32_t bec, uint32_t ber) {

    // Search for edge incident at 'bel' on the adjacent face. 
    uint32_t next_bel = get_cons_edge_on_adj_face(bel);
    // Search for edge incident at 'ber' on the adjacent face. 
    uint32_t next_ber = get_cons_edge_on_adj_face(ber);

    // print_3bridge_details(bel, bec, ber, next_bel, next_ber, false); // DEBUG

    bridge_simpl_t simpl = 
        is_3bridge_simplifiable(next_bel, bel, bec, ber, next_ber);

    if(simpl == bridge_simpl_type::no) return;

    uint32_t v0 =  edges[next_bel].notCommonVertex(edges[bel]);
    uint32_t v1 =  edges[bel].notCommonVertex(edges[bec]);
    uint32_t v2 =  edges[ber].notCommonVertex(edges[bec]);
    uint32_t v3 =  edges[next_ber].notCommonVertex(edges[ber]);
    
    uint32_t fi = edges[bec].inc_face[0];

    if(simpl == bridge_simpl_type::both) {

        // Replace 'bel', 'bec' and 'ber' with a unique new edge <'v1','v2'>.
        edges.push_back( CHAMedge(v1,v2,fi) ); mark_edges.push_back(0);
        size_t new_fbnd_size = faces[fi].bounding_edges.size()-2 ;
        faces[ fi ].make_last(ber);
        if(faces[fi].bounding_edges[new_fbnd_size] != bec) 
            faces[fi].make_last(bel);
        faces[fi].bounding_edges[ new_fbnd_size-1 ] = (uint32_t)edges.size()-1; 
        faces[fi].bounding_edges.resize(new_fbnd_size);
        edges[bel].isolate();
        edges[bec].isolate();
        edges[ber].isolate();
        return;
    }

    assert(simpl==bridge_simpl_type::left || simpl==bridge_simpl_type::right);

    // Initialize for 'left' simplification, if 'right' correct.
    uint32_t be_rem = bel, new_ep = v1, ext_ep = v0;  
    if(simpl==bridge_simpl_type::right){ be_rem = ber; new_ep = v2; ext_ep = v3;}

    uint32_t c = edges[bec].notCommonVertex(edges[be_rem]);

    // Replace 'bel' and 'bec' with a unique new edge <'new_ep', 'c'> if no
    // acute angles are introduced at 'new_ep' or 'c'.
                            
    if(!isAcuteAngle(vertices[v1],vertices[c],vertices[v2]) &&
        !isAcuteAngle(vertices[ext_ep],vertices[new_ep],vertices[c]) ){

        edges.push_back( CHAMedge(new_ep,c,fi) ); mark_edges.push_back(0);
        faces[ fi ].make_last(be_rem);
        size_t new_fbnd_size = faces[fi].bounding_edges.size()-1;
        if( faces[fi].bounding_edges[ new_fbnd_size-1 ] != bec){ 
            faces[ fi ].make_last(bec);
            assert(faces[fi].bounding_edges[ new_fbnd_size-1 ] == be_rem);
        }
        faces[fi].bounding_edges[ new_fbnd_size-1 ] = (uint32_t)edges.size()-1; 
        faces[fi].bounding_edges.resize(new_fbnd_size);
        edges[bec].isolate();
        edges[be_rem].isolate();
    }

}

void PLCc::simplify_2bridge(uint32_t bes, uint32_t bec) {

    // Search for edge incident at 'bes' on the adjacent face. 
    uint32_t next_bes = get_cons_edge_on_adj_face(bes);

    // print_2bridge_details(bec, bes, next_bes, false); // DEBUG
    
    if(next_bes == UINT32_MAX) return;

    uint32_t v0 = edges[next_bes].notCommonVertex(edges[bes]);
    uint32_t v1 = edges[next_bes].commonVertex(edges[bes]);
    uint32_t v2 = edges[bec].notCommonVertex(edges[bes]);

    if(!isAcuteAngle(vertices[v0], vertices[v1], vertices[v2])) {
        uint32_t fi = edges[bec].inc_face[0];
        edges.push_back( CHAMedge(v1,v2,fi) ); mark_edges.push_back(0);
        faces[ fi ].make_last(bes);
        size_t new_fbnd_size = faces[fi].bounding_edges.size()-1 ;
        if( faces[fi].bounding_edges[ new_fbnd_size-1 ] != bec){ 
            faces[ fi ].make_last(bec);
            assert(faces[fi].bounding_edges[ new_fbnd_size-1 ] == bes);
        }
        faces[fi].bounding_edges[ new_fbnd_size-1 ] = (uint32_t)edges.size()-1; 
        faces[fi].bounding_edges.resize(new_fbnd_size);
        edges[bec].isolate();
        edges[bes].isolate();
    }
}

void PLCc::chamfered_plc_simplification(){

    size_t old_num_edges = edges.size();

    // Collect bridge-edges
    std::vector<bridge> bridges;
    for(CHAMface& f : faces) {

        // On each face boundary, bridges (made of one, or two or three 
        // consecutive) edges are always alternated by one non-bridge edge.

        // Rotate bounding_edges vector untill there are no bridges divided by
        // the begin and the end of the vector.
        if(edges[f.bounding_edges[0]].isBridgeEdge() &&
            edges[f.bounding_edges.back()].isBridgeEdge()   ){
                do {
                    f.make_first( f.bounding_edges.back() );
                } while( edges[ f.bounding_edges.back() ].isBridgeEdge() );
        }        

        assert(!edges[f.bounding_edges[0]].isBridgeEdge() ||
               !edges[f.bounding_edges.back()].isBridgeEdge());

        const std::vector<uint32_t>& fbnd = f.bounding_edges;

        for(size_t i=0; i<fbnd.size(); i++) if(edges[ fbnd[i] ].isBridgeEdge()){
            bridges.push_back( bridge() );
            bridge& b = bridges.back();
            do{ 
                const CHAMedge& e = edges[ fbnd[i] ];
                b.set_bridge_edge(fbnd[i], vertices[e.ep[0]],vertices[e.ep[1]]);
                i++;
            } while( edges[ fbnd[i] ].same_bridge(edges[ fbnd[i-1] ])) ;

            if(edges[ fbnd[i] ].isBridgeEdge()) i--;
        }
    }

    if(verbose){ 
        uint32_t nbrs = 0;  
        for(const CHAMedge& e : edges) if(e.isBridgeEdge()) nbrs++;
        if(nbrs > 0) std::cout<<"There are "<< nbrs << " bridge-edges and "
                              << bridges.size() << " bridges.\n";
    }

    #ifdef PLCC_VERBOSE_DEBUG
    print_all_bridge_edges();
    #endif 
    
    assert( check_bridge_edges() );

    // Order bridges by increasing shortest length
    std::sort( bridges.begin(), bridges.end(), bridge::bridgeSortFunction );

    #ifdef PLCC_VERBOSE_DEBUG
    print_all_bridges(bridges, edges);
    #endif 

    // Simplify bridges if acute angles are not introducted
    for(size_t i=0; i<bridges.size(); i++){
        bridge& b = bridges[i];
        if( b.is_3bridge() ) simplify_3bridge(b.el, b.ec, b.er);
        else if(b.is_2bridge()) 
            simplify_2bridge( (b.er == EMPTY_BRIDGE_PIECE) ? b.el : b.er, b.ec);
    }

    cleanUp_edges(); // Remove superfluous edges

    if(verbose) std::cout << "Simplication COMPLETED: " 
                        << old_num_edges - edges.size() << " edges removed.\n";
}


// ------------------- //
// FINAL TRIANGULATION //
// ------------------- //

//
bool isConvexVrt(const pointType *p, const pointType *vrt, const pointType *q, int* ref_nv_sign)
{
    // DEBUG
    // std::cout<<"tri norm vect signs ("<<ref_nv_sign[0]<<","<<
    //                                     ref_nv_sign[1]<<","<<
    //                                     ref_nv_sign[2]<<")\n";
    //  std::cout<<"loc norm vect signs ("<<pointType::orient2Dyz(*p,*vrt,*q)<<","<<
    //                                     pointType::orient2Dzx(*p,*vrt,*q)<<","<<
    //                                     pointType::orient2Dxy(*p,*vrt,*q)<<")\n";
    
    if(pointType::orient2Dyz(*p,*vrt,*q) != ref_nv_sign[0]) return false;
    if(pointType::orient2Dzx(*p,*vrt,*q) != ref_nv_sign[1]) return false;
    if(pointType::orient2Dxy(*p,*vrt,*q) != ref_nv_sign[2]) return false;
    return true;
}

//
bool pointInTriangle(const pointType* pt, const pointType* v0, const pointType* v1, const pointType* v2){
    // Returns true if the coplanar point 'p' is either in the inner area of
    // 't' or on its border. Undetermined if p and t are not coplanar.
	if (pointType::pointInSegment(*pt, *v0, *v1)) return true;
	if (pointType::pointInSegment(*pt, *v1, *v2)) return true;
	if (pointType::pointInSegment(*pt, *v2, *v0)) return true;
	return pointType::pointInInnerTriangle(*pt, *v0, *v1, *v2);
}

// HEAR CLIPPING TRIANGULATION ALGORITHM
// Triangulation of a 2D polygon.
// Assumptions: the polygon is simple (no holes, closed, no self-intersections).
void PLCc::hear_clipping(uint32_t fi, std::vector<uint32_t>& out_tri_fv_list) const {

    // std::cout<<"HEAR CLIPPING\n"; // DEBUG

    const uint32_t* in_tri = plc.triangle_vertices.data() + (faces[fi].triangle*3);
    const pointType* ov0 = vertices[*in_tri];
    const pointType* ov1 = vertices[*(in_tri+1)];
    const pointType* ov2 = vertices[*(in_tri+2)];
    int signs_normVect_tri[] = {pointType::orient2Dyz(*ov0,*ov1,*ov2),
                                pointType::orient2Dzx(*ov0,*ov1,*ov2),
                                pointType::orient2Dxy(*ov0,*ov1,*ov2) };
            
    // Find first ear vertex. Create triangle. Remove ear vertex from list.
    // Do this while number of vertices (fv.size()) > 2.

    std::vector<uint32_t> fv; get_face_vertices(faces[fi], fv);
    size_t jl, jc = 0, jr;
    size_t fv_last = fv.size()-1;
    while(fv.size() > 2){

        if(jc == 0){ jl = fv_last; jr = 1; }
        else{ 
            if(jc == fv_last){ jl = jc-1; jr = 0; }
            else{ jl = jc-1; jr = jc+1; }
        }

        // DEBUG
        // std::cout<<"\njl = "<<jl<<", jc = "<<jc<<", jr = "<<jr<<"\n";
        // std::cout<<" fv["<<jl<<"] = "<<fv[jl]<<
        //            ", fv["<<jc<<"] = "<<fv[jc]<<
        //            ", fv["<<jr<<"] = "<<fv[jr]<<"\n";

        const pointType* vl = vertices[ fv[jl] ];
        const pointType* vc = vertices[ fv[jc] ];
        const pointType* vr = vertices[ fv[jr] ];

        bool is_convex = isConvexVrt(vl,vc,vr,signs_normVect_tri);

        // DEBUG
        // std::cout<<fv[jc]<<" is ";
        // if(!is_convex) std::cout<<"NOT ";
        // std::cout<<"convex vertex.\n";

        bool is_ear = true;

        if(is_convex){

            size_t jtest = 0;
            if(jr != fv_last) jtest = jr+1;

            // DEBUG
            // std::cout<<"jtest = "<<jtest<<"; fv["<<jtest<<"] = "<<fv[jtest]<<"\n";

            while( jtest != jl && is_ear ){
                const pointType* vt = vertices[ fv[jtest] ];
                is_ear = !pointInTriangle(vt,vl,vc,vr);

                // DEBUG
                // std::cout<<fv[jtest]<<" is ";
                // if(is_ear) std::cout<<"NOT ";
                // std::cout<<"inside <"<<fv[jl]<<","<<fv[jc]<<","<<fv[jr]<<">\n";

                if(jtest==fv_last) jtest = 0; else jtest++;

                // DEBUG
                // if(is_ear) std::cout<<"next jtest = "<<jtest<<"; fv["<<jtest<<"] = "<<fv[jtest]<<"\n";
            }

        }
        else{ is_ear = false; }

        if(is_ear){

            // DEBUG
            // std::cout<<"detraching triangle <"<<fv[jl]<<","<<fv[jc]<<","<<fv[jr]<<">\n";

            out_tri_fv_list.insert(out_tri_fv_list.end(), {fv[jl],fv[jc],fv[jr]} );
            fv.erase(std::remove(fv.begin(),fv.end(),fv[jc]));
            fv_last = fv.size()-1;
            jc = 0;
        }
        else{
            if(jc==fv_last) jc = 0; else jc++;
        }

    }
}

// 
void PLCc::get_triangles(std::vector<uint32_t>& tri_fv) const {

    for(size_t fi = 0; fi < faces.size(); fi++) {

        // std::cout<<"\ntriangualting "; print_face_edges(fi); // DEBUG

        if(faces[fi].bounding_edges.size() == 3){
            std::vector<uint32_t> fv;
            get_face_vertices(faces[fi], fv);
            tri_fv.insert(tri_fv.end(), {fv[0],fv[1],fv[2]} );
        }
        else{ hear_clipping((uint32_t)fi, tri_fv); }

    }
}

// Complmentar
void PLCc::get_complementar_tri(const std::vector<uint32_t>& out_tri, std::vector<uint32_t>& compl_tri){
    // make half-edges
    std::vector< std::vector<uint32_t> > he; // ep0, ep1, tri, occurences for each edge
    uint32_t v0,v1,v2,e0,e1;
    for(uint32_t i=0; i<out_tri.size()/3; i++){
        v0 = out_tri[3*i  ];
        v1 = out_tri[3*i+1];
        v2 = out_tri[3*i+2];
        e0 = v0; e1 = v1; if(v0 > v1) {e0 = v1; e1 = v0;} he.push_back({e0,e1,i,1});
        e0 = v1; e1 = v2; if(v1 > v2) {e0 = v2; e1 = v1;} he.push_back({e0,e1,i,1});
        e0 = v2; e1 = v0; if(v2 > v0) {e0 = v0; e1 = v2;} he.push_back({e0,e1,i,1});
    }
    sort(he.begin(), he.end(),
        [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b){ return (a[0] < b[0] || (a[0]==b[0] && a[1]<b[1])); } );

    for(size_t i=0; i<he.size()-1; i++){
        if(he[i][0] == he[i+1][0] && he[i][1] == he[i+1][1]){
            he[i+1][3] += he[i][3];
            he[i][3] = 0;
        }
    }

    he.erase(std::remove_if(he.begin(), he.end(), [](const std::vector<uint32_t> &a){ return (a[3] != 1); }), he.end());
    
    uint32_t r0, r1;
    for(size_t i=0; i<he.size(); i++){
        e0 = he[i][0]; e1 = he[i][1]; 

        if(vertices[e0]->isExplicit3D() || vertices[e1]->isExplicit3D() ) continue;

        r0 = ref_exp3D_vrt[e0]; r1 = ref_exp3D_vrt[e1];
        if( r0 == r1 ) {
            compl_tri.insert( compl_tri.end(), {e0, e1, r0} ); // orientation maybe have to be changed
        }
        else{
            if(e0 > e1) std::swap(e0, e1);
            compl_tri.insert( compl_tri.end(), {e0, e1, r0} );
            compl_tri.insert( compl_tri.end(), {r0, r1, e1} );
            // orientation maybe have to be changed
        }
    }
}
