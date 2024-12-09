#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#endif

#define EXIT_ON_THRESHOLD_NUMVERTICES //if (V.size() > 100000) { std::cout << "Limit num vertices reached!\nEXITING\n"; exit(11); }
 
// #define USE_TETGEN

// #define DISP_PROGRESS

#include <iostream>
#include <fstream>
#include "cdt.h"
#include "inputPLC.h"
#include "tetmesh.h"
#ifdef USE_TETGEN
#include "tetgen_interface.cpp"
#endif
#include <chrono>
#include "logger.h"
#include "getRSS.h"

using namespace std;

#ifdef ISVISUALSTUDIO
#ifndef NDEBUG
#define DEBUG
#endif
#endif

#include "cham.h"

typedef std::chrono::steady_clock chrono_clock;

// Export the data structure to optimizer,
// eturns true if the input surface enclses a volume, i.e. input edges have
// an even number of incident input triangles. 
bool chamferPLC(inputPLC& _plc,
	double _epsilon,
	std::vector<genericPoint*>& _vertices,
	std::vector<uint32_t>& _conn_vertices,
	std::vector<std::vector<uint32_t>>& _faces,
	bool simplify_cham_plc, bool print_surf,
	bool verbose) {
	
	PLCc* cut_plc = new PLCc(_plc, _epsilon, true, verbose);

	if (simplify_cham_plc) {
		size_t num_edges = cut_plc->edges.size();
		cut_plc->chamfered_plc_simplification();
		// assert( cut_plc->checkup() ); // DEBUG
		if (verbose) std::cout << "Chamfered PLC simplication COMPLETED: " << num_edges - cut_plc->edges.size() << " edges removed.\n";
	}

	std::vector<uint32_t> out_tri;
	cut_plc->get_triangles(out_tri); // computes a trinagulation of cut_plc (without modifing it) and saves into out_tris triangles indices wrt cut_plc->vertices
	if (verbose) std::cout << "Chamfered PLC triangulation COMPLETED\n\n";

	std::vector<uint32_t> used_vertex(cut_plc->vertices.size(), UINT32_MAX);
	for (size_t i = 0; i < cut_plc->n_in_vrts; i++) used_vertex[i] = 1; // all explicit points
	for (size_t i = 0; i < out_tri.size(); i++)	used_vertex[ out_tri[i] ] = 1;

	uint32_t idx = 0;
	for (size_t i = 0; i < used_vertex.size(); i++) if (used_vertex[i] != UINT32_MAX) {
		used_vertex[i] = idx++;
		_vertices.push_back(cut_plc->vertices[i]);
		_conn_vertices.push_back(cut_plc->ref_exp3D_vrt[i]); // explicit points do not change indexing
	}

	_faces.resize(out_tri.size() / 3);	
	uint32_t v0, v1, v2;
	for (size_t i = 0; i < out_tri.size() / 3; i++) {
		v0 = used_vertex[ out_tri[i * 3    ] ];
		v1 = used_vertex[ out_tri[i * 3 + 1] ];
		v2 = used_vertex[ out_tri[i * 3 + 2] ];
		_faces[i].assign({v0,v1,v2});
	}

	// cut_plc->save_rebuilded_input_after_chamfering(out_tri, "all_tris_chamf.off"); // DEBUG
	// cut_plc->saveFaces(); // save polygonal faces
	if(print_surf) cut_plc->saveTriangles(out_tri, "chamfered_plc.off"); // save triangulated faces

	return cut_plc->input_plc_defines_interior();
}

// createSteinerCDT
// 
// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!

TetMesh* createSteinerCDT(inputPLC& plc, double& min_PLC_dist, bool produce_output, bool verbose =false) {

	// Build a delaunay tetrahedrization of the vertices
	TetMesh* tin = new TetMesh;
	tin->init_vertices(plc.coordinates.data(), plc.numVertices());
	tin->tetrahedrize();

	if (verbose) printf("DT of the vertices built\n");

	// Build a structured PLC linked to the Delaunay tetrahedrization
	PLCx Steiner_plc(*tin, plc.triangle_vertices.data(), plc.numTriangles());
	Steiner_plc.segmentRecovery_HSi(!verbose);
	Steiner_plc.faceRecovery(!verbose);
	std::vector<bool> constr_tri_asCorners;
	Steiner_plc.markInnerTets_andGetConstrFaces(constr_tri_asCorners);
	//tin->saveConstrTrisToOFF("constrainedFaces.off", constr_tri_asCorners);
	if(min_PLC_dist != 0.0) min_PLC_dist = tin->compute_closest_features_dist(constr_tri_asCorners); // minimum distance between cdt elemnts contrained on the PLC

	return tin;
}

//
bool get_vrts_and_tris_for_cdt(Tetrahedrization& mesh, std::vector<uint32_t>& ref_vrts){
	
	// exctract triangles of the optimized mesh that are part of the chamfered surface (deltri)
	std::vector<uint32_t> constr_tri;
	size_t n_constr_tri = 0;
	uint32_t v0,v1,v2;
	for(const TetFace* tri : mesh.faces()) if(tri->deltri != NULL){
		if(tri->t1()==NULL || tri->t2()==NULL) continue; // do not use bounding box boundary faces
		v0 = (uint32_t) tri->deltri->v0()->getIndex();
		v1 = (uint32_t) tri->deltri->v1()->getIndex();
		v2 = (uint32_t) tri->deltri->v2()->getIndex();
		constr_tri.insert(constr_tri.end(), {v0,v1,v2});
		n_constr_tri++;
	}

	// make half-edges wrt constr_tri
    std::vector< std::vector<uint32_t> > be; // ep0, ep1, tri, occurences for each edge
    uint32_t e0,e1;
    for(uint32_t i=0; i<n_constr_tri; i++){
        v0 = constr_tri[3*i  ];
        v1 = constr_tri[3*i+1];
        v2 = constr_tri[3*i+2];
        e0 = v0; e1 = v1; if(v0 > v1) {e0 = v1; e1 = v0;} be.push_back({e0,e1,i,1});
        e0 = v1; e1 = v2; if(v1 > v2) {e0 = v2; e1 = v1;} be.push_back({e0,e1,i,1});
        e0 = v2; e1 = v0; if(v2 > v0) {e0 = v0; e1 = v2;} be.push_back({e0,e1,i,1});
    }
    sort(be.begin(), be.end(),
        [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b){ return (a[0] < b[0] || (a[0]==b[0] && a[1]<b[1])); } );

	// keep half-edges that occour only one time, i.e. those on the boundary of constr_tri
    for(size_t i=0; i<be.size()-1; i++){
        if(be[i][0] == be[i+1][0] && be[i][1] == be[i+1][1]){
            be[i+1][3] += be[i][3];
            be[i][3] = 0;
        }
    }
    be.erase(std::remove_if(be.begin(), be.end(), [](const std::vector<uint32_t> &a){ return (a[3] != 1); }), be.end());

	// make the vertex-edge relation of boundary edges
	std::vector< std::vector<uint32_t> > vbe_rel( mesh.num_vertices() );
	for(size_t i=0; i<be.size(); i++){
		vbe_rel[ be[i][0] ].push_back(i);
		vbe_rel[ be[i][1] ].push_back(i);
	}

	// complete ref_vrts by assigning to vertices of the boundary of constr_tri
	// introduced by Delauny refinement algorithm a ref_vrt to connect with while
	// closure triangles are created. 
	size_t n_cham_vrts = ref_vrts.size();
	ref_vrts.resize(mesh.num_vertices(), UINT32_MAX);

	for(size_t i=0; i<vbe_rel.size(); i++) if(ref_vrts[i]==UINT32_MAX) assert(vbe_rel[i].empty() || vbe_rel[i].size()==2); // DEBUG
	
	std::vector<uint32_t> queue;
	for(size_t vi=0; vi<n_cham_vrts; vi++) if(ref_vrts[vi]!=UINT32_MAX){
		// Consider each boundary edge incident at vi (at least 2)
		for(size_t j=0; j<vbe_rel[vi].size(); j++){
			uint32_t bei = vbe_rel[vi][j];
			e0 = be[bei][0]; e1 = be[bei][1]; // edge endpoints
			v0 = vi; // starting vertex
			// v0 is a vertex s.t. ref_vrt[v0] != UINT32_MAX, 
			// each boundary edge E = <b0,b1> of the chamfered PLC has a ref_vrt for both endpoint,
			// during Delaunay refinement new vertices may be inserted on E splitting it
			// producing a partition of E, each of these new vertices on E have exactly two incident
			// edges. 
			// Starting from b0 = v0 we wanto to rach b1 tracking all these new vertices using a queue to
			// assign them the correct ref_vrt.
			assert( ref_vrts[v0] != UINT32_MAX );
			v1 = e0; if(v1 == v0) v1 = e1;
			uint32_t curr_edge, last_edge = bei;
			while( ref_vrts[v1] == UINT32_MAX ){
				queue.push_back(v1);
				v0 = v1;
				assert(vbe_rel[v1].size() == 2);
				curr_edge = vbe_rel[v1][0]; if(curr_edge == last_edge) curr_edge = vbe_rel[v1][1];
				v1 = be[curr_edge][0]; if(v0 == v1) v1 = be[curr_edge][1];
				last_edge = curr_edge;
			}
			if( !queue.empty() ){
				v0 = vi;
				uint32_t refv = (ref_vrts[v0] < ref_vrts[v1]) ? ref_vrts[v0] : ref_vrts[v1]; // we have to choose one, the conventional rule is given by lower index
				for(uint32_t q : queue) ref_vrts[q] = refv;
				queue.clear();
			}
		}
	}

	// build the comlementary triangles of the chamfered optimized surface wrt input
    std::vector<uint32_t> compl_tri;
    uint32_t r0, r1;
    for(size_t i=0; i<be.size(); i++){
		e0 = be[i][0]; e1 = be[i][1]; 
		if(mesh.vrts()[e0]->getPoint()->isExplicit3D() || mesh.vrts()[e1]->getPoint()->isExplicit3D() ) continue;
		r0 = ref_vrts[e0]; r1 = ref_vrts[e1];
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

    // for(size_t i=0; i<compl_tri.size()/3; i++){ // DEBUG
    //     uint32_t a = compl_tri[3*i];
    //     uint32_t b = compl_tri[3*i+1];
    //     uint32_t c = compl_tri[3*i+2];
    //     assert( !vCollinear(a,b,c) );
    // }

	// DEBUG
	for(uint32_t i=0; i<mesh.num_vertices(); i++) {
		if( i != (uint32_t)mesh.vrts()[i]->getIndex() ){ 
			std::cout<<"# tot. vrts = "<<mesh.num_vertices()<<", # cham_vrts = "<<ref_vrts.size()<<"\n";
			std::cout<<"mismatching i = "<<i<<", vrt[i]_index = "<< (uint32_t)mesh.vrts()[i]->getIndex()<<"\n";
			exit(0);
		}
	}

	constr_tri.insert(constr_tri.end(), compl_tri.begin(), compl_tri.end());
	assert( (constr_tri.size() % 3) == 0 );

    FILE* fp = fopen("input_rebuild.off", "w");
    fprintf(fp, "OFF\n%u %u 0\n", mesh.num_vertices(), (uint32_t) constr_tri.size() / 3);

    for(uint32_t i=0; i<mesh.num_vertices(); i++) {
		assert( i == (uint32_t)mesh.vrts()[i]->getIndex() );
        double x, y, z;
        mesh.vrts()[i]->getPoint()->getApproxXYZCoordinates(x, y, z);
        fprintf(fp, "%f %f %f\n", x, y, z);
    }

    for(size_t i=0; i<constr_tri.size()/3; i++) {
        std::vector<uint32_t> v;
        v.assign(constr_tri.begin() + 3*i, constr_tri.begin() +3*i +3);
        fprintf(fp, "3 %u %u %u \n", v[0], v[1], v[2]);
    }
        
 	fclose(fp);

	return true;
}

// bool isTetInternal(Tetrahedron* t, TetMesh* cdt) {
// 	double ccc[3];
// 	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
// 	getTetBarycenter(v, ccc);
// 	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
// 	cdt->vertices.push_back(&p);
// 	static uint64_t tet = 0;
// 	tet = cdt->searchTetrahedron(tet, cdt->numVertices() - 1);
// 	cdt->vertices.pop_back();
// 	return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
// }

// Add by Lorenzo -- 06 Nov 2024 (old version above)
// mark as internal only tets that have an internal barycenter and all internal vertices.
bool isTetInternal(Tetrahedron* t, TetMesh* cdt) {
	double ccc[3];
	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
	getTetBarycenter(v, ccc);
	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
	
	static uint64_t tet = 0;
	tet = cdt->searchTetrahedron(tet, &p);
	//return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
	if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;

	t->mark<1>(); // TMP

	if (!(t->v0()->isMarked<0>())){
		tet = cdt->searchTetrahedron(tet, v[0]);
		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
	}

	if (!(t->v1()->isMarked<0>())){
		tet = cdt->searchTetrahedron(tet, v[1]);
		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
	}

	if (!(t->v2()->isMarked<0>())){
		tet = cdt->searchTetrahedron(tet, v[2]);
		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
	}

	if (!(t->v3()->isMarked<0>())){
		tet = cdt->searchTetrahedron(tet, v[3]);
		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
	}

	t->unmark<1>(); // TMP

	return true;
}

void smoothVertices(Tetrahedrization& mesh, uint32_t n_input_vrts, double optim_ratio){
	for (TetVertex* v : mesh.vrts()){ v->unmark<0>(); v->unmark<1>(); }; // Reset
	for (TetFace* f : mesh.faces()) if (f->isInterface()) {
		f->v0()->mark<1>();
		f->v1()->mark<1>();
		f->v2()->mark<1>();
	}
	// Tag explicit points whose index is greater than n_input_vrts: 
	// they are Steiner points that do not belong to the input surface
	for(uint32_t i=n_input_vrts; i<mesh.num_vertices(); i++){
		if(mesh.vrts()[i]->getPoint()->isExplicit3D()) mesh.vrts()[i]->mark<0>();
	}
	
	std::vector<double> print_coords; // TMP
	// The candidate optimal posotion for each tagged vertex is given 
	// by the meah of the coordinates of its neighbour vertices that belong to the surface.
	for (TetVertex* v : mesh.vrts()) if(v->isMarked<0>() && v->isMarked<1>()){
		TetVertices vv;  v->VV(vv);
		uint32_t n_constr_neighs = 0;
		double x_virt=0.0, y_virt=0.0, z_virt=0.0;
		for (TetVertex* w : vv) if(!w->isMarked<0>()){
			n_constr_neighs++;
			double x,y,z;  w->getPoint()->getApproxXYZCoordinates(x,y,z);
			x_virt += x;  y_virt += y;  z_virt += z;
		}
		// NOTE. what happens to movable vertices sorrounded only by movable vertices??
		// move the vertex to the new position only if the qulity of its incident tets
		// remain acceptable, if it is not the case try to get new potion closer to the original one.
		if(n_constr_neighs != 0){
			double den = (double) n_constr_neighs;
			x_virt /= den;  y_virt /= den;  z_virt /= den;
			explicitPoint3D* origin_v = (explicitPoint3D*) &v->getPoint()->toExplicit3D();
			double x,y,z; origin_v->getApproxXYZCoordinates(x,y,z);

			uint32_t n_it = 0;
			bool vt_quality_conserved = true;
			double t=1;
			while(n_it < 4){
				explicitPoint3D* virt_v = new explicitPoint3D( (1-t)*x + t*x_virt, (1-t)*y + t*y_virt, (1-t)*z + t*z_virt );
				v->setPoint(virt_v);
				Tetrahedra vt;  v->VT(vt);
				for(Tetrahedron* t : vt){
					if( !mesh.isTetPositive(t) || mesh.computeTetCost(t,false) > optim_ratio ){ 
						v->setPoint(origin_v);
						vt_quality_conserved = false; 
						delete virt_v;
						break;
					}
				}
				if(vt_quality_conserved){ 
					// std::cout<<"SMOOTHED("<<n_it<<")\n";
					// std::cout<<"ori: "<<origin_v->X()<<" "<<origin_v->Y()<<" "<<origin_v->Z()<<"\n";
					// std::cout<<"new: "<<virt_v->X()<<" "<<virt_v->Y()<<" "<<virt_v->Z()<<"\n";
					// std::cout<<"displ: "<<sqrt(vector3d(origin_v).dist_sq(vector3d(virt_v)))<<"\n\n";
					// delete origin_v; WHY MAKE THE PROGRAM CRASH??
					// TMP
					// if(sqrt(vector3d(origin_v).dist_sq(vector3d(virt_v))) > 0){
					// 	print_coords.push_back(virt_v->X());
					// 	print_coords.push_back(virt_v->Y());
					// 	print_coords.push_back(virt_v->Z());
					// }
					break;
				}
				t *= 0.5;
				vt_quality_conserved = true;
				n_it++;
			}
		}

	}

	// TMP
	if(!print_coords.empty()){
		FILE* fp = fopen("Mvrts.off", "w");
        fprintf(fp, "OFF\n%zu 0 0\n", print_coords.size()/3);
        for(size_t i=0; i<print_coords.size()/3; i++)
        fprintf(fp, "%f %f %f\n", print_coords[3*i], print_coords[3*i+1], print_coords[3*i+2]);
		fclose(fp);
	}

	for (TetVertex* v : mesh.vrts()){ v->unmark<0>(); v->unmark<1>();} // Reset
}

void markInternalTets(Tetrahedrization& mesh, TetMesh* cdt) {
	// mark all mesh vertices belonging to the PLC
	for (TetVertex* v : mesh.vrts()) v->unmark<0>();
	//for (PLC_Segment* s : mesh.get_PLCsegs()) { s->v0()->mark<0>(); s->v1()->mark<0>(); }	
	for (PLC_Face* f : mesh.get_PLCfaces()){
		for(DelEdge* e : f->getEdges() ) { e->v0()->mark<0>(); e->v1()->mark<0>(); }
	}
	
	//for (Tetrahedron* t : mesh.tets()) t->is_internal = isTetInternal(t, cdt);

	// Smarter method - exploit adjacencies to make location in CDT faster
	// for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
	for (Tetrahedron* t : mesh.tets()){ t->unmark<0>(); t->unmark<1>(); } // TMP
	Tetrahedra todo;
	todo.reserve(mesh.tets().size());
	Tetrahedron *s, *t = mesh.tets().front();
	todo.push_back(t); t->mark<0>();

	while (!todo.empty()) {
		t = todo.back();
		t->is_internal = isTetInternal(t, cdt);
		todo.pop_back();
		s = t->t0(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
		s = t->t1(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
		s = t->t2(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
		s = t->t3(); if (s != NULL && !s->isMarked<0>()) { todo.push_back(s); s->mark<0>(); }
	}



	for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
	for (TetVertex* v : mesh.vrts()) v->unmark<0>();
}

inline void logFinalStats(Tetrahedrization& mesh, uint64_t ms, bool input_encloses_vol){
	logUInteger((uint32_t)mesh.num_vertices());
	logUInteger((uint32_t)mesh.num_tetrahedra());
	logUInteger((uint32_t)ms);
	logDouble(mesh.minEdgeLength());
	
	double maxEneIN, maxEneEX;
	mesh.maxTetEnergy(maxEneIN, maxEneEX);
	logDouble(maxEneIN); 
	if(input_encloses_vol) logDouble(maxEneEX); else logEmpty();
	
	double minFAIN, maxFAIN, minFAEX, maxFAEX, minDAIN, maxDAIN, minDAEX, maxDAEX;
	mesh.minMaxTetAngle(minFAIN, maxFAIN, minFAEX, maxFAEX, minDAIN, maxDAIN, minDAEX, maxDAEX);
	logDouble(minFAIN); 
	logDouble(maxFAIN); 
	if(input_encloses_vol) logDouble(minFAEX); else logEmpty();
	if(input_encloses_vol) logDouble(maxFAEX); else logEmpty();
	logDouble(minDAIN); 
	logDouble(maxDAIN); 
	if(input_encloses_vol) logDouble(minDAEX); else logEmpty();
	if(input_encloses_vol) logDouble(maxDAEX); else logEmpty();
}


int main(int argc, char* argv[])
{
	initFPU();

	char filename[2048];
#ifndef DEBUG
	if (argc < 2) {
		std::cout << "Mesher - Create a well-shaped tetrahedral mesh out of a triangulated OFF file.\n";
		std::cout << "USAGE: ./delmesher [-v][-l][-o][-d 8] filename.off\n";
		std::cout << "OPTIONS:\n";
		std::cout << "\t[-a] -> extract from the Delaunay refined mesh a vaild input and gives it to a CDT algorithm that produces a conforming and quasi-optimal volume mesh (output see below) \n";
		std::cout << "\t[-s] -> produce as output the constrained surface of the quasi-optimal-CDT (see below)\n";
		std::cout << "\t[-i] -> (forces -a activation) if the input encloses a volume, uses only internal points and contrained faces as input for the quasi-optimal-CDT.\n";
		std::cout << "\t[-v] -> verbose mode\n";
		std::cout << "\t[-l] -> logging mode\n";
		std::cout << "\t[-d exp (pos. integer)] -> avoid Delaunay Refinement if the min pt.s dist after chamfering is < 10^-exp\n";
		std::cout << "\t[-b] -> computes the lower bound for minimum distance between mesh elements (point - segment - triangle), otherwise only point-point distances are computed.\n";
		std::cout << "\t[-o] -> produces outputs related to Delaunay refinement algorithm (see below)\n";
		std::cout << "\t[-c] -> produces an intermediate output (see below)\n";
		std::cout << "OUTPUT:\n";
		std::cout << "\t when [-a] is activated produces a volumetric (DRCDT_mesh.tet) meshe.\n";
		std::cout << "\t when [-s] is activated produces a surface (DRCDT_plcfaces.off) meshes.\n";
		std::cout << "\t when [-o] is activated produces a volumetric (DR_mesh.tet) and a surface (DR_plcfaces.off) meshes.\n";
		std::cout << "\t when [-c] is activated produces a surface (chamfered_plc.off) mesh of the chamfered plc (intermediate construction).\n";
		std::cout << "RETURNS:\n";
		std::cout << "\t0 when the whole execution terminates correctly (also when an iperror occours)\n";
		std::cout << "\t10 when option -d is activated and min dist. is violated\n";
		std::cout << std::endl;
		return 0;
	}
	else strcpy(filename, argv[1]);
#else
	strcpy(filename, "..\\Input_file\\acute\\cup_fixed_fixed.off");
	//strcpy(filename, "..\\Input_file\\subdcube.off");
	//char filename[2048] = "..\\Input_file\\twocubes.off";

#endif

	// Manage options
	std::string options = "";
	uint32_t min_dist_exp = UINT32_MAX;

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'd') { min_dist_exp = atoi(argv[++i]); continue; }
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i]) + 1);

	bool log_mode = (options.find('l') != std::string::npos);
	bool verbose_mode = (options.find('v') != std::string::npos);
	bool output_cdt = (options.find('a') != std::string::npos);
	bool only_interior_output_cdt = (options.find('i') != std::string::npos);
	if(only_interior_output_cdt) output_cdt = true;
	bool produce_outcdt_surf = (options.find('s') != std::string::npos);
	bool produce_cahm_surf = (options.find('c') != std::string::npos);
	bool produce_DR_output = (options.find('o') != std::string::npos);
	bool lowBnd_onMeshDist = (options.find('b') != std::string::npos);

	double optim_ratio = 2.0;

	if (log_mode) startLogging(filename);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, verbose_mode);

	if (!verbose_mode) std::cout<<"\ninput_file: "<<filename<<"\n";
	if (log_mode) advance_ProcessLogging("load_input");

#ifdef USE_TETGEN
	Tetrahedrization mesh;
	mesh.initWithTetgen(plc.numVertices(), plc.coordinates.data(), plc.numTriangles(), plc.triangle_vertices.data(), true, false);
	for (Tetrahedron* t : mesh.tets()) t->is_internal = true;
	mesh.printReport();
	std::cout << std::endl;
#else

	TetMesh tin;

	std::vector<genericPoint*> chamf_vertices;
	std::vector<uint32_t> to_close_ref_vrts;
	std::vector<std::vector<uint32_t>> chamf_faces;
	double epsilon = plc.bbDiag() / 1000.0;	// Use bounding-box-diagonal/1000 as chamfering distance

	bool simplify_chamferd_plc = true; // try to remove uncessary edges while keeping non-acute angles
	bool input_encloses_vol = chamferPLC(plc, epsilon, chamf_vertices, to_close_ref_vrts, chamf_faces, simplify_chamferd_plc, produce_cahm_surf, verbose_mode);

	if (log_mode) advance_ProcessLogging("chamfering");

	tin.init_vertices(chamf_vertices);

	std::cout << "Delaunizing vertices...\n";
	tin.addBoundingBoxVertices();
	tin.tetrahedrize();

	if (log_mode) advance_ProcessLogging("Del_vertices");

	Tetrahedrization mesh;
	chrono_clock::time_point time_point = chrono_clock::now();

	// Copy the DT to the new structure
	double closest_dist = mesh.initFromVerticesAndTets(tin.vertices, tin.tet_node);
	std::cout << "Distance of closest points relative to bb diagonal: " << closest_dist << "\n";

	TetMesh *cdt = NULL;  
	if(lowBnd_onMeshDist){
		double min_PLC_dist = DBL_MAX;
		cdt = createSteinerCDT(plc, min_PLC_dist, false); // needed to compute lower bound on generic mesh elements distances
		double mesh_BBox_len = euclideanDistance(mesh.vrts().back(), mesh.vrts()[mesh.num_vertices()-8]);
		min_PLC_dist = min_PLC_dist / (3.0 * mesh_BBox_len); // the lower bound is 1/3 * min_PLC_dist normalized wrt mesh bounding box diagonal
		std::cout << "Distance of closest generic mesh elements, relative to bb diagonal: " << min_PLC_dist << "\n";
		closest_dist = min(closest_dist, min_PLC_dist);
	}

	if(log_mode) logDouble(closest_dist);

	if (min_dist_exp != UINT32_MAX){
		double min_dist = std::pow(10.0, (double)min_dist_exp * (-1.0));
		if (closest_dist < min_dist){ 
			// ip_error("Closest points are too close. Optimization would produce too many tets!\nEXITING\n");
			if(log_mode) finishLogging();
			std::cout<<"\nPROGRAM ABORTED: Closest points are too close ( < "<< min_dist <<")\n\n\n";
			exit(10);
		}
	}
	
	if (log_mode) advance_ProcessLogging("DelRef_init");

	// Add PLC faces - at this stage faces are just collections of input triangles
	std::cout << "Adding PLC faces...\n";
	mesh.addPLCFaces(chamf_faces);

	if (log_mode) advance_ProcessLogging("add_PLCfaces");

	// Split all missing segments while no more remain.
	// This does not modify the triangles in PLC faces, but updates their edges with
	// the subedges being computed.
	std::cout << "Recovering segments...\n";
	mesh.recoverAllSegments();

	if (log_mode) advance_ProcessLogging("rec_seg");

	// Create a local 2D Delaunay triangulation for each PLC face
	std::cout << "Delaunizing faces...\n";
	mesh.delaunizePLCFaces();

	if (log_mode) advance_ProcessLogging("Del_faces");

	// Face recovery
	std::cout << "Recovering faces...\n";
	mesh.recoverAllFaces();

	if (log_mode) advance_ProcessLogging("rec_faces");

	// Tet optimization
	std::cout << "Optimizing tets...\n";
	mesh.optimizeTets(optim_ratio, false, true);

	if (log_mode) advance_ProcessLogging("tet_optim");

	// Creates a .off file with all optimized vertices and
	// an optimized triangulation of the input plc, to be used as input for cdt
	get_vrts_and_tris_for_cdt(mesh, to_close_ref_vrts);

	// Remove external tets after chamfering
	// In case the input edges has even number of incident input triangles,
	// i.e. the input surface defines an "interior",
	// we create a CDT of the input surface to decide wherever an "otimized mesh" tetrahedron 
	// is internal/external wrt the input surface.
	if(input_encloses_vol){
		std::cout<<"Input encloses a volume\n";
		if(cdt==NULL){ 
			double zero=0.0; 
			cdt = createSteinerCDT(plc, zero, false);
		}
		markInternalTets(mesh, cdt);

		if (log_mode) advance_ProcessLogging("IntExt_class");
	}
	else{ for (Tetrahedron* t : mesh.tets()) t->is_internal = true; }

	// smoothing vertices
	// smoothVertices(mesh, plc.numVertices()+8, optim_ratio);
	
	chrono_clock::time_point now = chrono_clock::now();
	uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	std::cout << "Elapsed time (ms): " << ms << "\n";
	uint64_t bmem = getPeakRSS();
	std::cout << "Peak memory RSS (byte): " << bmem << "\n";

	if (log_mode) {
		logFinalStats(mesh, ms, input_encloses_vol);
		advance_ProcessLogging("COMPLETED");
		finishLogging();
	}
	else {
		mesh.printReport(input_encloses_vol);
		std::cout << std::endl;
	}

	if(produce_DR_output){
		mesh.saveTET("DR_mesh.tet");
		mesh.saveOFFInterface("DR_plcfaces.off");
	}
	
	std::cout << "Execution correctly COMPLETED.\n\n\n";

#endif

	return 0;
}
