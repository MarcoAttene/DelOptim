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
// #include "tetmesh.h"
#ifdef USE_TETGEN
#include "tetgen_interface.cpp"
#endif
#include <chrono>
#include "logger.h"
#include "getRSS.h"
#include "optim_cdt.hpp"

using namespace std;

#ifdef ISVISUALSTUDIO
#ifndef NDEBUG
#define DEBUG
#endif
#endif

#include "cham.h"

typedef std::chrono::steady_clock chrono_clock;

// -------------------------- //
// -- CHAMFERING INTERFACE -- //
//--------------------------- //

// Export the data structure to optimizer throught _vertices, _conn_vertices, _faces.
// Returns true if the input surface enclses a volume,
// i.e. input edges have an even number of incident input triangles. 
bool chamferPLC(inputPLC& _plc, double _epsilon,
	std::vector<genericPoint*>& _vertices,
	std::vector<uint32_t>& _conn_vertices,
	std::vector<std::vector<uint32_t>>& _faces,
	bool& is_manifold,
	bool simplify_cham_plc, bool print_surf, bool verbose) {
	
	PLCc* cut_plc = new PLCc(_plc, _epsilon, verbose);

	if (simplify_cham_plc) {
		// Try to blend consecutive non-necesary edges
		size_t num_edges = cut_plc->edges.size();
		cut_plc->chamfered_plc_simplification();
		// assert( cut_plc->checkup() ); // DEBUG
		if (verbose) std::cout << "Chamfered PLC simplication COMPLETED: " << num_edges - cut_plc->edges.size() << " edges removed.\n";
	}

	// Fills out_tri with triples of vertices indices defining a triangulation of the chamfered PLC.
	std::vector<uint32_t> out_tri;
	cut_plc->get_triangles(out_tri); // computes a trinagulation of cut_plc (without modifing it) and saves into out_tris triangles indices wrt cut_plc->vertices
	if (verbose) std::cout << "Chamfered PLC triangulation COMPLETED\n\n";

	// Collect data to export the triangulated chamfered PLC to the optimizer
	std::vector<uint32_t> uv(cut_plc->vertices.size(), UINT32_MAX); // Used Vertex 
	for (size_t i = 0; i < cut_plc->n_in_vrts; i++) uv[i] = 1; // all explicit points
	for (size_t i = 0; i < out_tri.size(); i++)	uv[ out_tri[i] ] = 1; // uv used as marker

	uint32_t idx = 0;
	for (size_t i = 0; i < uv.size(); i++) if (uv[i] != UINT32_MAX) {
		uv[i] = idx++; // now uv stores new indexing
		_vertices.push_back(cut_plc->vertices[i]);
		_conn_vertices.push_back(cut_plc->ref_exp3D_vrt[i]); 
		// explicit points are all used and stored at the beginning of 
		// the vertices vector, so they do not change indexing.
	}

	_faces.resize(out_tri.size() / 3);	
	for (size_t i = 0; i < out_tri.size() / 3; i++) {
		const uint32_t* otv = out_tri.data() + i*3;
		_faces[i].assign({uv[ *otv ], uv[ *(otv+1) ], uv[ *(otv+2) ]});
	}

	// Produce some intermediate output mesh and save .off files 
	//cut_plc->save_rebuilded_input_after_chamfering(out_tri, "all_tris_chamf.off"); // Save chamfered input + complementar triangles to rebuild the input PLC
	// cut_plc->saveFaces(); // save polygonal faces
	if(print_surf) cut_plc->saveTriangles(out_tri, "chamfered_plc.off"); // save triangulated faces

	is_manifold = cut_plc->input_plc_is_manifold();
	return cut_plc->input_plc_defines_interior();
}

// ------------- //
// CDT interface //
// ------------- //

// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!
TetMesh* createSteinerCDT(inputPLC& plc, bool min_PLC_dist, bool produce_output, uint32_t& nct, bool verbose =false) {

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
	if (produce_output) tin->saveConstrTrisToOFF("constrainedFaces.off", constr_tri_asCorners);
	nct = tin->countConstrTris(constr_tri_asCorners);
	if (min_PLC_dist) tin->set_min_inputPLC_dist(constr_tri_asCorners); // minimum distance between cdt elemnts contrained on the PLC

	return tin;
}

// OPTIONAL: Compute the minum distance between any two mesh elements (vertices, edges, triangles)
//			 A cdt is needed to make this computation efficient.
double get_lower_bound_on_delRef_mesh_min_dist(inputPLC& plc, const Tetrahedrization& mesh, double closest_dist, TetMesh* cdt, uint32_t& nct){
	cdt = createSteinerCDT(plc, true, false, nct); // needed to compute lower bound on generic mesh elements distances
	double mesh_BBox_len = euclideanDistance(mesh.vrts().back(), mesh.vrts()[mesh.num_vertices()-8]);
	double min_PLC_dist = cdt->get_min_inputPLC_dist() / (3.0 * mesh_BBox_len); // the lower bound is 1/3 * min_PLC_dist normalized wrt mesh bounding box diagonal
	// std::cout << "Distance of closest generic mesh elements, relative to bb diagonal: " << min_PLC_dist << "\n";
	return min(closest_dist, min_PLC_dist);
}

bool isTetInternal(Tetrahedron* t, TetMesh* cdt) {
	double ccc[3];
	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
	getTetBarycenter(v, ccc);
	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
	cdt->vertices.push_back(&p);
	static uint64_t tet = 0;
	tet = cdt->searchTetrahedron(tet, cdt->numVertices() - 1);
	cdt->vertices.pop_back();
	return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
}

// Add by Lorenzo -- 06 Nov 2024 (old version above)
// mark as internal only tets that have an internal barycenter and all internal vertices.
// bool isTetInternal(Tetrahedron* t, TetMesh* cdt) {
// 	double ccc[3];
// 	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
// 	getTetBarycenter(v, ccc);
// 	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
	
// 	static uint64_t tet = 0;
// 	tet = cdt->searchTetrahedron(tet, &p);
// 	//return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
// 	if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;

// 	t->mark<1>(); // TMP

// 	if (!(t->v0()->isMarked<0>())){
// 		tet = cdt->searchTetrahedron(tet, v[0]);
// 		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
// 	}

// 	if (!(t->v1()->isMarked<0>())){
// 		tet = cdt->searchTetrahedron(tet, v[1]);
// 		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
// 	}

// 	if (!(t->v2()->isMarked<0>())){
// 		tet = cdt->searchTetrahedron(tet, v[2]);
// 		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
// 	}

// 	if (!(t->v3()->isMarked<0>())){
// 		tet = cdt->searchTetrahedron(tet, v[3]);
// 		if (cdt->mark_tetrahedra[tet>>2] != DT_IN) return false;
// 	}

// 	t->unmark<1>(); // TMP

// 	return true;
// }

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

void make_histogram(Tetrahedrization& mesh, const char* name){
	std::vector<uint32_t> fah, dah;
	mesh.makeFaceAngleHistogram(fah, 12);
	mesh.makeDihedralAngleHistogram(dah, 12);
	std::cout << name << " FACE ANGLE HISTOGRAM\n";
	for (auto x : fah) std::cout << x << "\t";
	std::cout << "\n";
	std::cout << name << " DIHEDRAL ANGLE HISTOGRAM\n";
	for (auto x : dah) std::cout << x << "\t";
	std::cout << "\n";
}

inline void logFinalStats(Tetrahedrization& mesh, uint64_t ms, bool input_encloses_vol){
	//logUInteger((uint32_t)mesh.num_vertices());
	//logUInteger((uint32_t)mesh.num_tetrahedra());
	//logUInteger((uint32_t)ms);
	//logDouble(mesh.minEdgeLength());
	
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

	double avMinPlan, avMaxPlan, avMinDihed, avMaxDihed;
	mesh.averageAngles(avMinPlan, avMaxPlan, avMinDihed, avMaxDihed);
	logDouble(avMinPlan);
	logDouble(avMaxPlan);
	logDouble(avMinDihed);
	logDouble(avMaxDihed);
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
		std::cout << "\t[-m max_vrt] -> set to max_vrt the maximum number of vertices of the mesh during the Delaunay refinement stage.\n";
		std::cout << "\t[-a] -> extract from the Delaunay refined mesh a vaild input and gives it to a CDT algorithm that produces a conforming and quasi-optimal volume mesh (output see below) \n";
		std::cout << "\t[-s] -> produce as output the constrained surface of the quasi-optimal-CDT (see below)\n";
		std::cout << "\t[-i] -> (forces -a activation) if the input encloses a volume, uses only internal points and constrained faces as input for the quasi-optimal-CDT.\n";
		std::cout << "\t[-v] -> verbose mode\n";
		std::cout << "\t[-l] -> logging mode\n";
		std::cout << "\t[-d exp (pos. integer)] -> avoid Delaunay Refinement if the min pt.s dist after chamfering is < 10^-exp\n";
		std::cout << "\t[-b] -> computes the lower bound for minimum distance between mesh elements (point - segment - triangle), otherwise only point-point distances are computed.\n";
		std::cout << "\t[-o] -> produces outputs related to Delaunay refinement algorithm (see below)\n";
		std::cout << "\t[-c] -> produces an intermediate output (see below)\n";
		std::cout << "OUTPUT:\n";
		std::cout << "\t when [-a] is activated produces a volumetric (DRCDT_mesh.tet) mesh.\n";
		std::cout << "\t when [-s] is activated produces a surface (constrainedFaces.off) meshes.\n";
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
	uint32_t max_vrts = UINT32_MAX;

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'd') { min_dist_exp = atoi(argv[++i]); continue; }
			if (argv[i][1] == 'm') { max_vrts = atoi(argv[++i]); continue; }
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

	chrono_clock::time_point time_point = chrono_clock::now();

	if (log_mode) startLogging(filename);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, verbose_mode);

	chrono_clock::time_point now = chrono_clock::now();
	uint64_t time_init = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (!verbose_mode) std::cout<<"\ninput_file: "<<filename<<"\n";
	if (log_mode){ 
		logUInteger(plc.numVertices());
		logUInteger(plc.numTriangles());
		logDouble(time_init);
		advance_ProcessLogging("load_input");
	}

#ifdef USE_TETGEN
	Tetrahedrization mesh;
	mesh.initWithTetgen(plc.numVertices(), plc.coordinates.data(), plc.numTriangles(), plc.triangle_vertices.data(), true, false);
	for (Tetrahedron* t : mesh.tets()) t->is_internal = true;
	mesh.printReport();
	std::cout << std::endl;
#else

	// Chamfering of the input PLC
	double epsilon = plc.bbDiag() / 1000.0;	// Use bounding-box-diagonal/1000 as chamfering distance
	bool simplify_chamferd_plc = true; // try to remove uncessary edges while keeping non-acute angles
	std::vector<genericPoint*> chamf_vertices;
	std::vector<uint32_t> to_close_ref_vrts;
	std::vector<std::vector<uint32_t>> chamf_faces;
	bool is_manifold;
	uint32_t input_encloses_vol = chamferPLC(plc, epsilon, chamf_vertices, to_close_ref_vrts, chamf_faces, is_manifold, simplify_chamferd_plc, produce_cahm_surf, verbose_mode);

	now = chrono_clock::now();
	uint64_t time_cham = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (log_mode){ 
		logBoolean(is_manifold);
		logBoolean(!input_encloses_vol);
		logUInteger((uint32_t)chamf_vertices.size());
		logUInteger((uint32_t)chamf_faces.size());
		logDouble(time_cham);
		advance_ProcessLogging("chamfering");
	}

	TetMesh tin; // Dealunay Tetrahedrization of the set of vertices of the chamfered PLC
	tin.init_vertices(chamf_vertices);

	if (verbose_mode) std::cout << "Delaunizing vertices...\n";
	tin.addBoundingBoxVertices(); // Adds to tin.vertices 8 new vertices defining a bounding box, by default dist = 1.0
	tin.tetrahedrize();

	if (log_mode) advance_ProcessLogging("Del_vertices");

	Tetrahedrization mesh; // Delaunay Refinement data structure

	// Copy the DT to the new structure
	double closest_dist = mesh.initFromVerticesAndTets(tin.vertices, tin.tet_node);
	if (verbose_mode) std::cout << "Distance of closest points relative to bb diagonal: " << closest_dist << "\n";

	now = chrono_clock::now();
	uint64_t time_DRinit = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	TetMesh *input_cdt = NULL;  
	uint64_t time_INcdt = 0;
	uint32_t INnct;
	if (lowBnd_onMeshDist){ 
		closest_dist = get_lower_bound_on_delRef_mesh_min_dist(plc, mesh, closest_dist, input_cdt, INnct);
		if (verbose_mode) std::cout << "Distance of closest elems relative to bb diagonal: " << closest_dist << "\n";
	}


	now = chrono_clock::now();
	uint64_t time_closedist = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (log_mode){ 
		logDouble(closest_dist);
		logDouble(time_closedist);
	}

	if (min_dist_exp != UINT32_MAX){
		double min_dist = std::pow(10.0, (double)min_dist_exp * (-1.0));
		if (closest_dist < min_dist){ 
			if(log_mode) finishLogging();
			std::cout<<"\nPROGRAM ABORTED: Closest points are too close ( < "<< min_dist <<")\n\n\n";
			exit(10);
		}
	}
	
	if (log_mode) advance_ProcessLogging("DelRef_init");

	// Add PLC faces - at this stage faces are just collections of input triangles
	if (verbose_mode) std::cout << "Adding PLC faces...\n";
	mesh.addPLCFaces(chamf_faces);

	now = chrono_clock::now();
	time_DRinit += std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (log_mode){ 
		logDouble(time_DRinit);
		advance_ProcessLogging("add_PLCfaces");
	}

	// Delaunay refinement algorithm
	// Split all missing segments while no more remain.
	// This does not modify the triangles in PLC faces, but updates their edges with
	// the subedges being computed.
	if (verbose_mode) std::cout << "Recovering segments...\n";
	
	mesh.recoverAllSegments();

	now = chrono_clock::now();
	uint64_t time_DRsr = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (log_mode){ 
		logDouble(time_DRsr);
		advance_ProcessLogging("rec_seg");
	}

	// Create a local 2D Delaunay triangulation for each PLC face
	if (verbose_mode) std::cout << "Delaunizing faces...\n";
	mesh.delaunizePLCFaces();
	if (log_mode) advance_ProcessLogging("Del_faces");

	// Face recovery
	if (verbose_mode) std::cout << "Recovering faces...\n";
	mesh.recoverAllFaces();

	now = chrono_clock::now();
	uint64_t time_DRfr = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (log_mode){ 
		logDouble(time_DRfr);
		advance_ProcessLogging("rec_faces");
	}

	
	// Tet optimization
	if (verbose_mode) {
		if(max_vrts==UINT32_MAX) std::cout << "Optimizing tets...\n";
		else std::cout << "Optimizing tets (max number of vertices = "<<max_vrts<<") ...\n";
	}
	mesh.optimizeTets(optim_ratio, true, true, max_vrts);

	now = chrono_clock::now();
	uint64_t time_DRopt = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if (log_mode){ 
		logDouble(time_DRopt);
		advance_ProcessLogging("tet_optim");
	}

	// Classify internal/external tets after optimization
	// In case the input edges has even number of incident input triangles,
	// i.e. the input surface defines an "interior",
	// we create a CDT (or use an existing one) of the input surface to decide wherever
	// an "otimized mesh" tetrahedron is internal/external wrt the input surface.
	if(input_encloses_vol){
		if (verbose_mode) std::cout<<"Input encloses a volume\n";
		if (input_cdt == NULL){ 
			input_cdt = createSteinerCDT(plc, false, false, INnct);

			now = chrono_clock::now();
			time_INcdt = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
			time_point = now;
		}
		markInternalTets(mesh, input_cdt);

		if (log_mode) advance_ProcessLogging("IntExt_class");
	}
	else{ for (Tetrahedron* t : mesh.tets()) t->is_internal = true; }

	now = chrono_clock::now();
	uint64_t time_DRie = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	time_point = now;

	if(log_mode) logDouble(time_DRie);

	// smoothing vertices TO REMOVE
	// smoothVertices(mesh, plc.numVertices()+8, optim_ratio);
	
	// chrono_clock::time_point now = chrono_clock::now();
	// uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	uint64_t ms = time_init + time_cham + time_closedist + time_DRinit + time_DRsr + time_DRfr + time_DRopt + time_DRie;
	// std::cout << "Elapsed time (ms): " << ms << "\n";
	// uint64_t bmem = getPeakRSS();
	// std::cout << "Peak memory RSS (byte): " << bmem << "\n";

	if (log_mode) {
		logUInteger((uint32_t)mesh.num_vertices());
		logUInteger((uint32_t)mesh.num_tetrahedra());
		logUInteger((uint32_t)mesh.count_DelTris(true)); // only constrained triangles
		logDouble(mesh.minEdgeLength());
		logFinalStats(mesh, ms, input_encloses_vol);
		advance_ProcessLogging("Optim");
		//finishLogging();
	}
	else {
		mesh.printReport(input_encloses_vol, "DelRef Mesh");
		std::cout << std::endl;
	}
	make_histogram(mesh, "DR");

	if(produce_DR_output){
		mesh.saveTET("DR_mesh.tet");
		mesh.saveOFFInterface("DR_plcfaces.off");
	}
	
	if (verbose_mode) std::cout << "\nClosing and computing CDT.\n";

	// Creates a .off file with all optimized vertices and
	// an optimized triangulation of the input plc, to be used as input for cdt
	std::vector<double> cdt_vrts;
	std::vector<uint32_t> cdt_tris;
	char empty_str[] = "";
	get_vrts_and_tris_for_cdt(mesh, to_close_ref_vrts, cdt_vrts, cdt_tris, empty_str, verbose_mode);

	// computes a quasi-optimized CDT
	if(true){
		inputPLC qo_plc; 
		qo_plc.initFromVectors(cdt_vrts.data(), cdt_vrts.size()/3, cdt_tris.data(), cdt_tris.size()/3, false);
		uint32_t nct=0;
		TetMesh* qo_cdt = createSteinerCDT(qo_plc, false, produce_outcdt_surf, nct);
		Tetrahedrization stat_mesh;
		stat_mesh.initFromVerticesAndTets(qo_cdt->vertices, qo_cdt->tet_node);
		for (Tetrahedron* t : stat_mesh.tets()) t->is_internal = true;

		if(log_mode){
			now = chrono_clock::now();
			uint64_t time_CDT = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
			time_point = now;
			logDouble(time_CDT);
			logUInteger( stat_mesh.num_vertices() );
			logUInteger( stat_mesh.num_tetrahedra() );
			logUInteger( nct );
			logFinalStats(stat_mesh, time_CDT, false);
			advance_ProcessLogging("final_CDT");
			// input CDT stats
			Tetrahedrization stat_INmesh;
			stat_INmesh.initFromVerticesAndTets(input_cdt->vertices, input_cdt->tet_node);
			for (Tetrahedron* t : stat_INmesh.tets()) t->is_internal = true;
			logDouble(time_INcdt);
			logUInteger( stat_INmesh.num_vertices() );
			logUInteger( stat_INmesh.num_tetrahedra() );
			logUInteger( INnct );
			logFinalStats(stat_INmesh, time_INcdt, false);
			advance_ProcessLogging("reg_inCDT_stats");
		}
		else { stat_mesh.printReport(false, "CDT of partially Delaunay refined mesh"); std::cout << std::endl; }
		// stat_mesh.saveOFFInterface("DRCDT_plcfaces.off");
		make_histogram(mesh, "CDT");
		make_histogram(mesh, "DR+CDT");


	}
	else{
		if(log_mode){
			logDouble(0.0);
			logUInteger( 0 );
			logUInteger( 0 );
			logUInteger( 0 );
		}
	}
	
	if(log_mode) finishLogging();
	std::cout << "Execution correctly COMPLETED.\n\n\n";

#endif

	return 0;
}
