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
	std::vector<std::vector<uint32_t>>& _faces,
	bool simplify_cham_plc, bool verbose) {
	
	PLCc* cut_plc = new PLCc(_plc, _epsilon, true, verbose);

	if (simplify_cham_plc) {
		size_t num_edges = cut_plc->edges.size();
		cut_plc->chamfered_plc_simplification();
		// assert( cut_plc->checkup() );
		// assert( cut_plc->check_acuteness() );
		if (verbose) std::cout << "Chamfered PLC simplication COMPLETED: " << num_edges - cut_plc->edges.size() << " edges removed.\n";
	}

	std::vector<uint32_t> used_vertex(cut_plc->vertices.size(), UINT32_MAX);
	//std::vector<uint32_t> used_vertex(cut_plc->vertices.size(), 1);

	std::vector<uint32_t> out_tri;
	cut_plc->get_triangles(out_tri);
	if (verbose) std::cout << "Chamfered PLC triangulation COMPLETED\n\n";

	for (size_t i = 0; i < out_tri.size(); i++)	used_vertex[out_tri[i]] = 1;

	uint32_t idx = 0;
	for (size_t i = 0; i < used_vertex.size(); i++) if (used_vertex[i] != UINT32_MAX) {
		used_vertex[i] = idx++;
		_vertices.push_back(cut_plc->vertices[i]);
	}

	_faces.resize(out_tri.size() / 3);
	for (size_t i = 0; i < out_tri.size() / 3; i++) {
		_faces[i].push_back(used_vertex[out_tri[i * 3]]);
		_faces[i].push_back(used_vertex[out_tri[i * 3 + 1]]);
		_faces[i].push_back(used_vertex[out_tri[i * 3 + 2]]);
	}

	// cut_plc->saveFaces(); // save polygonal faces
	// cut_plc->saveTriFaces(); // triangulates and save faces

	return cut_plc->input_plc_defines_interior();
}

// createSteinerCDT
// 
// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!

TetMesh* createSteinerCDT(inputPLC& plc, double& min_PLC_dist, bool verbose =false) {

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

	return true;
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
	for (Tetrahedron* t : mesh.tets()) t->unmark<0>();
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
		std::cout << "USAGE: ./delmesher [-v][-l][-o][-t 60][-m 32000][-d 8] filename.off\n";
		std::cout << "OPTIONS:\n";
		std::cout << "\t[-v] -> verbose mode\n";
		std::cout << "\t[-l] -> logging mode\n";
		std::cout << "\t[-d exp (pos. integer)] -> avoid Delaunay Refinement if the min pt.s dist after chamfering is < 10^-exp\n";
		std::cout << "\t[-a] -> computes the lower bound for minimum distance between mesh elements (point - segment - triangle), otherwise only point-point distances are computed.\n";
		std::cout << "\t[-o] -> produces outputs (see below)\n";
		std::cout << "OUTPUT:\n";
		std::cout << "\t when [-o] is activated produces a volumetric (mesh.tet) and a surface (plcfaces.off) meshes.\n";
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
	bool produce_output = (options.find('o') != std::string::npos);
	bool lowBnd_onMeshDist = (options.find('a') != std::string::npos);

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
	std::vector<std::vector<uint32_t>> chamf_faces;
	double epsilon = plc.bbDiag() / 1000.0;	// Use bounding-box-diagonal/1000 as chamfering distance

	bool simplify_chamferd_plc = true; // try to remove uncessary edges while keeping non-acute angles
	bool input_encloses_vol = chamferPLC(plc, epsilon, chamf_vertices, chamf_faces, simplify_chamferd_plc, options.find('v') != std::string::npos);

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
		cdt = createSteinerCDT(plc, min_PLC_dist); // needed to compute lower bound on generic mesh elements distances
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
	mesh.optimizeTets(2.0, false, true);

	if (log_mode) advance_ProcessLogging("tet_optim");

	// Remove external tets after chamfering
	// In case the input edges has even number of incident input triangles,
	// i.e. the input surface defines an "interior",
	// we create a CDT of the input surface to decide wherever an "otimized mesh" tetrahedron 
	// is internal/external wrt the input surface.
	if(input_encloses_vol){
		std::cout<<"Input encloses a volume\n";
		if(cdt==NULL){ 
			double zero=0.0; 
			cdt = createSteinerCDT(plc, zero);
		}
		markInternalTets(mesh, cdt);

		if (log_mode) advance_ProcessLogging("IntExt_class");
	}
	else{ for (Tetrahedron* t : mesh.tets()) t->is_internal = true; }
	
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

	if(produce_output){
		mesh.saveTET("mesh.tet");
		mesh.saveOFFInterface("plcfaces.off");
	}
	
	std::cout << "Execution correctly COMPLETED.\n\n\n";

#endif

	return 0;
}
