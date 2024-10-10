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
		if (verbose) std::cout << "Chamfered PLC simplication COMPLETED: " << num_edges - cut_plc->edges.size() << " edges removed.\n\n";
	}

	std::vector<uint32_t> used_vertex(cut_plc->vertices.size(), UINT32_MAX);

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

	// cut_plc->saveFaces();

	return cut_plc->input_plc_defines_interior();
}

// createSteinerCDT
// 
// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!

TetMesh* createSteinerCDT(inputPLC& plc, bool verbose =false) {

	// Build a delaunay tetrahedrization of the vertices
	TetMesh* tin = new TetMesh;
	tin->init_vertices(plc.coordinates.data(), plc.numVertices());
	tin->tetrahedrize();

	if (verbose) printf("DT of the vertices built\n");

	// Build a structured PLC linked to the Delaunay tetrahedrization
	PLCx Steiner_plc(*tin, plc.triangle_vertices.data(), plc.numTriangles());
	Steiner_plc.segmentRecovery_HSi(!verbose);
	Steiner_plc.faceRecovery(!verbose);
	Steiner_plc.markInnerTets();

	return tin;
}

bool isTetInternal(Tetrahedron* t, TetMesh* cdt) {
	double ccc[3];
	const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
	getTetBarycenter(v, ccc);
	explicitPoint3D p(ccc[0], ccc[1], ccc[2]);
	cdt->vertices.push_back(&p);
	uint64_t tet = cdt->searchTetrahedron(0, cdt->numVertices() - 1);
	cdt->vertices.pop_back();
	return (cdt->mark_tetrahedra[tet>>2] == DT_IN);
}


int main(int argc, char* argv[])
{
	initFPU();

	char filename[2048];
#ifndef DEBUG
	if (argc < 2) {
		std::cout << "Mesher - Create a well-shaped tetrahedral mesh out of a triangulated OFF file.\n";
		std::cout << "USAGE: ./delmesher [-v][-t60] filename.off\n";
		std::cout << "OUTPUT:\n";
		std::cout << "plcfaces.off, mesh.tet\n";
		std::cout << "OPTIONS:\n";
		std::cout << "\t[-v] -> verbose mode\n";
		std::cout << "\t[-l] -> logging mode\n";
		std::cout << "\t[-t max time in minutes] -> time out mode";
		std::cout << "\t[-m max allocatable memory in Mb] -> memory out mode";
		std::cout << std::endl;
		return 0;
	}
	else strcpy(filename, argv[1]);
#else
	strcpy(filename, "..\\Input_file\\acute\\cup_fixed_fixed.off");
	//strcpy(filename, "..\\Input_file\\subdcube.off");
	//char filename[2048] = "..\\Input_file\\twocubes.off";

#endif

	std::string options = "";

	uint64_t time_out = 0;
	uint64_t mem_out = 0;
	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			if (argv[i][1] == 't') { time_out = atoi(argv[++i]); continue; }
			else if (argv[i][1] == 'm') { mem_out = atoi(argv[++i]); continue; }
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i]) + 1);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, options.find('v') != std::string::npos);

#ifndef DISP_PROCESS
	std::cout<<"\n\ninput_file: "<<filename<<"\n";
#endif
	
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

	tin.init_vertices(chamf_vertices);

	std::cout << "Delaunizing vertices...\n";
	tin.addBoundingBoxVertices();
	tin.tetrahedrize();

	Tetrahedrization mesh;
	std::chrono::steady_clock::time_point time_point = std::chrono::steady_clock::now();

	if (time_out > 0) {
		std::cout<<"Time out set at "<<time_out<<" min\n";
		time_out *= 60000; // convert to millisecond
		mesh.set_optimization_time_out(time_point, time_out);
	}

	if (mem_out > 0) {
		std::cout<<"Max memory set at "<<mem_out<<" Mb\n";
		mem_out *= 1000000; // convert to byte
		mesh.set_optimization_mem_out(mem_out);
	}

	// Copy the DT to the new structure
	double closest_pts_dist = mesh.initFromVerticesAndTets(tin.vertices, tin.tet_node);
	std::cout << "Distance of closest points relative to bb diagonal: " << closest_pts_dist << "\n";

	//if (closest_pts_dist < 1.0e-8) ip_error("Closest points are too close. Optimization would produce too many tets!\nEXITING\n");

	// Add PLC faces - at this stage faces are just collections of input triangles
	std::cout << "Adding PLC faces...\n";
	mesh.addPLCFaces(chamf_faces);

	// Split all missing segments while no more remain.
	// This does not modify the triangles in PLC faces, but updates their edges with
	// the subedges being computed.
	std::cout << "Recovering segments...\n";
	mesh.recoverAllSegments();

	// Create a local 2D Delaunay triangulation for each PLC face
	std::cout << "Delaunizing faces...\n";
	mesh.delaunizePLCFaces();

	// Face recovery
	std::cout << "Recovering faces...\n";
	mesh.recoverAllFaces();

	// Tet optimization
	std::cout << "Optimizing tets...\n";
	mesh.optimizeTets(2.0, false, true);

	// Remove external tets after chamfering
	// In case the input edges has even number of incident input triangles,
	// i.e. the input surface define an "interior",
	// we create a CDT of the input surface to decide wherever an "otimized mesh" tetrahedron 
	// is internal/external wrt the input surface.
	if(input_encloses_vol){
		std::cout<<"Input enclose a volume\n";
		TetMesh *cdt = createSteinerCDT(plc);
		for (Tetrahedron* t : mesh.tets()) t->is_internal = isTetInternal(t, cdt);
	}
	else{ for (Tetrahedron* t : mesh.tets()) t->is_internal = true; }
	
	std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
	uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
	std::cout << "Elapsed time (ms): " << ms << "\n";
	uint64_t bmem = getPeakRSS();
	std::cout << "Peak memory RSS (byte): " << bmem << "\n";

	if (options.find('l') != std::string::npos) {
		// append a line to logfile.csv

		// LOGfile
		startLogging(filename);
		logInteger((uint32_t)mesh.num_vertices());
		logInteger((uint32_t)mesh.num_tetrahedra());
		logInteger((uint32_t)ms);
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
		finishLogging();
	}
	else {
		mesh.printReport(input_encloses_vol);
		std::cout << std::endl;

		mesh.saveOFFInterface("plcfaces.off");
		//mesh.saveTET("mesh.tet");
	}

#endif
	//mesh.saveTET("mesh.tet");

	return 0;
}
