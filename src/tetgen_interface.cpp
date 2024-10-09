#include "tetgen.h"
#include "predicates.cxx"
#include "tetgen.cxx"

void Tetrahedrization::initWithTetgen(int num_vertices, double vertices[], int num_triangles, uint32_t triangles[], bool quality, bool no_erosion) {
	/* Copy triangulation vertices to tetgen data structures */
	int i;
	tetgenio in, out;
	tetgenio::polygon* pol;

	in.numberofpoints = num_vertices;
	in.pointlist = vertices;
	in.numberoffacets = num_triangles;
	in.facetlist = new tetgenio::facet[in.numberoffacets];

	for (i = 0; i < num_triangles; i++)
	{
		in.init(&(in.facetlist[i]));
		in.facetlist[i].numberofpolygons = 1;
		in.facetlist[i].polygonlist = new tetgenio::polygon[1];
		pol = &(in.facetlist[i].polygonlist[0]);
		in.init(pol);
		pol->numberofvertices = 3;
		pol->vertexlist = new int[3];
		pol->vertexlist[0] = triangles[i * 3];
		pol->vertexlist[1] = triangles[i * 3 + 1];
		pol->vertexlist[2] = triangles[i * 3 + 2];
	}

	/* Run TetGen to create Delaunay tetrahedrization of vertices */
	tetgenbehavior b;
	char cmdline[24];
	strcpy(cmdline, "zpv");
	if (quality) strcat(cmdline, "q");
	if (no_erosion) strcat(cmdline, "c");
	b.parse_commandline(cmdline);
	tetrahedralize(&b, &in, &out, NULL, NULL);

	 
	/* Create the tet mesh starting from tetgen data structure */
	V.reserve(out.numberofpoints);
	for (i = 0; i < out.numberofpoints; i++) {
		V.push_back(new TetVertex(out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]));
		initFullVE(V.back());
	}

	for (uint32_t t = 0; t < (uint32_t)out.numberoftetrahedra; t++)
		createTet(V[out.tetrahedronlist[t * 4]], V[out.tetrahedronlist[t * 4 + 1]], V[out.tetrahedronlist[t * 4 + 2]], V[out.tetrahedronlist[t * 4 + 3]]);

	for (auto* v : V) deleteFullVE(v);
	for (auto* e : E) deleteFullEF(e);

	in.pointlist = NULL; // To avoid TetGen to deallocate 'vertices' when deleting 'in'
}

