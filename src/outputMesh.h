#ifndef _OUTMESH_
#define _OUTMESH_

#include <iostream>
#include <fstream>

class output_mesh {

    const double* coords; // pointer to linearized vector of 3D coordinates
    uint32_t np; // number of points (3 coordinates for each point)
    const uint32_t* inds; // pointer to linearized vector of indices of vertices
    uint32_t ne; // number of elements: triangles, polygons or tetrahedra
    const char* filename;

    double xp, yp, zp;

public:
    output_mesh(const double* _coords, uint32_t _np,
                const uint32_t* _inds, uint32_t _ne,
                const char* _filename) : 
        coords(_coords), np(_np), inds(_inds), ne(_ne), filename(_filename) {};

    void get_point_coords(uint32_t i, double& x, double& y, double& z) {
        const double* ci = coords + i*3; 
        x = *ci;  y = *(ci + 1);  z = *(ci + 2);
    }

    void get_tri_vrts(uint32_t tri_i, uint32_t& i0, uint32_t& i1, uint32_t& i2) {
        const uint32_t* tv = inds + tri_i*3; 
        i0 = *tv;  i1 = *(tv + 1);  i2 = *(tv + 2);
    }

    void save_triangle_mesh() {
        
        FILE* file = fopen(filename, "w");
        
        fprintf(file, "OFF\n%zu %zu 0\n", np, ne);
        
        for(uint32_t i = 0; i < np; i++) {
            get_point_coords(i, xp, yp, zp);
            fprintf(file, "%f %f %f\n", xp, yp, zp);
        }
        
        for(uint32_t tri_i = 0; tri_i < ne; tri_i++) {
            uint32_t i0, i1, i2; 
            get_tri_vrts(tri_i, i0, i1, i2);
            fprintf(file, "3 %d %d %d\n", i0, i1, i2);
        }

        fclose(file);
    }

    void save_polygon_mesh(const uint32_t* fnv) {

        FILE* file = fopen(filename, "w");
        
        fprintf(file, "OFF\n%zu %zu 0\n", np, ne);
        
        for(uint32_t i = 0; i < np; i++) {
            get_point_coords(i, xp, yp, zp);
            fprintf(file, "%f %f %f\n", xp, yp, zp);
        }
        
        uint32_t pos = 0;
        for(uint32_t i = 0; i < ne; i++) {
            uint32_t n_pol_vrts = *(fnv + i);
            fprintf(file, "%zu ", n_pol_vrts);
            for(size_t j = 0; j < n_pol_vrts; j++) 
                fprintf(file, "%u ", *(inds + (pos++)));
            fprintf(file, "\n");
        }
        fclose(file);
    }
    
    void save_tetrahedral_mesh() {
        FILE* fp = fopen(filename, "w");
        fprintf(fp, "OFF\n%zu %zu 0\n", np, ne);
        // MISSING
    }

    // Minimal saves use only necessary vertices, i.e. those that compare
    // at least one as an element vretex.

    void save_minimal_triangle_mesh() {

        FILE* file = fopen(filename, "w");

        std::vector<uint32_t> used_vrt(np, UINT32_MAX);
	    for (size_t i = 0; i < ne; i++){ 
            const uint32_t* tri = inds + i*3;
            used_vrt[ *tri ] = 1;
            used_vrt[ *(tri +1) ] = 1;
            used_vrt[ *(tri +2) ] = 1;
        }

	    uint32_t idx = 0;
	    for (size_t i=0; i < used_vrt.size(); i++) 
            if (used_vrt[i] != UINT32_MAX) used_vrt[i] = idx++; // new indexing

        fprintf(file, "OFF\n%u %u 0\n", idx, ne);

        for(uint32_t i = 0; i < np; i++) if(used_vrt[i] != UINT32_MAX){
            get_point_coords(i, xp, yp, zp);
            fprintf(file, "%f %f %f\n", xp, yp, zp);
        }

        for(size_t tri_i = 0; tri_i < ne; tri_i++) {
            uint32_t i0, i1, i2; 
            get_tri_vrts(tri_i, i0, i1, i2);
            fprintf(file, "3 %d %d %d\n", 
                        used_vrt[i0], used_vrt[i1], used_vrt[i2]);
        }

        fclose(file);
    }

    void save_minimal_polygon_mesh(const uint32_t* fnv) {

        FILE* file = fopen(filename, "w");

        std::vector<uint32_t> used_vrt(np, UINT32_MAX);
        uint32_t pos = 0;
	    for (size_t i = 0; i < ne; i++){
            uint32_t n_pol_vrts = *(fnv + i);
            for(size_t j = 0; j < n_pol_vrts; j++) 
                used_vrt[  *(inds + (pos++)) ] = 1;
        }

	    uint32_t idx = 0;
	    for (size_t i=0; i < used_vrt.size(); i++) 
            if (used_vrt[i] != UINT32_MAX) used_vrt[i] = idx++; // new indexing

        fprintf(file, "OFF\n%u %u 0\n", idx, ne);

        for(uint32_t i = 0; i < np; i++) if(used_vrt[i] != UINT32_MAX){
            get_point_coords(i, xp, yp, zp);
            fprintf(file, "%f %f %f\n", xp, yp, zp);
        }

        pos = 0;
        for(uint32_t poly_i = 0; poly_i < ne; poly_i++) {
            uint32_t n_poly_vrts = *(fnv + poly_i);
            fprintf(file, "%zu ", n_poly_vrts);
            for(size_t j = 0; j < n_poly_vrts; j++) 
                fprintf(file, "%u ", used_vrt[ *(inds + (pos++)) ] );
            fprintf(file, "\n");
        }

        fclose(file);
    }
};

#endif