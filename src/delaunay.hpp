using namespace std;

void TetMesh::init_vertices(const double* coords, uint32_t num_v) {
    vertices.reserve(num_v);
    for (uint32_t i = 0; i < num_v; i++)
        vertices.push_back(new explicitPoint(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
    inc_tet.resize(num_v, UINT64_MAX);
    marked_vertex.resize(num_v, 0);
}

void TetMesh::init_vertices(std::vector<genericPoint*>& pts) {
    vertices.assign(pts.begin(), pts.end());
    inc_tet.resize(vertices.size(), UINT64_MAX);
    marked_vertex.resize(vertices.size(), 0);
}

void TetMesh::addBoundingBoxVertices(const std::vector<double>& bbox_coords) {
    assert(bbox_coords.size() == 24);
    vertices.insert( vertices.end(), {
        new explicitPoint3D(bbox_coords[0], bbox_coords[1], bbox_coords[2]),
        new explicitPoint3D(bbox_coords[3], bbox_coords[4], bbox_coords[5]),
        new explicitPoint3D(bbox_coords[6], bbox_coords[7], bbox_coords[8]),
        new explicitPoint3D(bbox_coords[9], bbox_coords[10], bbox_coords[11]),
        new explicitPoint3D(bbox_coords[12], bbox_coords[13], bbox_coords[14]),
        new explicitPoint3D(bbox_coords[15], bbox_coords[16], bbox_coords[17]),
        new explicitPoint3D(bbox_coords[18], bbox_coords[19], bbox_coords[20]),  
       new explicitPoint3D(bbox_coords[21], bbox_coords[22], bbox_coords[23])});
    inc_tet.resize(vertices.size(), UINT64_MAX);
    marked_vertex.resize(vertices.size(), 0);
}

void TetMesh::init(uint32_t& unswap_k, uint32_t& unswap_l){
  const uint32_t n = numVertices();

  // Find non-coplanar vertices (we assume that no coincident vertices exist)
  int ori=0;
  uint32_t i=0, j=1, k=2, l=3;

  for (; ori == 0 && k < n - 1; k++)
      for (l = k + 1; ori == 0 && l < n; l++)
          ori = vOrient3D(i, j, k, l);

  l--; k--;

  if(ori==0){
    // ip_error("TetMesh::init() - Input vertices do not define a volume.\n");
    std::cout<<"[delaunay.hpp - init()] Input vertices do not define a volume.\n"; exit(1);
  }

  unswap_k = k;
  unswap_l = l;
  std::swap(vertices[k], vertices[2]); k=2;
  std::swap(vertices[l], vertices[3]); l=3;

  if(ori<0) std::swap(i, j); // Tets must have positive volume

  const uint32_t base_tet[] = { l, k, j, i, l, j, k, INFINITE_VERTEX, l, k, i, INFINITE_VERTEX, l, i, j, INFINITE_VERTEX, k, j, i, INFINITE_VERTEX };
  const uint64_t base_neigh[] = { 19, 15, 11, 7, 18, 10, 13, 3, 17, 14, 5, 2, 16, 6, 9, 1, 12, 8, 4, 0 };

  resizeTets(5);
  std::memcpy(getTetNodes(0), base_tet, 20 * sizeof(uint32_t));
  std::memcpy(getTetNeighs(0), base_neigh, 20 * sizeof(uint64_t));

  // set the vertex-(one_of_the)incident-tetrahedron relation
  inc_tet[i] = inc_tet[j] = inc_tet[k] = inc_tet[l] = 0;
}


void TetMesh::tetrahedrize() {
    uint32_t uk, ul;
    init(uk, ul); // First tet is made of vertices 0, 1, uk, ul

    // Need to unswap immediately to keep correct indexing and
    // ensure symbolic perturbation is coherent
    if (ul != 3) {
        std::swap(vertices[ul], vertices[3]);
        std::swap(inc_tet[ul], inc_tet[3]);
        for (uint32_t& tn : tet_node) if (tn == 3) tn = ul; else if (tn == ul) tn = 3;
    }

    if (uk != 2) {
        std::swap(vertices[uk], vertices[2]);
        std::swap(inc_tet[uk], inc_tet[2]);
        for (uint32_t& tn : tet_node) if (tn == 2) tn = uk; else if (tn == uk) tn = 2;
    }

    uint64_t ct = 0;
    for (uint32_t i = 2; i < numVertices(); i++) if (i != uk && i != ul) insertExistingVertex(i, ct);

    removeDelTets();
}


bool TetMesh::saveTET(const char* filename, bool inner_only) const
{
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveTET: Can't open file for writing.\n";
        return false;
    }

    f << numVertices() << " vertices\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    if (inner_only) {
        f << ngnt << " tets\n";
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    else {
        f << ngnt << " inner tets\n";
        f << countNonGhostTets()-ngnt << " outer tets\n";
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    
    f.close();

    return true;
}


bool TetMesh::saveMEDIT(const char* filename, bool inner_only) const
{
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveMEDIT: Can't open file for writing.\n";
        return false;
    }

    f << "MeshVersionFormatted 2\nDimension\n3\n";

    f << "Vertices\n" << numVertices() << "\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    if (inner_only) {
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << " 1\n";
        f << "Tetrahedra\n" << ngnt << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << tet_node[i * 4]+1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 1\n";
    }
    else {
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << " 1\n";
        f << "Tetrahedra\n" << countNonGhostTets() << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 1\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 2\n";
    }

    f.close();

    return true;
}


bool TetMesh::saveBinaryTET(const char* filename, bool inner_only) const
{
    ofstream f(filename, ios::binary);

    if (!f) {
        std::cerr << "\nTetMesh::saveBinaryTET: Can't open file for writing.\n";
        return false;
    }

    uint32_t num_v = numVertices(), num_t = 0;

    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) num_t++;

    f << num_v << " vertices\n";

    if (inner_only) {
        f << num_t << " tets\n";
    }
    else {
        f << num_t << " inner tets\n";
        f << countNonGhostTets() - num_t << " outer tets\n";
    }

    double c[3];
    for (uint32_t i = 0; i < numVertices(); i++) {
        vertices[i]->getApproxXYZCoordinates(c[0], c[1], c[2], true);
        f.write((const char*)(&c), sizeof(double) * 3);
    }

    const uint32_t* tnd = tet_node.data();

    if (inner_only) {
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
    }
    else {
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
    }

    f.close();

    return true;
}

bool TetMesh::saveBoundaryToOFF(const char* filename) const {
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveBoundaryToOFF: Can't open file for writing.\n";
        return false;
    }

    f << "OFF\n" << numVertices() << " ";

    size_t num_tris = 0;
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) num_tris++;

    f << num_tris << " 0\n";

    for (uint32_t i = 0; i < numVertices(); i++)
        f << *vertices[i] << "\n";

    uint32_t fv[3];
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) {
            getFaceVertices(i, fv);
            f << "3 " << fv[0] << " " << fv[1] << " " << fv[2] << "\n";
        }
    f.close();

    return true;
}

bool TetMesh::saveConstrTrisToOFF(const char* filename, const std::vector<bool>& constrTris) const {
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveBoundaryToOFF: Can't open file for writing.\n";
        return false;
    }

    f << "OFF\n" << numVertices() << " ";

    size_t num_tris = 0;
    for (uint64_t i = 0; i < tet_node.size(); i++) 
        if (i > tet_neigh[i] && constrTris[i] ) num_tris++;

    f << num_tris << " 0\n";

    for (uint32_t i = 0; i < numVertices(); i++)
        f << *vertices[i] << "\n";

    uint32_t fv[3];
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && constrTris[i]) {
            getFaceVertices(i, fv);
            f << "3 " << fv[0] << " " << fv[1] << " " << fv[2] << "\n";
        }
    f.close();

    return true;
}

uint32_t TetMesh::countConstrTris(const std::vector<bool>& constrTris) const {

    size_t num_tris = 0;
    for (uint64_t i = 0; i < tet_node.size(); i++) 
        if (i > tet_neigh[i] && constrTris[i] ) num_tris++;

    return (uint32_t)num_tris;
}

bool TetMesh::saveRationalTET(const char* filename, bool inner_only)
{
#ifdef USE_INDIRECT_PREDS
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveRationalTET: Can't open file for writing.\n";
        return false;
    }

    f << numVertices() << " vertices\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    if (inner_only) {
        f << ngnt << " tets\n";
        for (uint32_t i = 0; i < numVertices(); i++) {
            bigrational c[3];
            vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
            f << c[0] << " " << c[1] << " " << c[2] << "\n";
        }
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    else {
        f << ngnt << " inner tets\n";
        f << countNonGhostTets() - ngnt << " outer tets\n";
        for (uint32_t i = 0; i < numVertices(); i++) {
            bigrational c[3];
            vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
            f << c[0] << " " << c[1] << " " << c[2] << "\n";
        }
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }

    f.close();
#endif

    return true;
}

// Swap t and l while assuming that l is a valid tet (not to be deleted)
inline void TetMesh::moveDeletedToTail(uint64_t t, uint64_t l) {
    const uint64_t t2 = t >> 2, l2 = l >> 2;
    uint32_t* tnt = tet_node.data() + t, *lnt = tet_node.data() + l;
    uint64_t* tnn = tet_neigh.data() + t, *lnn = tet_neigh.data() + l;
    const uint32_t* te = tnt + 4;

    do {
        *tnt++ = *lnt;

        uint64_t neigh = *lnn++;
        *tnn++ = neigh;
        tet_neigh[neigh] = t++;

        if (*lnt != INFINITE_VERTEX && inc_tet[*lnt] == l2)
            inc_tet[*lnt] = t2;
        lnt++;
    } while (tnt != te);

    mark_tetrahedra[t2] = mark_tetrahedra[l2];
}

void TetMesh::removeDelTets() {
  uint64_t tn = numTets();
  uint64_t* dp = Del_deleted.data();
  uint64_t* de = dp + Del_deleted.size();

  while (dp != de) {
      uint64_t last = (--tn) * 4;

      if (isToDelete(last))
      {
          uint64_t* dj = dp;
          while (++dj != de && *dj != last);
          *dj = *dp;
      }
      else moveDeletedToTail(*dp, last);

      dp++;
  }

  resizeTets(tn);
  Del_deleted.clear();
}

bool TetMesh::tetHasVertex(uint64_t t, uint32_t v) const {
    t <<= 2;
    return tet_node[t] == v || tet_node[t + 1] == v || tet_node[t + 2] == v || tet_node[t + 3] == v;
}

void TetMesh::oppositeTetEdge(const uint64_t tet, const uint32_t v[2], uint32_t ov[2]) const {
    int i = 0, j = 0;
    while (i < 4) {
        const uint32_t w = tet_node[tet + i];
        if (w != v[0] && w != v[1]) ov[j++] = w;
        i++;
    }
    assert(j == 2);
}

uint64_t TetMesh::getCornerFromOppositeTet(uint64_t t, uint64_t n) const {
    t <<= 2;
    for (int i = 0; i < 4; i++)
        if ((tet_neigh[t + i] >> 2) == n)
            return tet_neigh[t + i];
    assert(0);
    return UINT64_MAX;
}

void TetMesh::getFaceVertices(uint64_t t, uint32_t v[3]) const {
    uint64_t tv = t & 3;
    const uint32_t* Node = tet_node.data() + (t - tv);
    v[0] = Node[(++tv) & 3];
    v[1] = Node[(++tv) & 3];
    v[2] = Node[(++tv) & 3];
}

bool TetMesh::getTetsFromFaceVertices(uint32_t v1, uint32_t v2, uint32_t v3, uint64_t* nt) const {
    static std::vector<uint64_t> vt; // Static to avoid reallocation at each call
    VTfull(v1, vt);
    int i = 0;
    for (uint64_t t : vt) if (tetHasVertex(t, v2) && tetHasVertex(t, v3)) nt[i++] = t;
    vt.clear();
    return (i == 2);
}

uint64_t TetMesh::tetOppositeCorner(uint64_t t, uint32_t v1, uint32_t v2, uint32_t v3) const {
    const uint64_t tb = t << 2;
    const uint32_t* n = tet_node.data() + tb;
    for (int i = 0; i < 3; i++)
        if (n[i] != v1 && n[i] != v2 && n[i] != v3)
            return tet_neigh[tb + i];
    assert(n[3] != v1 && n[3] != v2 && n[3] != v3);
    return tet_neigh[tb + 3];
}

void TetMesh::resizeTets(uint64_t new_size) {
    mark_tetrahedra.resize(new_size, 0);
    new_size <<= 2;
    tet_node.resize(new_size);
    tet_neigh.resize(new_size);
}

void TetMesh::reserveTets(uint64_t new_capacity) {
    mark_tetrahedra.reserve(new_capacity);
    new_capacity <<= 2;
    tet_node.reserve(new_capacity);
    tet_neigh.reserve(new_capacity);
}

uint64_t TetMesh::searchTetrahedron(uint64_t tet, const uint32_t v_id)
{
    if (tet_node[tet + 3] == INFINITE_VERTEX)
        tet = getIthNeighbor(getTetNeighs(tet), 3);

    uint64_t i, f0 = 4;
    uint64_t count = 0, maxcnt = numTets();

    do {
        const uint32_t* Node = getTetNodes(tet);
        if (Node[3] == INFINITE_VERTEX || count++ == maxcnt) return tet;

        const uint64_t* Neigh = getTetNeighs(tet);
        for (i = 0; i < 4; i++)
            if (i != f0 && vOrient3D(Node[tetON1(i)], Node[tetON2(i)], Node[tetON3(i)], v_id) < 0) {
                tet = getIthNeighbor(Neigh, i);
                f0 = Neigh[i] & 3;
                break;
            }
    } while (i != 4);

    return tet;
}

uint64_t TetMesh::searchTetrahedron(uint64_t tet, const pointType* pt)
{
    if (tet_node[tet + 3] == INFINITE_VERTEX)
        tet = getIthNeighbor(getTetNeighs(tet), 3);

    uint64_t i, f0 = 4;
    do {
        const uint32_t* Node = getTetNodes(tet);
        if (Node[3] == INFINITE_VERTEX) return tet;

        const uint64_t* Neigh = getTetNeighs(tet);
        for (i = 0; i < 4; i++)
            if (i != f0 && (-pointType::orient3D(*vertices[Node[tetON1(i)]],*vertices[Node[tetON2(i)]],*vertices[Node[tetON3(i)]],*pt)) < 0) {
                tet = getIthNeighbor(Neigh, i);
                f0 = Neigh[i] & 3;
                break;
            }
    } while (i != 4);

    return tet;
}


int TetMesh::symbolicPerturbation(uint32_t indices[5]) const {
    int swaps = 0;
    int n = 5;
    int count;
    do {
        count = 0;
        n--;
        for (int i = 0; i < n; i++) {
            if (indices[i] > indices[i + 1]) {
                std::swap(indices[i], indices[i + 1]);
                count++;
            }
        }
        swaps += count;
    } while (count);

    n = vOrient3D(indices[1], indices[2], indices[3], indices[4]);
    if (n) return (swaps % 2) ? (-n) : n;

    n = vOrient3D(indices[0], indices[2], indices[3], indices[4]);
    return (swaps % 2) ? (n) : (-n);
}

int TetMesh::vertexInTetSphere(const uint32_t Node[4], uint32_t v_id) const {
    int det = vInSphere(Node[0], Node[1], Node[2], Node[3], v_id);
    if (det) return det;
    uint32_t nn[5] = { Node[0],Node[1],Node[2],Node[3],v_id };
    det = symbolicPerturbation(nn);
    if (det == 0.0){ 
        //ip_error("Symbolic perturbation failed! Should not happen.\n");
        std::cout<<"[delaunay.hpp - vertexInSphere()] Symbolic perturbation failed! Should not happen\n"; exit(1);
    }
    return det;
}

int TetMesh::vertexInTetSphere(uint64_t tet, uint32_t v_id) const
{
  const uint32_t* Node = getTetNodes(tet);
  int det;

  if (Node[3] == INFINITE_VERTEX) {
      if ((det = vOrient3D(Node[0], Node[1], Node[2], v_id)) != 0) return det;
      const uint32_t nn[4] = {Node[0], Node[1], Node[2], tet_node[tet_neigh[tet + 3]]};
      return -vertexInTetSphere(nn, v_id);
  }
  else return vertexInTetSphere(Node, v_id);
}

#ifdef USE_MAROTS_METHOD
void TetMesh::deleteInSphereTets(uint64_t tet, const uint32_t v_id)
{
  pushAndMarkDeletedTets(tet);

  for(uint64_t t = Del_deleted.size() - 1; t < Del_deleted.size(); t++) {
    uint64_t tet = Del_deleted[t];
    uint64_t* Neigh = getTetNeighs(tet);
    uint32_t* Node = getTetNodes(tet);

    uint64_t neigh = getIthNeighbor(Neigh, 0);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 1);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 2);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 3);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0){
        if(Node[1]<Node[2])
          bnd_push(v_id, Node[0], Node[2], Node[1], Neigh[3]);
        else
          bnd_push(v_id, Node[1], Node[0], Node[2], Neigh[3]);
      }
      else pushAndMarkDeletedTets(neigh);
    }
  }
}


void TetMesh::tetrahedrizeHole(uint64_t* tet){
  uint64_t clength = Del_deleted.size(); // Num tets removed
  uint64_t blength = numDelTmp(); // Num tets to insert

  uint64_t tn = numTets();

  if(blength > clength){
    for (uint64_t i = clength; i<blength; i++, tn++)
        Del_deleted.push_back(tn<<2);

    clength = blength;
    resizeTets(tn);
  }

  uint64_t start = clength - blength;

  for (uint64_t i=0; i<blength; i++)
  {
    const uint64_t tet = Del_deleted[i + start];
    uint32_t* Node = getTetNodes(tet);

    Node[0] = Del_tmp[i].node[0];
    Node[1] = Del_tmp[i].node[1];
    Node[2] = Del_tmp[i].node[2];
    Node[3] = Del_tmp[i].node[3];

    uint64_t bnd = Del_tmp[i].bnd;
    tet_neigh[tet] = bnd;
    tet_neigh[bnd] = tet;
    Del_tmp[i].bnd = tet;

    mark_tetrahedra[tet >> 2] = 0;

    if(tet_node[tet+3]!=INFINITE_VERTEX)
      for(uint32_t j=0; j<4; j++)
          inc_tet[tet_node[tet + j]] = tet>>2;
  }

  uint64_t tlength = 0;
  const uint64_t middle = blength * 3 / 2;

  uint64_t* Tmp = delTmpVec();
  const unsigned index[4] = { 2,3,1,2 };

  for (uint64_t i = 0; i < blength; i++)
  {
      uint64_t tet = Del_deleted[start + i];
      const uint32_t* Node = getTetNodes(tet);

      for (uint64_t j = 0; j < 3; j++)
      {
          uint64_t key = ((uint64_t)Node[index[j]] << 32) + Node[index[j + 1]];
          tet++;

          uint64_t k;
          for (k = 0; k < tlength; k++) if (Tmp[k] == key) break;

          if (k == tlength) {
              Tmp[tlength] = (key >> 32) + (key << 32);
              Tmp[middle + tlength] = tet;
              tlength++;
          }
          else {
              uint64_t pairValue = Tmp[middle + k];
              tet_neigh[tet] = pairValue;
              tet_neigh[pairValue] = tet;
              tlength--;
              if (k < tlength) {
                  Tmp[k] = Tmp[tlength];
                  Tmp[middle + k] = Tmp[middle + tlength];
              }
          }
      }
  }

  flushDelTmp();
  *tet = Del_deleted[start];
  Del_deleted.resize(start);
}

void TetMesh::insertExistingVertex(const uint32_t vi, uint64_t& ct)
{
    ct = searchTetrahedron(ct, vi);
    deleteInSphereTets(ct, vi);
    tetrahedrizeHole(&ct);
    uint64_t lt = ct;
    if (tet_node[lt + 3] == INFINITE_VERTEX) lt = tet_neigh[lt + 3];
    inc_tet[vi] = lt >> 2;
}

#else
// Start from c and turn around v1-v2 as long as adjacencies are well defined.
// When an invalid adjacency is found, reinit it and exit.
void TetMesh::seekAndSetMutualAdjacency(int p_o0, int p_o1, int p_o2, const uint32_t* v, uint64_t c, uint64_t o, const uint32_t* tet_node_data, uint64_t* tet_neigh_data) {
    const uint32_t ov = v[p_o0], v1 = v[p_o1], v2 = v[p_o2];
    o += p_o0;

    c &= (~3);
    while (tet_node_data[c] != ov) c++;

    for (;;) {
        uint64_t t = c;
        if ((c = tet_neigh_data[c]) == UINT64_MAX) {
            tet_neigh_data[t] = o;
            tet_neigh_data[o] = t;
            return;
        }
        const uint32_t w = tet_node_data[c];
        c &= (~3);
        while (tet_node_data[c] == v1 || tet_node_data[c] == v2 || tet_node_data[c] == w) c++;
    }
}

// Rebuild internal adjacencies for the cavity tet opposite to c
void TetMesh::restoreLocalConnectivty(uint64_t c, const uint32_t* tet_node_data, uint64_t* tet_neigh_data) {
    const uint64_t o = tet_neigh_data[c];
    const uint32_t* v = tet_node_data + o;
    const uint64_t* n = tet_neigh_data + o;
    if (n[1] == UINT64_MAX) seekAndSetMutualAdjacency(1, 2, 3, v, c, o, tet_node_data, tet_neigh_data);
    if (n[2] == UINT64_MAX) seekAndSetMutualAdjacency(2, 1, 3, v, c, o, tet_node_data, tet_neigh_data);
    if (n[3] == UINT64_MAX) seekAndSetMutualAdjacency(3, 1, 2, v, c, o, tet_node_data, tet_neigh_data);
}

// Collect all tets whose circumsphere contains v_id and replace them
// with a star of new tets originating at v_id
void TetMesh::insertExistingVertex(const uint32_t v_id, uint64_t& tet)
{
    static std::vector<uint64_t> cavityCorners; // Static to avoid reallocation on each call
    static const int fi[4][3] = { {2, 1, 3} ,{0, 2, 3} ,{1, 0, 3} ,{0, 1, 2} };
    uint32_t* tet_node_data = tet_node.data();
    uint64_t* tet_neigh_data = tet_neigh.data();

    // Move by adjacencies to find the tet containing v_id
    if (tet_node_data[tet + 3] == INFINITE_VERTEX)
        tet = tet_neigh_data[tet + 3] & (~3);

    uint64_t i, f0 = 4;
    do {
        const uint32_t* Node = tet_node_data + tet;
        if (Node[3] == INFINITE_VERTEX) break;

        for (i = 0; i < 4; i++)
            if (i != f0 && vOrient3D(Node[tetON1(i)], Node[tetON2(i)], Node[tetON3(i)], v_id) < 0) {
                const uint64_t ni = tet_neigh_data[tet + i];
                tet = ni & (~3);
                f0 = ni & 3;
                break;
            }
    } while (i != 4);

    tet >>= 2;

    // Expand by adjacencies to collect all tets whose circumsphere contains v_id
    size_t first = Del_deleted.size();
    pushAndMarkDeletedTets(tet << 2);

    for (size_t i = first; i < Del_deleted.size(); i++) {
        const uint64_t* nb = tet_neigh_data + Del_deleted[i];
        const uint64_t* nl = nb + 4;

        for (; nb < nl; nb++)
        {
            const uint64_t n0 = *nb >> 2;
            uint32_t& mtn0 = mark_tetrahedra[n0];
            if (mtn0 == 0) {
                if (vertexInTetSphere(n0 << 2, v_id) < 0) {
                    mtn0 = 2;
                    cavityCorners.push_back(*nb);
                }
                else {
                    pushAndMarkDeletedTets(n0 << 2);
                }
            }
            else if (mtn0 == 2) cavityCorners.push_back(*nb);
        }
    }


    // Resize the mesh to host the new tets
    uint64_t ntb, newpos = tet_node.size();
    if (cavityCorners.size() > Del_deleted.size()) {
        resizeTets(numTets() + (cavityCorners.size() - Del_deleted.size()));
        tet_node_data = tet_node.data();
        tet_neigh_data = tet_neigh.data();
    }

    // Create the new tets
    for (const uint64_t c : cavityCorners) {
        mark_tetrahedra[c >> 2] = 0;
        if (Del_deleted.empty()) {
            ntb = newpos;
            newpos += 4;
        }
        else {
            ntb = Del_deleted.back();
            Del_deleted.pop_back();
        }
        const uint64_t cb = c & 3;
        const uint32_t* cr = tet_node_data + (c - cb);
        uint32_t* cn = tet_node_data + ntb;
        *cn++ = v_id;
        *cn++ = cr[fi[cb][0]];
        *cn++ = cr[fi[cb][1]];
        *cn++ = cr[fi[cb][2]];

        tet_neigh_data[ntb] = c; tet_neigh_data[c] = ntb;
        tet_neigh_data[ntb + 1] = tet_neigh_data[ntb + 2] = tet_neigh_data[ntb + 3] = UINT64_MAX;

        ntb >>= 2;
        if ((*(--cn)) != INFINITE_VERTEX) {
            inc_tet[*cn] = ntb;
            inc_tet[*(--cn)] = ntb;
            inc_tet[*(--cn)] = ntb;
            inc_tet[v_id] = ntb;
        }
        mark_tetrahedra[ntb] = 0;
    }

    // Restore the connectivity within the cavity
    for (uint64_t c : cavityCorners) restoreLocalConnectivty(c, tet_node_data, tet_neigh_data);

    tet = tet_neigh_data[cavityCorners.back()];

    cavityCorners.clear();
}
#endif
void TetMesh::VT(uint32_t v, std::vector<uint64_t>& vt) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v];

    vt_queue.push_back(tetCornerAtVertex(t << 2, v));
    mark_Tet_31(t);

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb+j)&3];
            const uint64_t tbb = tb >> 2;
            if (tet_node[tb] != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) { 
                vt_queue.push_back(tetCornerAtVertex(tb & (~3), v)); 
                mark_Tet_31(tbb); 
            }
        }
    }

    for (uint64_t t : vt_queue) {
        t >>= 2;
        unmark_Tet_31(t);
        vt.push_back(t);
    }
    vt_queue.clear();
}

void TetMesh::VV(uint32_t v, std::vector<uint32_t>& vv) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v];
    const uint64_t tb = t << 2;

    const uint64_t s = tetCornerAtVertex(tb, v);
    vt_queue.push_back(s);
    mark_Tet_31(t);

    const uint32_t* tn = tet_node.data() + tb;
    const uint64_t sb = s & 3;
    for (int j = 1; j < 4; j++) {
        const uint32_t w = tn[(sb + j) & 3];
        marked_vertex[w] |= 128;
        vv.push_back(w);
    }

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb + j) & 3];
            const uint64_t tbb = tb >> 2;
            const uint32_t w = tet_node[tb];
            if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                vt_queue.push_back(tetCornerAtVertex(tb & (~3), v));
                mark_Tet_31(tbb);
                if (!(marked_vertex[w] & 128)) {
                    marked_vertex[w] |= 128;
                    vv.push_back(w);
                }
            }
        }
    }

    for (uint64_t t : vt_queue) unmark_Tet_31(t>>2);
    vt_queue.clear();
    for (uint32_t w : vv) marked_vertex[w] &= 127;
}

void TetMesh::ET(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    VT(v1, et);
    for (size_t i = 0; i < et.size();)
        if (!tetHasVertex(et[i], v2)) {
            std::swap(et[i], et[et.size() - 1]);
            et.pop_back();
        }
        else i++;
}

void TetMesh::ETfull(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    VTfull(v1, et);
    for (size_t i = 0; i < et.size();)
        if (!tetHasVertex(et[i], v2)) {
            std::swap(et[i], et[et.size() - 1]);
            et.pop_back();
        }
        else i++;
}

void TetMesh::ETcorners(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    uint64_t t;
    VTfull(v1, et);
    for (uint64_t s : et) if (tetHasVertex(s, v2)) { t = (s<<2); break; }
    
    while (tet_node[t] == v1 || tet_node[t] == v2) t++;

    et.clear();

    uint64_t c0 = t;
    do {
        et.push_back(t); // Add tet
        uint64_t oc = tet_neigh[t] & (~3); // Get next base
        uint32_t cv = tet_node[t];
        t &= (~3);
        while (tet_node[t] == v1 || tet_node[t] == v2 || tet_node[t] == cv) t++;
        t = tetCornerAtVertex(oc, tet_node[t]); // Get corresp corner at opposite tet
    } while (t != c0);
}

void TetMesh::VTfull(uint32_t v, std::vector<uint64_t>& vt) const {

    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t s, t = inc_tet[v];
    vt_queue.push_back(t);
    mark_Tet_31(t);

    while (!vt_queue.empty()) {
        t = vt_queue.back();
        vt_queue.pop_back();
        vt.push_back(t);
        t <<= 2;
        s = tet_neigh[t] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 1] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 2] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 3] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
    }

    for (uint64_t t : vt) unmark_Tet_31(t);
}


bool TetMesh::hasEdge(uint32_t v1, uint32_t v2) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v1];
    const uint64_t tb = t << 2;
    if (tet_node[tb] == v2 || tet_node[tb + 1] == v2 || tet_node[tb + 2] == v2 || tet_node[tb + 3] == v2) return true;

    vt_queue.push_back(tetCornerAtVertex(tb, v1));
    mark_Tet_31(t);

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb + j) & 3];
            const uint64_t tbb = tb >> 2;
            const uint32_t w = tet_node[tb];
            if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                vt_queue.push_back(tetCornerAtVertex(tbb << 2, v1));
                mark_Tet_31(tbb);
                if (w == v2) {
                    for (uint64_t t : vt_queue) unmark_Tet_31(t >> 2);
                    vt_queue.clear();
                    return true;
                }
            }
        }
    }

    for (uint64_t t : vt_queue) unmark_Tet_31(t >> 2);
    vt_queue.clear();
    return false;
}


void TetMesh::swapTets(const uint64_t t1, const uint64_t t2) 
{
    if (t1 == t2) return;

    const uint64_t t1_id = t1<<2;
    const uint64_t t2_id = t2<<2;

    // update VT base relation
    for (int i = 0; i < 3; i++) if (inc_tet[tet_node[t1_id + i]] == t1) inc_tet[tet_node[t1_id + i]] = t2;
    if (tet_node[t1_id + 3] != INFINITE_VERTEX && inc_tet[tet_node[t1_id + 3]] == t1) inc_tet[tet_node[t1_id + 3]] = t2;

    for (int i = 0; i < 3; i++) if (inc_tet[tet_node[t2_id + i]] == t2) inc_tet[tet_node[t2_id + i]] = t1;
    if (tet_node[t2_id + 3] != INFINITE_VERTEX && inc_tet[tet_node[t2_id + 3]] == t2) inc_tet[tet_node[t2_id + 3]] = t1;

    // Update nodes and marks
    for (int i = 0; i < 4; i++) std::swap(tet_node[t1_id + i], tet_node[t2_id + i]);
    std::swap(mark_tetrahedra[t1], mark_tetrahedra[t2]);

    // update neigh-neigh relations
    const uint64_t ng1[] = { tet_neigh[t1_id + 0], tet_neigh[t1_id + 1], tet_neigh[t1_id + 2], tet_neigh[t1_id + 3] };
    const uint64_t ng2[] = { tet_neigh[t2_id + 0], tet_neigh[t2_id + 1], tet_neigh[t2_id + 2], tet_neigh[t2_id + 3] };

    for (int i = 0; i < 4; i++) if ((ng2[i] >> 2) != t1) tet_neigh[ng2[i]] = t1_id + i;
    for (int i = 0; i < 4; i++) if ((ng1[i] >> 2) != t2) tet_neigh[ng1[i]] = t2_id + i;

    for (int i = 0; i < 4; i++)
        if ((ng2[i] >> 2) != t1) tet_neigh[t1_id + i] = tet_neigh[t2_id + i];
        else tet_neigh[t1_id + i] = (tet_neigh[t2_id + i] & 3) + (t2 << 2);

    for (int i = 0; i < 4; i++)
        if ((ng1[i] >> 2) != t2) tet_neigh[t2_id + i] = ng1[i];
        else tet_neigh[t2_id + i] = (ng1[i] & 3) + (t1 << 2);
}

size_t TetMesh::markInnerTets(const std::vector<bool>& cornerMask, 
                                uint64_t single_start) {
                                    
    std::vector<uint64_t> C;

    // All ghosts are DT_OUT
    for (size_t i = 0; i < numTets(); i++)
        mark_tetrahedra[i] = (isGhost(i)) ? DT_OUT : DT_UNKNOWN;

    if (single_start != UINT64_MAX) C.push_back(single_start);
    else for (size_t i = 0; i < numTets(); i++)
        if (mark_tetrahedra[i] == DT_OUT) C.push_back(i);

    for (size_t i = 0; i < C.size(); i++) {
        uint64_t t = C[i], t4 = C[i] << 2;
        for (int j = 0; j < 4; j++) {
            const uint64_t neigh_corn = tet_neigh[t4 + j];
            const uint64_t n = neigh_corn >> 2;
            if (mark_tetrahedra[n] == DT_UNKNOWN) {
                if (!cornerMask[ neigh_corn ]) {
                    // constrained surface is not crossed
                    mark_tetrahedra[n] = mark_tetrahedra[t];
                }
                else {
                    // constrained surface is crossed
                    mark_tetrahedra[n] = 
                        (mark_tetrahedra[t] == DT_IN) ? DT_OUT : DT_IN;
                }
                C.push_back(n);
            }
        }
    }

    return std::count(mark_tetrahedra.begin(), mark_tetrahedra.end(), DT_IN);
}

bool TetMesh::hasBadSnappedOrientations(size_t& num_flipped, size_t& num_flattened) const {
    const uint32_t* tn = tet_node.data();
    const uint32_t* end = tn + tet_node.size();
    num_flipped = num_flattened = 0;
    explicitPoint v[4];
    while (tn < end) {
        if (tn[3] != INFINITE_VERTEX) {
            for (int i = 0; i < 4; i++) {
                const pointType* p = vertices[tn[i]];
                if (p->isExplicit3D()) v[i] = p->toExplicit3D();
                else p->apapExplicit(v[i]);
            }
            const int o = pointType::orient3D(v[0], v[1], v[2], v[3]);
            if (o > 0) num_flipped++;
            else if (o == 0) num_flattened++;
        }
        tn += 4;
    }

    return (num_flipped || num_flattened);
}

void TetMesh::checkMesh(bool checkDelaunay) const {
    size_t i;
    const uint32_t num_vertices = (uint32_t)vertices.size();
    // Check tet nodes	
    for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
        const uint32_t* tn = tet_node.data() + i * 4;
        if (tn[0] >= num_vertices){ 
            //ip_error("Wrong tet node!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong tet node!\n"; exit(1);
        }
        if (tn[1] >= num_vertices){ 
            //ip_error("Wrong tet node!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong tet node!\n"; exit(1);
        }
        if (tn[2] >= num_vertices){ 
            //ip_error("Wrong tet node!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong tet node!\n"; exit(1);
        }
        if (tn[3] != INFINITE_VERTEX && tet_node[i * 4 + 3] >= num_vertices){ 
            //ip_error("Wrong tet node!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong tet node!\n"; exit(1);
        }
        if (tn[0] == tn[1] || tn[0] == tn[2] || tn[0] == tn[3]
            || tn[1] == tn[2] || tn[1] == tn[3] || tn[2] == tn[3]){ 
            //ip_error("Wrong tet node indexes!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong tet node!\n"; exit(1);
        }
    }

    // Check neighbors	
    for (i = 0; i < numTets() * 4; i++) if (!isToDelete(i))
        if (tet_neigh[i] >= tet_neigh.size() || tet_neigh[tet_neigh[i]] != i){
            //ip_error("Wrong neighbor!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong neighbor!\n"; exit(1);
        }

    // Check neighbor-node coherence
    for (i = 0; i < numTets() * 4; i++) if (!isToDelete(i)) {
        if (tetHasVertex(tet_neigh[i] >> 2, tet_node[i])){
            //ip_error("Incoherent neighbor!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Incoherent neighbor!\n"; exit(1);
        }
        else {
            uint32_t v[3];
            getFaceVertices(i, v);
            if (!tetHasVertex(tet_neigh[i] >> 2, v[0])){ 
                //ip_error("Incoherent face at neighbors!\n");
                std::cout<<"[delaunay.hpp - checkMesh()] Incoherent face at neighbors!\n"; exit(1);
            }
        }
    }

    // Check vt*	
    for (i = 0; i < num_vertices; i++) if (inc_tet[i]!=UINT64_MAX) {
        if (inc_tet[i] >= numTets()){
            //ip_error("Wrong vt* (out of range)!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong vt* (out of range)!\n"; exit(1);
        }
        if (isGhost(inc_tet[i])){
            ip_error("Wrong vt* (ghost tet)!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong vt* (ghost tet)!\n"; exit(1);
        }
        const uint32_t* tn = tet_node.data() + inc_tet[i] * 4;
        if (tn[0] != i && tn[1] != i && tn[2] != i && tn[3] != i){
            //ip_error("Wrong vt*!\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Wrong vt*!\n"; exit(1);
        }
    }

    // Check marks
    //for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2))
    //    if (mark_tetrahedra[i])
    //        ip_error("Marked tet\n");

    // Check geometry
    for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
        const uint32_t* tn = tet_node.data() + i * 4;
        if (tn[3] != INFINITE_VERTEX && vOrient3D(tn[0], tn[1], tn[2], tn[3]) <= 0){ 
            //ip_error("Inverted/degn tet\n");
            std::cout<<"[delaunay.hpp - checkMesh()] Inverted/degn tet!\n"; exit(1);
        }
    }

    if (checkDelaunay) {
        for (size_t i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
            const uint32_t* n = tet_node.data() + (i * 4);
            if (n[3] == INFINITE_VERTEX) continue;
            for (int j = 0; j < 4; j++) {
                uint32_t ov = tet_node[tet_neigh[i * 4 + j]];
                if (ov != INFINITE_VERTEX && vertexInTetSphere(n, ov) > 0){ 
                    //ip_error("Non delaunay\n");
                    std::cout<<"[delaunay.hpp - checkMesh()] Non delaunay!\n"; exit(1);
                }
            }
        }
    }

    printf("checkMesh passed\n");
}

uint32_t TetMesh::findEncroachingPoint(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    static std::vector<uint64_t> enc_queue; // Static to avoid reallocation upon each call

    // Start collecting tetrahedra incident at the endpoints
    VT(ep0, enc_queue);

    for (uint64_t j : enc_queue) mark_Tet_1(j);

    const vector3d p0 = vertices[ep0];
    const vector3d p1 = vertices[ep1];
    const double eslen = (p0 - p1).sq_length();

    vector3d ep;
    uint32_t enc_pt_i = UINT32_MAX;

    marked_vertex[ep0] = marked_vertex[ep1] = 1;

    // Collect all encroaching points while expanding around insphere vertices
    for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
        const uint64_t tet = enc_queue[ti];
        const uint64_t tb = tet << 2;

        // Check each tet vertex for 'isphereness' and keep track of the one with largest sphere
        const uint32_t* tn = tet_node.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint32_t ui = tn[i];
            if (!marked_vertex[ui]) {              
                const vector3d& pui = vertices[ui];
                if (((pui - p0).sq_length() + (pui - p1).sq_length()) <= eslen) {
                    marked_vertex[ui] = 1;
                    if (enc_pt_i == UINT32_MAX || vector3d::hasLargerSphere(p0, p1, pui, ep)) {
                        ep = pui; enc_pt_i = ui;
                        tet_e = tb;
                    }
                } 
                else marked_vertex[ui] = 2;
            }
        }

        const int nvmask[] = { (marked_vertex[tn[0]] == 1), (marked_vertex[tn[1]] == 1), (marked_vertex[tn[2]] == 1), (marked_vertex[tn[3]] == 1) };
        const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

        // Expand on adjacent tets if at least one common vertex is insphere
        const uint64_t* tg = tet_neigh.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint64_t nc = tg[i];
            const uint64_t n = nc >> 2;
            if (is_marked_Tet_1(n)==2 || tet_node[nc] == INFINITE_VERTEX) continue;
            const int totmarked = totmarkeda - nvmask[i];
            if (totmarked) {
                mark_Tet_1(n);
                enc_queue.push_back(n);
            }
        }
    }

    // Clear all marks
    marked_vertex[ep0] = marked_vertex[ep1] = 0;
    for (uint64_t j : enc_queue) {
        unmark_Tet_1(j);
        j <<= 2;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j]] = 0;
    }
    enc_queue.clear();

    return enc_pt_i;
}

bool TetMesh::is_constrained_edge(uint32_t ep0, uint32_t ep1, const std::vector<bool>& constr_tri_asCorners){
    std::vector<uint64_t> vt; VT(ep0,vt);
    for(uint64_t t : vt) {
        // search for constrained faces incident at ep0 
        // (i.e. with opposite node different from ep0),
        // that have also ep1 as face vertex.
        for(size_t k=0; k<4; k++){
            uint64_t c = t*4+k;
            if(tet_node[c] != ep0 && constr_tri_asCorners[c]){
                uint32_t fv[3]; getFaceVertices(c, fv);
                if(fv[0]==ep1 || fv[1]==ep1 || fv[2]==ep1) return true;
            }
        }
    }
    return false;
}

// it is ASSUMED that input vertices are all explicit while 
// steiner points are all LNC (of explicit vrts)
inline bool TetMesh::are_vrts_connected_in_inputPLC(uint32_t v, uint32_t w,
                                                    const inputPLC& plc) const {

    const pointType* v_pt = vertices[v];
    const pointType* w_pt = vertices[w];

    bool v_is_steiner = !v_pt->isExplicit3D();
    bool w_is_steiner = !w_pt->isExplicit3D();

    assert( (v_is_steiner && v_pt->isLNC()) || v_pt->isExplicit3D() );
    assert( (w_is_steiner && w_pt->isLNC()) || w_pt->isExplicit3D() );
    
    // Two input vertices are not connected only if no input triangle edges
    // connect them.
    if(!v_is_steiner && !w_is_steiner) return plc.are_vertices_connected(v, w);
    if( v_is_steiner && !w_is_steiner ){
        // If 'v' is a point inserted on a segment incident at 'w',
        // then they are "connected", i.e. belong to incident entities.
        const explicitPoint3D* pv = &v_pt->toLNC().P().toExplicit3D();
        const explicitPoint3D* qv = &v_pt->toLNC().Q().toExplicit3D();
        const explicitPoint3D* we = &w_pt->toExplicit3D();
        return (*we == *pv || *we == *qv);
    }
    
    if(!v_is_steiner && w_is_steiner){
        // If 'w' is a point inserted on a segment incident at 'v',
        // then they are "connected", i.e. belong to incident entities.
        const explicitPoint3D* pw = &w_pt->toLNC().P().toExplicit3D();
        const explicitPoint3D* qw = &w_pt->toLNC().Q().toExplicit3D();
        const explicitPoint3D* ve = &v_pt->toExplicit3D();
        return (*ve == *pw || *ve == *qw);
    }
        
    // both are Steiner vertices
    // If 'v' and 'w' have been inserted on two segments sharing at least an
    // endpoint, then they are "connected", i.e. belong to incident entities.
    const explicitPoint3D* pv = &v_pt->toLNC().P().toExplicit3D();
    const explicitPoint3D* qv = &v_pt->toLNC().Q().toExplicit3D();
    const explicitPoint3D* pw = &w_pt->toLNC().P().toExplicit3D();
    const explicitPoint3D* qw = &w_pt->toLNC().Q().toExplicit3D();
    return (*pv == *pw || *pv == *qw || *qv == *pw || *qv == *qw);

}

bool TetMesh::tet_intersects_ballCR(uint64_t tet, uint32_t vi, double rsq){
    
    const vector3d c( vertices[vi] );

    uint64_t tb = tet << 2;
    vector3d p0(vertices[tet_node[tb  ]]);
    vector3d p1(vertices[tet_node[tb+1]]);
    vector3d p2(vertices[tet_node[tb+2]]);
    vector3d p3(vertices[tet_node[tb+3]]);

    // if a tet-vertex is closer to 'vi' less than sqrt(rsq)
    if( c.dist_sq(p0) <= rsq ) return true;
    if( c.dist_sq(p1) <= rsq ) return true;
    if( c.dist_sq(p2) <= rsq ) return true;
    if( c.dist_sq(p3) <= rsq ) return true;
    
    // if a tet-edge is closer to 'vi' less than sqrt(rsq) 
    if( c.sq_dist_segment(p0, p1) <= rsq ) return true;
    if( c.sq_dist_segment(p0, p2) <= rsq ) return true;
    if( c.sq_dist_segment(p0, p3) <= rsq ) return true;
    if( c.sq_dist_segment(p1, p2) <= rsq ) return true;
    if( c.sq_dist_segment(p1, p3) <= rsq ) return true;
    if( c.sq_dist_segment(p2, p3) <= rsq ) return true;

    // if a tet-face is closer to 'vi' less than sqrt(rsq)
    if( c.sq_dist_triangle(p0, p1, p2) <= rsq ) return true;
    if( c.sq_dist_triangle(p0, p1, p3) <= rsq ) return true;
    if( c.sq_dist_triangle(p0, p2, p3) <= rsq ) return true;
    if( c.sq_dist_triangle(p1, p2, p3) <= rsq ) return true;

    // all the sub-simplex are out from the ball
    return false;
}

// UNEFFICIENT
double TetMesh::compute_inputPLC_LFS(inputPLC& plc,
                            const std::vector<bool>& constr_tri_asCorners ) {
    const std::vector<double>& input_coords = plc.get_coordinates();
    const std::vector<uint32_t>& input_tris = plc.get_triangle_vertices();

    // Preliminary: comput input triangulation VT relation
    plc.init_vt_rel();
    
    // to compute efficiently the following minimum distances we need the CDT 
    // struture related to the PLC. The adjacencies of the PLC itself are not 
    // sufficient since closest features may not be adjacent.

    // Here we ASSUME that input vertices are all explicit while steiner points 
    // are all LNC (of explicit vrts)

    // First thing is to "measure" all input triangles, to detect the minimum 
    // between the shortest input triangle edge and the shortest input triangle 
    // heigh. This will be the initialization value of the LFS.

    // 'md' is stores the minimum measured distance (must be initialized)
    double md = vector3d(vertices[0]).dist_sq(vector3d(vertices[1]));

    uint32_t va3, vb3, vc3;
    for(size_t i = 0; i < (size_t)plc.numTriangles(); i++) {
        va3 = input_tris[i*3   ] * 3;
        vb3 = input_tris[i*3 +1] * 3;
        vc3 = input_tris[i*3 +2] * 3;
        vector3d Oa(input_coords[va3], input_coords[va3+1], input_coords[va3+2]);
        vector3d Ob(input_coords[vb3], input_coords[vb3+1], input_coords[vb3+2]);
        vector3d Oc(input_coords[vc3], input_coords[vc3+1], input_coords[vc3+2]);
        md = std::min(md, Oa.dist_sq(Ob));
        md = std::min(md, Ob.dist_sq(Oc));
        md = std::min(md, Oc.dist_sq(Oa));
        md = std::min(md, Oa.sq_dist_segment(Ob, Oc));
        md = std::min(md, Ob.sq_dist_segment(Oc, Oa));
        md = std::min(md, Oc.sq_dist_segment(Oa, Ob));
    }

    // Now have to consider all the features that are not connected by input
    // triangles.
    // Given a mesh-vertex 'vi' we have to consider the set of all its incident 
    // tetrahedra (VT) and we have to add all the other tetrahdera that ar 
    // progressivelly adjacent and inside a ball of radius R (which is the 
    // current minimum distance 'md') centered in V.
    // Then we have to chek the distance between constrained simplexes incident 
    // at 'vi' and "non-connected" and constrained simplexes contained in the  
    // set of tetrahedra.
    
    // Consider each mesh vertex
    for(uint32_t vi = 0; vi < numVertices(); vi++) {
        std::vector<uint64_t> vt; VT(vi, vt); // tets incident at vi
        vector3d Ov(vertices[vi]);

        // Growing VT untill all tetrahedra near to 'vi' than distance 'md'
        // are included. 
        for(uint64_t t : vt) mark_Tet_30(t);
        for(size_t i=0; i<vt.size(); i++){
            uint64_t t = vt[i];
            // get adjacents tets and collect them if intersects the ball
            const uint64_t* adj_tet_corn = getTetNeighs(t << 2);
            for(size_t j=0; j<4; j++) if( !isGhost(adj_tet_corn[j] >> 2) ) {
                uint64_t tj = adj_tet_corn[j] >> 2;
                if(!is_marked_Tet_30(tj) && tet_intersects_ballCR(tj, vi, md)){
                    vt.push_back(tj); 
                    mark_Tet_30(tj);
                }
            }
        }
        for(uint64_t t : vt) unmark_Tet_30(t);
        
        // Computing the minimum distance relative to vi
        uint64_t tb;
        for(uint64_t t : vt) {

            mark_Tet_30(t); // visited: avoid evaluate faces two times
            tb = t << 2;
            for(uint64_t c = tb; c < tb+4; c++) if(constr_tri_asCorners[c]){ 
                // c is tet-corner identifing a constrained tri

                // When the adjacent tetrahedron has been already visited: 
                // skip to next corner
                if(is_marked_Tet_30( tet_neigh[c] >> 2 )) continue;

                uint32_t fv[3]; getFaceVertices(c, fv);

                if(fv[0]!=vi && fv[1]!=vi && fv[2]!=vi) { 

                    // If 'vi' and a face vertex belong to two incident 
                    // input PLC segments (are "connected"), then they both 
                    // belong to the same triangle. 
                    // (Measured during initialization) 
                    if( are_vrts_connected_in_inputPLC(vi, fv[0], plc) ||
                        are_vrts_connected_in_inputPLC(vi, fv[1], plc) ||
                        are_vrts_connected_in_inputPLC(vi, fv[2], plc) ) 
                            continue;

                    vector3d Of0(vertices[fv[0]]);
                    vector3d Of1(vertices[fv[1]]);
                    vector3d Of2(vertices[fv[2]]);
                    vector3d Oc = Ov.closestPointTriangle(Of0, Of1, Of2);
                    md = std::min(md, Ov.dist_sq( Oc ));
                    // contribute to vrt-edge, vrt-face distances
                }
                else {                  
                    // tet-faces incident at vi

                    // We have to check both edges <vi,u[0]>, <vi,u[1]>
                    uint32_t u[] = {fv[0], fv[1]};
                    if(     fv[0] == vi){ u[0] = fv[2]; }
                    else if(fv[1] == vi){ u[1] = fv[2]; }

                    for(size_t k=0; k<2; k++){
                        // consider the contrained edge <vi, uk>
                        uint32_t e_vuk[] = { vi , u[k] }, oe_vuk[2];
                        std::vector<uint64_t> et; ET(vi,u[k], et);
                        for(size_t x : et){
                            oppositeTetEdge( (x << 2), e_vuk, oe_vuk);

                            if(is_constrained_edge(oe_vuk[0], oe_vuk[1],
                                                   constr_tri_asCorners)){

                                if( are_vrts_connected_in_inputPLC(vi, oe_vuk[0], plc) ||
                                    are_vrts_connected_in_inputPLC(vi, oe_vuk[1], plc) ||
                                    are_vrts_connected_in_inputPLC(u[k], oe_vuk[0], plc) ||
                                    are_vrts_connected_in_inputPLC(u[k], oe_vuk[1], plc) ) continue;

                                double d = Ov.sq_dist_seg_seg(vertices[u[k]], 
                                                          vertices[oe_vuk[0]],
                                                          vertices[oe_vuk[1]]);
                                md = std::min(md, d);
                            }
                        }
                    }
                }
            }
        }
        for(uint64_t t : vt) unmark_Tet_30(t);
    }

    plc.clear_vt_rel();
    return sqrt(md);
}

void TetMesh::set_inputPLC_LFS( inputPLC& plc,
                                const std::vector<bool>& constr_tri_asCorners){
    inputPLC_LFS = compute_inputPLC_LFS(plc, constr_tri_asCorners);
}

void TetMesh::mark_small_tris(  inputPLC& plc,
                                const std::vector<bool>& constr_tri_asCorners,
                                const std::vector<uint32_t>& corner_PLCface_map,
                                const std::vector<PLCface>& PLC_faces,
                                double bbox_rel_toll ) {

    const std::vector<double>& input_coords = plc.get_coordinates();
    const std::vector<uint32_t>& input_tris = plc.get_triangle_vertices();
    double abs_toll = bbox_rel_toll * plc.get_BBox_diag();

    // Preliminary: compute input triangulation VT relation
    plc.init_vt_rel();

    // to spot efficiently which input triangles are undersized we need the CDT 
    // struture related to the PLC. The adjacencies of the PLC itself are not 
    // sufficient since closest features may not be adjacent.

    // Here we ASSUME that input vertices are all explicit while steiner points 
    // are all LNC (of explicit vrts)

    double toll_sq = abs_toll * abs_toll;
    uint32_t va, vb, vc, va3, vb3, vc3;
    for(size_t i = 0; i < (size_t)plc.numTriangles(); i++){ 
        if(!plc.is_tri_ignored(i)) continue;
        va = input_tris[i*3   ];
        vb = input_tris[i*3 +1];
        vc = input_tris[i*3 +2];
        va3 = va * 3;  vb3 = vb * 3;  vc3 = vc * 3;
        vector3d Oa(input_coords[va3], input_coords[va3+1], input_coords[va3+2]);
        vector3d Ob(input_coords[vb3], input_coords[vb3+1], input_coords[vb3+2]);
        vector3d Oc(input_coords[vc3], input_coords[vc3+1], input_coords[vc3+2]);
        if(Oa.dist_sq(Ob) < toll_sq) plc.ignore_tri(i);
        if(Oa.dist_sq(Oc) < toll_sq) plc.ignore_tri(i);
        if(Ob.dist_sq(Oc) < toll_sq) plc.ignore_tri(i);
        if(Oa.sq_dist_segment(Ob, Oc) < toll_sq){ 
            plc.ignore_all_tris_at(va);
            plc.ignore_all_tris_at(vb, vc);
        } 
        if(Ob.sq_dist_segment(Oc, Oa) < toll_sq){ 
            plc.ignore_all_tris_at(vb);
            plc.ignore_all_tris_at(vc, va);
        }  
        if(Oc.sq_dist_segment(Oa, Ob) < toll_sq){ 
            plc.ignore_all_tris_at(vc);
            plc.ignore_all_tris_at(va, vb);
        } 
    }

    // Now consider all the features that are not connected by input triangles.
    // Given a mesh-vertex 'vi' (it is always an input point or a Steiner point 
    // on an input edge):
    // (1) collect the set "S" of all its incident tetrahedra (VT), 
    // (2) add to "S" all the other tetrahdera inside a ball of radius  
    //     'abs_toll' centered in 'vi',
    // (3) chek the distances between constrained simplexes incident at 'vi' and
    //     "non-connected" and constrained simplexes contained in "S".
    
    for(uint32_t vi=0; vi<numVertices(); vi++) {
        std::vector<uint64_t> vt; VT(vi,vt); // tets incident at vi
        vector3d Ov( vertices[vi] );

        // Growing VT untill all tetrahedra near to 'vi' than distance 
        // 'abs_toll' are included. 
        for(uint64_t t : vt) mark_Tet_30(t);
        for(size_t i=0; i<vt.size(); i++){
            // get adjacents tets and collect them if intersects the ball
            const uint64_t* adj_tet_corn = getTetNeighs(vt[i] << 2);
            for(size_t j=0; j<4; j++) {
                uint64_t tj = adj_tet_corn[j] >> 2;
                if( isGhost(tj) || is_marked_Tet_30(tj) ) continue;
                if( tet_intersects_ballCR(tj, vi, toll_sq) ) {
                    vt.push_back(tj); 
                    mark_Tet_30(tj);
                }
            }
        }
        for(uint64_t t : vt) unmark_Tet_30(t);
        
        // Mark as "to ignore" all input triagles that has a simplex not 
        // connected to 'vi' inside growed-vt.
        uint64_t tb;
        for(uint64_t t : vt) {

            mark_Tet_30(t); // visited: avoid evaluate faces two times
            tb = t << 2;
            for(uint64_t c = tb; c < tb+4; c++) if(constr_tri_asCorners[c]){ 
                // c is tet-corner identifing a constrained tri

                // When the adjacent tetrahedron has been already visited: 
                // skip to next corner
                if(is_marked_Tet_30( tet_neigh[c] >> 2 )) continue;

                uint32_t fv[3]; getFaceVertices(c, fv);

                if(fv[0]!=vi && fv[1]!=vi && fv[2]!=vi) { 
                    // If 'vi' and a face vertex belong to two incident input
                    // PLC segment (are "connected"), then they both belong
                    // to the same input triangle. (Already measured) 
                    if(are_vrts_connected_in_inputPLC(vi, fv[0], plc)) continue;
                    if(are_vrts_connected_in_inputPLC(vi, fv[1], plc)) continue;
                    if(are_vrts_connected_in_inputPLC(vi, fv[2], plc)) continue;

                    vector3d Of0(vertices[fv[0]]);
                    vector3d Of1(vertices[fv[1]]);
                    vector3d Of2(vertices[fv[2]]);
                    if(Ov.sq_dist_triangle(Of0, Of1, Of2) < toll_sq){ 
                        const PLCface& f = PLC_faces[ corner_PLCface_map[c] ];
                        for(uint32_t tfi : f.triangles) plc.ignore_tri(tfi); 
                    }
                }
                else {                  
                    // tet-faces incident at vi

                    // We have to check both edges <vi,u[0]>, <vi,u[1]>
                    uint32_t u[] = {fv[0], fv[1]};
                    if(     fv[0] == vi){ u[0] = fv[2]; }
                    else if(fv[1] == vi){ u[1] = fv[2]; }

                    for(size_t k=0; k<2; k++){
                        // consider the contrained edge 'e' = <vi, uk>
                        uint32_t e[] = { vi , u[k] }, oe[2];
                        // and each oppsite edge 'oe' on tets incident at 'e'. 
                        std::vector<uint64_t> et; ET(vi,u[k], et);
                        for(size_t x : et){
                            oppositeTetEdge( (x << 2), e, oe);

                            if(!is_constrained_edge(oe[0], oe[1],
                                                constr_tri_asCorners)) continue;

                            if(are_vrts_connected_in_inputPLC(vi, oe[0], plc) ||
                               are_vrts_connected_in_inputPLC(vi, oe[1], plc) ||
                               are_vrts_connected_in_inputPLC(u[k], oe[0], plc) ||
                               are_vrts_connected_in_inputPLC(u[k], oe[1], plc) ) continue;

                            double d_sq = Ov.sq_dist_seg_seg(vertices[u[k]], 
                                                             vertices[oe[0]],
                                                             vertices[oe[1]]);

                            if(d_sq < toll_sq){ 
                                const PLCface& f = PLC_faces[ corner_PLCface_map[c] ];
                                for(uint32_t tfi : f.triangles) plc.ignore_tri(tfi); 
                            }
                        }
                    }
                }
            }
        }
        for(uint64_t t : vt) unmark_Tet_30(t);
    }
    plc.clear_vt_rel();
}