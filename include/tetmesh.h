#include "numeric_wrapper.h"
#include "fp_geometry.h"
#include <assert.h>
#include <functional>
#include <queue>
#include "getRSS.h"

typedef unsigned char maskType;

#define TETMESH_STATIC static

// Function to delete pointed elements within a vector based on a predicate
template <class T>
void destroyUnusedElements(std::vector<T>& v, const std::function<bool(T)>& f) {
	for (T& w : v) if (f(w)) { delete w; w = NULL; }
	std::erase(v, (T)NULL);
}

// Element in a graph or complex
class TetElement
{
	friend class Tetrahedrization;

private:
	maskType mask;
	const void* infos;

public:
	TetElement() : mask(0), infos(NULL) { }

	template <int i> void mark() { mask |= (1 << i); }
	template <int i> void unmark() { mask &= (~((maskType)(1 << i))); }
	template <int i> bool isMarked() const { return ((mask & ((maskType)(1 << i))) != 0); }
	maskType getMask() const { return mask; }

	void setInfo(const void* i) { infos = i; }
	const void* getInfo() const { return infos; }

protected:
	void setMask(maskType m) { mask = m; }
};



class TetEdge;
class TetFace;
class Tetrahedron;
typedef std::vector<class PLC_Segment*> PLC_Segments;

typedef std::vector<TetEdge*> TetEdges;
typedef std::vector<TetFace*> TetFaces;
typedef std::vector<Tetrahedron*> Tetrahedra;

//////////////////////////////////////////////////////////////////////////
// 
// Vertex in an Euclidean manifold complex
//
//////////////////////////////////////////////////////////////////////////

#define USE_MEMORY_POOLS


// This is used to assign a unique index to each vertex to guarantee a
// coherent and reproducible symbolic perturbation.
// Should be used for debugging. Disable in normal conditions.
#define USE_INDEX_BASED_PERTURBATION

#ifdef USE_INDEX_BASED_PERTURBATION
uint32_t _TetVertexGlobalIndex = 0;
#define GLOBAL_INDEX_DECL uint32_t _index
#define INIT_GLOBAL_INDEX _index = _TetVertexGlobalIndex++
#define DEINIT_GLOBAL_INDEX _index = _TetVertexGlobalIndex--
#define TETVERTEX_INDEX _index
#else
#define GLOBAL_INDEX_DECL
#define INIT_GLOBAL_INDEX
#define DEINIT_GLOBAL_INDEX
#define TETVERTEX_INDEX this
#endif

class TetVertex : public TetElement
{
	pointType *_p;
	
protected:
	TetEdge* e0;
	PLC_Segments* vs_rel;

	GLOBAL_INDEX_DECL;

public:
	TetVertex() : TetElement(), _p(NULL), e0(NULL), vs_rel(NULL) { INIT_GLOBAL_INDEX; }
	TetVertex(pointType* a) : TetElement(), _p(a), e0(NULL), vs_rel(NULL) { INIT_GLOBAL_INDEX; }
	TetVertex(double x, double y, double z) : TetElement(), _p(new explicitPoint(x,y,z)), e0(NULL), vs_rel(NULL) { INIT_GLOBAL_INDEX; }
	TetVertex(pointType* a, bool dont_index) : TetElement(), _p(a), e0(NULL), vs_rel(NULL) { }
	//~TetVertex(){ DEINIT_GLOBAL_INDEX; }
#ifdef USE_MEMORY_POOLS
	void* operator new (std::size_t count);
	void operator delete(void* pointer);
#endif


	uint64_t getIndex() const { return (uint64_t)TETVERTEX_INDEX; }

	void setIncidentEdge(TetEdge* e) { e0 = e; }
	TetEdge* getIncidentEdge() const { return e0; }

	void VV(std::vector<TetVertex*>&) const;
	void VE(TetEdges&) const;
	void VF(TetFaces&) const;
	void VT(Tetrahedra&) const;

	PLC_Segments* getIncidentSegments() { return vs_rel; }
	void setIncidentSegments(PLC_Segments *ss) { vs_rel = ss; }

	TetEdge* getEdge(const TetVertex*) const;

	void unlink() { e0 = NULL; }
	bool isLinked() const { return (e0 != NULL); }

	const pointType* getPoint() const { return _p; }
	void setPoint(pointType* a) { _p = a; }
};

typedef std::vector<TetVertex*> TetVertices;

inline double sqEuclideanDistance(const TetVertex* a, const TetVertex* b) {
	return sqEuclideanDistance(a->getPoint(), b->getPoint());
}

inline double euclideanDistance(const TetVertex* a, const TetVertex* b) {
	return euclideanDistance(a->getPoint(), b->getPoint());
}

inline double tetrahedronEnergy(const TetVertex* w1, const TetVertex* w2, const TetVertex* w3, const TetVertex* w4) {
	return tetrahedronEnergy(w1->getPoint(), w2->getPoint(), w3->getPoint(), w4->getPoint());
}

//////////////////////////////////////////////////////////////////////////
// 
// Abstract graphs / complexes
//
//////////////////////////////////////////////////////////////////////////

// Generic unoriented edge connecting two vertices
class AbstractUnorientedEdge {
protected:
	TetVertex* _v[2];

public:
	AbstractUnorientedEdge(TetVertex* a, TetVertex* b) : _v{ a, b } { }

	const TetVertex* v0() const { return _v[0]; }
	const TetVertex* v1() const { return _v[1]; }
	TetVertex* v0() { return _v[0]; }
	TetVertex* v1() { return _v[1]; }
	TetVertex *const *vertices() const { return _v; }

	void unlink() { _v[0] = _v[1] = NULL; }
	bool isLinked() const { return (_v[0] != NULL); }

	void invert() { std::swap(_v[0], _v[1]); }

	void setV0(TetVertex* v) { _v[0] = v; }
	void setV1(TetVertex* v) { _v[1] = v; }

	void replaceVertex(const TetVertex* a, TetVertex* b) {
		if (_v[0] == a) _v[0] = b;
		else {
			assert(_v[1] == a);
			_v[1] = b;
		}
	}

	bool hasVertex(const TetVertex* a) const { return (_v[0] == a || _v[1] == a); }

	const TetVertex* oppositeVertex(const TetVertex* a) const {
		if (_v[0] == a) return _v[1];
		assert(_v[1] == a);
		return _v[0];
	}

	bool haveCommonVertex(const AbstractUnorientedEdge* f) const {
		return f->v0() == v0() || f->v1() == v0() || f->v0() == v1() || f->v1() == v1();
	}

	TetVertex* commonVertex(const AbstractUnorientedEdge* f)
	{
		if (f->v0() == v0() || f->v1() == v0()) return v0();
		assert(f->v0() == v1() || f->v1() == v1());
		return v1();
	}

	bool operator==(const AbstractUnorientedEdge& e) const {
		return (v0() == e.v0() && v1() == e.v1()) || (v1() == e.v0() && v0() == e.v1());
	}

	// TRUE if this segment is lexicographically less than e
	bool operator<(const AbstractUnorientedEdge& e) const {
		const TetVertex* tv1 = v0(), * tv2 = v1();
		if (tv1 > tv2) std::swap(tv1, tv2);
		const TetVertex* ov1 = e.v0(), * ov2 = e.v1();
		if (ov1 > ov2) std::swap(ov1, ov2);
		return tv1 < ov1 || (tv1 == ov1 && tv2 < ov2);
	}

	bool isEncroachedBy(const pointType* p) const {
		return (pointType::dotProductSign3D(*v0()->getPoint(), *v1()->getPoint(), *p) <= 0);
	}
};

// Generic oriented triangle defined by its three edges
template <class edge_type>
class AbstractTriangle
{
	friend class DelEdge;

protected:
	edge_type _e[3];

public:
	AbstractTriangle(edge_type a, edge_type b, edge_type c) : _e{ a,b,c } { }

	edge_type e0() const { return _e[0]; }
	edge_type e1() const { return _e[1]; }
	edge_type e2() const { return _e[2]; }

	const TetVertex* v0() const { return e1()->commonVertex(e2()); }
	const TetVertex* v1() const { return e2()->commonVertex(e0()); }
	const TetVertex* v2() const { return e0()->commonVertex(e1()); }
	TetVertex* v0() { return e1()->commonVertex(e2()); }
	TetVertex* v1() { return e2()->commonVertex(e0()); }
	TetVertex* v2() { return e0()->commonVertex(e1()); }

	void invert() { std::swap(_e[0], _e[1]); }

	edge_type* const* edges() const { return _e; }

	edge_type nextEdge(const edge_type e) const
	{
		if (e == e0()) return e1();
		if (e == e1()) return e2();
		assert(e == e2());
		return e0();
	}
	edge_type prevEdge(const edge_type e) const
	{
		if (e == e0()) return e2();
		if (e == e1()) return e0();
		assert(e == e2());
		return e1();
	}

	const TetVertex* nextVertex(const TetVertex* v) const	// Slight optimization possible here.
	{
		if (v == v0()) return v1();
		if (v == v1()) return v2();
		assert(v == v2());
		return v0();
	}

	void unlink() { _e[0] = NULL; }
	bool isLinked() const { return (_e[0] != NULL); }

	void replaceEdge(const edge_type a, edge_type b)
	{
		if (a == e0()) _e[0] = b;
		else if (a == e1()) _e[1] = b;
		else {
			assert(a == e2());
			_e[2] = b;
		}
	}

	bool hasEdge(const edge_type e) const { return (e0() == e || e1() == e || e2() == e); }
	bool hasEdges(const edge_type a, const edge_type b) const { return (hasEdge(a) && hasEdge(b)); }
	bool hasVertex(const TetVertex* v) const { return (e0()->hasVertex((TetVertex*)v) || e1()->hasVertex((TetVertex*)v)); }

	const TetVertex* oppositeVertex(const edge_type e) const
	{
		if (e == e0()) return v0();
		if (e == e1()) return v1();
		assert(e == e2());
		return v2();
	}

	edge_type oppositeEdge(const TetVertex* v) const
	{
		if (v == v0()) return e0();
		if (v == v1()) return e1();
		assert(v == v2());
		return e2();
	}

	edge_type commonEdge(const AbstractTriangle* f) const
	{
		if (f->hasEdge(e0())) return e0();
		if (f->hasEdge(e1())) return e1();
		assert(f->hasEdge(e2()));
		return e2();
	}

	const TetVertex* commonVertex(const AbstractTriangle* f) const
	{
		if (f->hasVertex(v0())) return v0();
		if (f->hasVertex(v1())) return v1();
		assert(f->hasVertex(v2()));
		return v2();
	}


	// TRUE if faces have the same three vertices
	bool operator==(const AbstractTriangle& f) const
	{
		const TetVertex* a1 = v0(), * a2 = v1(), * a3 = v2();
		const TetVertex* b1 = f.v0(), * b2 = f.v1(), * b3 = f.v2();
		return ((a1 == b1 || a1 == b2 || a1 == b3) && (a2 == b1 || a2 == b2 || a2 == b3) && (a3 == b1 || a3 == b2 || a3 == b3));
	}

	bool isEncroachedBy(const pointType* p) {
		return pointType::inGabrielSphere(*p, *v0()->getPoint(), *v1()->getPoint(), *v2()->getPoint()) <= 0;
	}
};


//////////////////////////////////////////////////////////////////////////
// 
// Tetrahedrization
//
//////////////////////////////////////////////////////////////////////////

class PLC_Segment;

// Edge in a manifold tetrahedrization
class TetEdge : public TetElement, public AbstractUnorientedEdge
{
protected:
	TetFace* f0;

public:
	PLC_Segment* segment;

	TetEdge(TetVertex* a, TetVertex* b) : TetElement(), AbstractUnorientedEdge(a, b), f0(NULL), segment(NULL) { }

#ifdef USE_MEMORY_POOLS
	void* operator new (std::size_t count);
	void operator delete(void* pointer);
#endif

	void EF(TetFaces&) const;
	void ET(Tetrahedra&) const;

	void setIncidentFace(TetFace* f) { f0 = f; }
	TetFace* getIncidentFace() const { return f0; }
	TetFace* getFace(TetVertex* v) const;

	bool isEncroached() const;
};


// Triangular face in a manifold tetrahedrization
class TetFace : public TetElement, public AbstractTriangle<TetEdge *>
{
protected:
	Tetrahedron* it1, * it2;

public:
	class DelTriangle* deltri;

	TetFace(TetEdge* a, TetEdge* b, TetEdge* c) : TetElement(), AbstractTriangle(a,b,c), it1(NULL), it2(NULL), deltri(NULL) { }

#ifdef USE_MEMORY_POOLS
	void* operator new (std::size_t count);
	void operator delete(void* pointer);
#endif

	Tetrahedron* t1() const { return it1; }
	Tetrahedron* t2() const { return it2; }

	bool isInterface() const;

	void setFirstTet(Tetrahedron* t) { it1 = t; }
	void setSecondTet(Tetrahedron* t) { it2 = t; }

	bool isOnBoundary() const { return (it2 == NULL); }

	Tetrahedron* oppositeTet(const Tetrahedron* t) const {
		if (t == it1) return it2;
		assert(t == it2);
		return it1;
	}

	void replaceTet(Tetrahedron* t, Tetrahedron* n) {
		if (t == it1) it1 = n;
		else {
			assert(t == it2);
			it2 = n;
		}
	}

	bool isEncroached() const;
	bool isWeaklyEncroached() const;
};


// Tetrahedron in a manifold tetrahedrization
class Tetrahedron : public TetElement
{
protected:
	TetFace* _f[4];

public:
	bool is_internal;
	Tetrahedron(TetFace* a, TetFace* b, TetFace* c, TetFace* d) : TetElement(), _f{ a, b, c, d } { }

#ifdef USE_MEMORY_POOLS
	void* operator new (std::size_t count);
	void operator delete(void* pointer);
#endif

	TetFace* f0() const { return _f[0]; }
	TetFace* f1() const { return _f[1]; }
	TetFace* f2() const { return _f[2]; }
	TetFace* f3() const { return _f[3]; }

	TetEdge* e0() const { return _f[0]->commonEdge(_f[1]); }
	TetEdge* e1() const { return _f[0]->commonEdge(_f[2]); }
	TetEdge* e2() const { return _f[0]->commonEdge(_f[3]); }
	TetEdge* e3() const { return _f[1]->commonEdge(_f[2]); }
	TetEdge* e4() const { return _f[1]->commonEdge(_f[3]); }
	TetEdge* e5() const { return _f[2]->commonEdge(_f[3]); }

	std::pair<TetFace*, TetFace*> facesAtEdge(const TetEdge* e) const {
		if (e == e0()) return std::pair<TetFace*, TetFace*>(_f[0], _f[1]);
		if (e == e1()) return std::pair<TetFace*, TetFace*>(_f[0], _f[2]);
		if (e == e2()) return std::pair<TetFace*, TetFace*>(_f[0], _f[3]);
		if (e == e3()) return std::pair<TetFace*, TetFace*>(_f[1], _f[2]);
		if (e == e4()) return std::pair<TetFace*, TetFace*>(_f[1], _f[3]);
		assert(e == e5());
		return std::pair<TetFace*, TetFace*>(_f[2], _f[3]);
	}

	TetFace* const* faces() const { return _f; }

	const TetVertex* v0() const { return e4()->commonVertex(e5()); }
	const TetVertex* v1() const { return e1()->commonVertex(e5()); }
	const TetVertex* v2() const { return e0()->commonVertex(e4()); }
	const TetVertex* v3() const { return e0()->commonVertex(e1()); }
	TetVertex* v0() { return e4()->commonVertex(e5()); }
	TetVertex* v1() { return e1()->commonVertex(e5()); }
	TetVertex* v2() { return e0()->commonVertex(e4()); }
	TetVertex* v3() { return e0()->commonVertex(e1()); }

	Tetrahedron* t0() const { return _f[0]->oppositeTet(this); }
	Tetrahedron* t1() const { return _f[1]->oppositeTet(this); }
	Tetrahedron* t2() const { return _f[2]->oppositeTet(this); }
	Tetrahedron* t3() const { return _f[3]->oppositeTet(this); }

	TetFace* nextFace(const TetFace* f) const
	{
		if (f == _f[0]) return _f[1];
		if (f == _f[1]) return _f[2];
		if (f == _f[2]) return _f[3];
		assert(f == _f[3]);
		return _f[0];
	}

	TetFace* prevFace(const TetFace* f) const
	{
		if (f == _f[0]) return _f[3];
		if (f == _f[1]) return _f[0];
		if (f == _f[2]) return _f[1];
		assert(f == _f[3]);
		return _f[2];
	}

	TetEdge* nextEdge(const TetEdge* e) const
	{
		if (e == e0()) return e1();
		if (e == e1()) return e2();
		if (e == e2()) return e3();
		if (e == e3()) return e4();
		if (e == e4()) return e5();
		assert(e == e5());
		return e0();
	}

	const TetVertex* nextVertex(const TetVertex* v) const
	{
		if (v == v0()) return v1();
		if (v == v1()) return v2();
		if (v == v2()) return v3();
		assert(v == v3());
		return v0();
	}

	const TetVertex* oppositeVertex(const TetFace* f) const
	{
		if (f == _f[0]) return v0();
		if (f == _f[1]) return v1();
		if (f == _f[2]) return v2();
		assert(f == _f[3]);
		return v3();
	}

	TetFace* oppositeFace(const TetVertex* v) const
	{
		if (v == v0()) return _f[0];
		if (v == v1()) return _f[1];
		if (v == v2()) return _f[2];
		assert(v == v3());
		return _f[3];
	}

	TetEdge* oppositeEdge(const TetEdge* e) const
	{
		if (e == e0()) return e5();
		if (e == e1()) return e4();
		if (e == e2()) return e3();
		if (e == e3()) return e2();
		if (e == e4()) return e1();
		assert(e == e5());
		return e0();
	}

	TetFace* oppositeEdgeFace(const TetEdge* e, const TetFace* f) const {
		const std::pair<TetFace*, TetFace*> efs = facesAtEdge(e);
		if (f == efs.first) return efs.second;
		assert(f == efs.second);
		return efs.first;
	}

	bool hasFace(const TetFace* f) const { return (f == _f[0] || f == _f[1] || f == _f[2] || f == _f[3]); }

	bool hasEdge(const TetEdge* e) const { return (e == e0() || e == e1() || e == e2() || e == e3() || e == e4() || e == e5()); }

	bool hasVertex(const TetVertex* v) const { return (v == v0() || v == v1() || v == v2() || v == v3()); }

	TetFace* commonFace(const Tetrahedron* t) const
	{
		return (t->hasFace(_f[0])) ? (_f[0]) : ((t->hasFace(_f[1])) ? (_f[1]) : ((t->hasFace(_f[2])) ? (_f[2]) : ((t->hasFace(_f[3])) ? (_f[3]) : (NULL))));
	}

	void replaceFace(const TetFace* of, TetFace* nf)
	{
		if (_f[0] == of) _f[0] = nf; 
		else if (_f[1] == of) _f[1] = nf; 
		else if (_f[2] == of) _f[2] = nf; 
		else {
			assert(_f[3] == of);
			_f[3] = nf;
		}
	}

	void invert() { std::swap(_f[0], _f[1]); }

	void unlink() { _f[0] = NULL; }
	bool isLinked() const { return (_f[0] != NULL); }


	bool sameOrientation(const TetFace* f) const
	{
		TetFace* ef1 = oppositeEdgeFace(f->e0(), f);
		TetFace* nf1 = nextFace(ef1); if (nf1 == f) nf1 = nextFace(nf1);
		TetFace* ef2 = oppositeEdgeFace(f->e1(), f);
		TetFace* nf2 = nextFace(ef2); if (nf2 == f) nf2 = nextFace(nf2);
		TetFace* ef3 = oppositeEdgeFace(f->e2(), f);
		TetFace* nf3 = nextFace(ef3); if (nf3 == f) nf3 = nextFace(nf3);

		bool aso = (nf1 == ef2 && nf2 == ef3 && nf3 == ef1);

		return (f == _f[0] || f == _f[2]) ? (aso) : (!aso);
	}

	// Energy to be minimized for mesh optimization - no need to be exact
//
// This returns DBL_MAX for flipped or degenerate tets
// For a regular tetrahedron returns 3.
// For generic non-degenerate tetrahedra returns a value in the range [3, DBL_MAX]
	double tetEnergy() const { return tetrahedronEnergy(v0(), v1(), v2(), v3()); }
	void minMaxDihedral(double& min, double& max) const {
		const pointType* v[] = { v0()->getPoint(), v1()->getPoint(), v2()->getPoint(), v3()->getPoint() };
		getMinMaxTetDihedralAngles(min, max, v);
	}
	double minDihedral() const { return (180/M_PI)*acos(minTetDihedralAngleCos(v0()->getPoint(), v1()->getPoint(), v2()->getPoint(), v3()->getPoint())); }
};

#ifdef USE_MEMORY_POOLS
inline N_memory_pool TetVertex_base_pool(4096, sizeof(TetVertex));
inline N_memory_pool TetEdge_base_pool(4096, sizeof(TetEdge));
inline N_memory_pool TetFace_base_pool(4096, sizeof(TetFace));
inline N_memory_pool Tetrahedron_base_pool(4096, sizeof(Tetrahedron));

inline void* TetVertex::operator new(std::size_t count) { return TetVertex_base_pool.alloc(); }
inline void TetVertex::operator delete(void* pointer) { TetVertex_base_pool.release(pointer); }

inline void* TetEdge::operator new(std::size_t count) { return TetEdge_base_pool.alloc(); }
inline void TetEdge::operator delete(void* pointer) { TetEdge_base_pool.release(pointer); }

inline void* TetFace::operator new(std::size_t count) { return TetFace_base_pool.alloc(); }
inline void TetFace::operator delete(void* pointer) { TetFace_base_pool.release(pointer); }

inline void* Tetrahedron::operator new(std::size_t count) { return Tetrahedron_base_pool.alloc(); }
inline void Tetrahedron::operator delete(void* pointer) { Tetrahedron_base_pool.release(pointer); }

#endif


//////////////////////////////////////////////////////////////////////////
// 
// P L C
//
//////////////////////////////////////////////////////////////////////////



class PLC_Face;

class PLC_Segment : public TetElement, public AbstractUnorientedEdge
{
public:
	std::vector<PLC_Face*> inc_PLCfaces;

	PLC_Segment(TetVertex* a, TetVertex* b) : AbstractUnorientedEdge(a, b) {
//		assert(a->getPoint()->isExplicit3D() && b->getPoint()->isExplicit3D());
	}

	void setDirty() { mark<4>(); } // Set as possibly missing
	void unsetDirty() { unmark<4>(); } // Set as not missing
	bool isDirty() { return isMarked<4>(); } // Might be missing

	void setAcute() { mark<5>(); }
	bool isAcute() const { return isMarked<5>(); }

	void addIncidentFace(PLC_Face* f) { inc_PLCfaces.push_back(f); }

	bool isFlat() const;

	bool isEncroached() const;

	PLC_Segment* split(TetVertex* vm);
};

//
// Local 2D Delaunay triangulation for PLC faces
//
class DelEdge;

class DelTriangle : public TetElement, public AbstractTriangle<DelEdge*>
{
	friend class DelEdge;

public:
	bool is_internal;

	void setDirty() { mark<4>(); }
	void unsetDirty() { unmark<4>(); }
	bool isDirty() const { return isMarked<4>(); }

	DelTriangle(DelEdge* a, DelEdge* b, DelEdge* c, PLC_Face* f) : TetElement(), AbstractTriangle(a, b, c), is_internal(true) { setInfo(f); }

#ifdef USE_MEMORY_POOLS
	void* operator new (std::size_t count);
	void operator delete(void* pointer);
#endif
};

class DelEdge : public TetElement, public AbstractUnorientedEdge
{
protected:
	DelTriangle* t1, *t2;

public:
	DelEdge(TetVertex* a, TetVertex* b) : AbstractUnorientedEdge(a, b), t1(NULL), t2(NULL) { }

#ifdef USE_MEMORY_POOLS
	void* operator new (std::size_t count);
	void operator delete(void* pointer);
#endif

	DelTriangle* firstTriangle() const { return t1; }
	DelTriangle* secondTriangle() const { return t2; }

	void setFirstTriangle(DelTriangle* t) { t1 = t; }
	void setSecondTriangle(DelTriangle* t) { t2 = t; }

	bool hasTriangle(DelTriangle* t) const { return t1 == t || t2 == t; }

	void replaceTriangle(DelTriangle* ot, DelTriangle* nt) {
		if (t1 == ot) t1 = nt;
		else {
			assert(t2 == ot);
			t2 = nt;
		}
	}

	void swap();

	DelTriangle* oppositeTriangle(const DelTriangle* t) const { if (t == t1) return t2; assert(t == t2); return t1; }

	bool isOnBoundary() const { return t2 == NULL; }

	bool isOnFaceBoundary() const {
		return (t2 == NULL && t1->is_internal) || (t2 != NULL && (t1->is_internal != t2->is_internal));
	}

	DelEdge* nextOnFaceBoundary() const {
		assert(isOnFaceBoundary());
		if (t1->is_internal) {
			DelTriangle* t = t1;
			DelEdge* e = t->nextEdge((DelEdge*)this);
			while (!e->isOnFaceBoundary()) {
				t = e->oppositeTriangle(t);
				e = t->nextEdge(e);
			}
			return e;
		}
		else {
			DelTriangle* t = t2;
			DelEdge* e = t->nextEdge((DelEdge*)this);
			while (!e->isOnFaceBoundary()) {
				t = e->oppositeTriangle(t);
				e = t->nextEdge(e);
			}
			return e;
		}
	}

	DelEdge* nextOnBoundary() const {
		assert(isOnBoundary());
		DelTriangle* t = t1;
		DelEdge* e = t->nextEdge((DelEdge*)this);

		while (!e->isOnBoundary()) {
			t = e->oppositeTriangle(t);
			e = t->nextEdge(e);
		}

		return e;
	}

	DelEdge* prevOnBoundary() const {
		assert(isOnBoundary());
		DelTriangle* t = t1;
		DelEdge* e = t->prevEdge((DelEdge*)this);

		while (!e->isOnBoundary()) {
			t = e->oppositeTriangle(t);
			e = t->prevEdge(e);
		}
		return e;
	}

	bool isEncroachedBy(const pointType* p) const {
		return (pointType::dotProductSign3D(*v0()->getPoint(), *v1()->getPoint(), *p) <= 0);
	}

	bool isEncroached() const {
		bool t2in = (t2 != NULL && t2->is_internal);
		if (t1->is_internal != t2in) { // Boundary of inner tris
			if (t1->is_internal) return isEncroachedBy(t1->oppositeVertex((DelEdge *)this)->getPoint());
			if (t2 != NULL) return isEncroachedBy(t2->oppositeVertex((DelEdge*)this)->getPoint());
		}
		return false;
	}
};

typedef std::vector<DelEdge*> DelEdges;
typedef std::vector<DelTriangle*> DelTriangles;

#ifdef USE_MEMORY_POOLS
inline N_memory_pool DelEdge_base_pool(4096, sizeof(DelEdge));
inline void* DelEdge::operator new(std::size_t count) { return DelEdge_base_pool.alloc(); }
inline void DelEdge::operator delete(void* pointer) { DelEdge_base_pool.release(pointer); }

inline N_memory_pool DelTriangle_base_pool(4096, sizeof(DelTriangle));
inline void* DelTriangle::operator new(std::size_t count) { return DelTriangle_base_pool.alloc(); }
inline void DelTriangle::operator delete(void* pointer) { DelTriangle_base_pool.release(pointer); }

#endif

// All faces are initially triangles, then they are merged across flat edges
// Possible internal vertices are implicitly defined by the local 2D Delaunay triangulation.
// 
// Before merging faces, we make sure that the bounding segments are strongly Deluanay.
// This is achieved by the segment recovery process (in 3D).
// 
// Once all segments are strongly Delaunay, we may merge coplanar feces and make a local
// 2D Delaunay triangulation of the convex hull. Because segments are strongly Delaunay,
// this triangulation will be boundary conformal.
// 
// Cocircularity is handled by making sybolic perturbation in 2D coherent with that in 3D
//
//
// To merge triangular faces into one single maximal flat face
// 1) Create a triangulation out of the input triangles
// 2) Keep a map from edges to segments (edges' info field)
// 3) Navigate the boundary of the triangulation and collect segments
//
// 4) Delaunize the 2D triangulation
// 5) Recover segments in 3D
// 6) Mark segment-edges in 2D (no segment should be missing)
// 7) Tag internal triangles in 2D

class PLC_Face : public TetElement
{
	friend class PLC_Segment;
	friend class Tetrahedrization;

protected:
	// Unordered segments bounding the face
	std::vector<PLC_Segment*> segments;

	// Local 2D triangulation
	DelEdges    E;
	DelTriangles    T;

	int max_normal_component; // 1 = X, 2 = Y, 2 = Z, -1 = -X, -2 = -Y, -3 = -Z

	int orient2D(const pointType* p1, const pointType* p2, const pointType* p3) const;
	int vOrient2D(const TetVertex* v1, const TetVertex* v2, const TetVertex* v3) const;
	int symbolicPerturbation(const TetVertex* indices[4]) const;
	DelTriangle* searchTriangle(TetVertex* v) const;
	void splitAndDelaunizeEdge(TetVertex* v, DelEdge* e);
	void splitAndDelaunizeTriangle(TetVertex* v, DelTriangle* t);
	void addExternalPoint(TetVertex* v);
	void iterativeSwap(DelEdges& toswap);
	bool swapIfNotDelaunay(DelEdge* e, DelEdge* diamond[4]);
	void insertExistingVertex(TetVertex* v);
	void insertExistingVertexOnSegment(TetVertex* v, TetVertex* sv0, TetVertex* sv1);
	void delaunizeInitialTriangles();
	void check();

	bool nearlyCollinear(const pointType* ap, const pointType* bp, const pointType* cp);

public:
	bool vInDiametralSphere(const TetVertex* v1, const TetVertex* v2, const TetVertex* v3, const TetVertex* v4) const;

	const explicitPoint* ref_t[3]; // Reference triangle to create BPTs

	PLC_Face(PLC_Segment* a, PLC_Segment* b, PLC_Segment* c) : TetElement(),
		E{ new DelEdge(a->v0(), a->v1()), new DelEdge(b->v0(), b->v1()), new DelEdge(c->v0(), c->v1()) },
		T{ new DelTriangle(E[0], E[1], E[2], this) }
	{
		DelTriangle* t = T[0];
		E[0]->setFirstTriangle(t); E[0]->setInfo(a);
		E[1]->setFirstTriangle(t); E[1]->setInfo(b);
		E[2]->setFirstTriangle(t); E[2]->setInfo(c);

		assert(t->v0()->getPoint()->isExplicit3D() && t->v1()->getPoint()->isExplicit3D() && t->v2()->getPoint()->isExplicit3D());
		ref_t[0] = &(t->v0()->getPoint()->toExplicit3D());
		ref_t[1] = &(t->v1()->getPoint()->toExplicit3D());
		ref_t[2] = &(t->v2()->getPoint()->toExplicit3D());

		double x1, y1, z1, x2, y2, z2, x3, y3, z3;
		ref_t[0]->getApproxXYZCoordinates(x1, y1, z1);
		ref_t[1]->getApproxXYZCoordinates(x2, y2, z2);
		ref_t[2]->getApproxXYZCoordinates(x3, y3, z3);
		max_normal_component = pointType::maxComponentInTriangleNormal(x1, y1, z1, x2, y2, z2, x3, y3, z3) + 1;
		if (pointType::orient2D(*ref_t[0], *ref_t[1], *ref_t[2], max_normal_component - 1) < 0) max_normal_component = -max_normal_component;
	}

	PLC_Face(PLC_Segments& c) : TetElement(), ref_t{ NULL, NULL, NULL } {
		TetVertex* apex = c.front()->commonVertex(c.back());
		TetVertex* lastv = (TetVertex*)c.front()->oppositeVertex(apex);
		E.push_back(new DelEdge(apex, lastv)); E.back()->setInfo(c.front());
		for (size_t i = 1; i < c.size()-1; i++) {
			TetVertex *cur = (TetVertex*)c[i]->oppositeVertex(lastv);
			size_t pi = E.size() - 1;
			DelEdge* prev = E[pi];
			E.push_back(new DelEdge(lastv, cur)); E.back()->setInfo(c[i]);
			E.push_back(new DelEdge(cur, apex));
			DelTriangle* t = new DelTriangle(E[pi], E[pi + 1], E[pi + 2], this);
			T.push_back(t);
			if (E[pi]->firstTriangle() == NULL) E[pi]->setFirstTriangle(t); else E[pi]->setSecondTriangle(t);
			E[pi + 1]->setFirstTriangle(t);
			E[pi + 2]->setFirstTriangle(t);
			lastv = cur;
		}
		E.back()->setInfo(c.back());

		// Search a valid triangle to be used as reference
		// We assume that each triangle is either made of three explicit points
		// or has at least one BPT
		const pointType* v0, * v1, * v2;
		for (DelTriangle* t : T) {
			v0 = t->v0()->getPoint();
			v1 = t->v1()->getPoint();
			v2 = t->v2()->getPoint();
			if (genericPoint::misaligned(*v0, *v1, *v2)) {
				if (v0->isBPT()) {
					const implicitPoint3D_BPT& b = v0->toBPT();
					ref_t[0] = &b.P(); ref_t[1] = &b.R(); ref_t[2] = &b.Q();
				}
				else if (v1->isBPT()) {
					const implicitPoint3D_BPT& b = v1->toBPT();
					ref_t[0] = &b.P(); ref_t[1] = &b.R(); ref_t[2] = &b.Q();
				}
				else if (v2->isBPT()) {
					const implicitPoint3D_BPT& b = v2->toBPT();
					ref_t[0] = &b.P(); ref_t[1] = &b.R(); ref_t[2] = &b.Q();
				}
				else if (v0->isExplicit3D() && v1->isExplicit3D() && v2->isExplicit3D()) {
					ref_t[0] = &v0->toExplicit3D();
					ref_t[1] = &v1->toExplicit3D();
					ref_t[2] = &v2->toExplicit3D();
				}

				// start ADDED by Lorenzo 20/08/2024
				else if (v0->isLNC() || v1->isLNC() || v2->isLNC()) {

					const explicitPoint3D* p[6];
					size_t i = 0;
					if (v0->isLNC()) { p[i++] = &v0->toLNC().P(); p[i++] = &v0->toLNC().Q(); }
					else { p[i++] = &v0->toExplicit3D(); }
					if (v1->isLNC()) { p[i++] = &v1->toLNC().P(); p[i++] = &v1->toLNC().Q(); }
					else { p[i++] = &v1->toExplicit3D(); }
					if (v2->isLNC()) { p[i++] = &v2->toLNC().P(); p[i++] = &v2->toLNC().Q(); }
					else { p[i++] = &v2->toExplicit3D(); }
					for (size_t j = i; j < 6; j++) p[j] = NULL;

					size_t k = 1;
					for (size_t j = k; j < i; j++) {
						size_t l = 0; for (; l < k; l++) if ((p[j] == p[l])) break;
						if (l == k) { if (j != k) { std::swap(p[j], p[k]); }  k++; }
						if (k == 3) break;
					}
					assert(k == 3);
					assert(!(p[0] == p[1]) && !(p[0] == p[2]) && !(p[2] == p[1]));

					if (pointType::orient2Dyz(*p[0], *p[1], *p[2]) != pointType::orient2Dyz(*v0, *v1, *v2) ||
						pointType::orient2Dzx(*p[0], *p[1], *p[2]) != pointType::orient2Dzx(*v0, *v1, *v2) ||
						pointType::orient2Dxy(*p[0], *p[1], *p[2]) != pointType::orient2Dxy(*v0, *v1, *v2)) {
						std::swap(p[0], p[1]);
					}

					ref_t[0] = p[0];
					ref_t[1] = p[1];
					ref_t[2] = p[2];

				}
				// end ADDED by Lorenzo 20/08/2024

				else continue;
				break;
			}
		}

		assert(ref_t[0]!=NULL);

		double x1, y1, z1, x2, y2, z2, x3, y3, z3;
		ref_t[0]->getApproxXYZCoordinates(x1, y1, z1);
		ref_t[1]->getApproxXYZCoordinates(x2, y2, z2);
		ref_t[2]->getApproxXYZCoordinates(x3, y3, z3);
		max_normal_component = pointType::maxComponentInTriangleNormal(x1, y1, z1, x2, y2, z2, x3, y3, z3) + 1;
		int oref = pointType::orient2D(*ref_t[0], *ref_t[1], *ref_t[2], max_normal_component - 1);
		if (oref < 0) max_normal_component = -max_normal_component;
	}

	const DelTriangles& getTriangles() const { return T; }
	const DelEdges& getEdges() const { return E; }

	uint32_t numEdges() const { return (uint32_t)E.size(); }

	// Steal triangles and edges from f and leave it empty
	void absorb(PLC_Face* f) {
		for (DelEdge* e : f->E) if (e->getInfo()!=NULL) {
			std::vector<PLC_Face *>& ff = ((PLC_Segment*)e->getInfo())->inc_PLCfaces;
			std::replace(ff.begin(), ff.end(), f, this);
		}
		E.insert(E.end(), f->E.begin(), f->E.end());
		T.insert(T.end(), f->T.begin(), f->T.end());
		f->E.clear();
		f->T.clear();
	}

	bool empty() const { return E.empty(); }

	// Merge edges
	void mergeEdges() {
		std::sort(E.begin(), E.end(), [](DelEdge* e1, DelEdge* e2) { return *e1 < *e2; });

		for (size_t i = 1; i < E.size(); i++) {
			DelEdge* e0 = E[i - 1], *e1 = E[i];
			if (*e0 == *e1) {
				// assert e0 and e1 point to same segment
				assert(e0->getInfo() == e1->getInfo());

				// mark such segment to be deleted
				PLC_Segment* tds = ((PLC_Segment*)e0->getInfo());
				if (tds != NULL) tds->inc_PLCfaces.clear();

				// NULLize info for both e0 and e1
				e0->setInfo(NULL); e1->setInfo(NULL);

				// e1 steals e0's incident triangle t
				DelTriangle* t = e0->firstTriangle();
				assert(t->e0() == e0 || t->e1() == e0 || t->e2() == e0);

				// replace e0 with e1 in t
				t->replaceEdge(e0, e1);

				e1->setSecondTriangle(t);
				e0->setFirstTriangle(NULL);

				i++;
			}
		}
		// Remove (and delete) redundant edges
		destroyUnusedElements<DelEdge *>(E, [](DelEdge* x) { return x->firstTriangle() == NULL; });
	}

	void makeSegments() {
		PLC_Segment* s;
		for (DelEdge* e : E) if ((s = (PLC_Segment*)e->getInfo()) != NULL) {
			assert(!s->inc_PLCfaces.empty());
			segments.push_back(s);
		}
	}

	// TRUE if v0 and v1 are in this order in the face boundary
	bool hasCoherentOrientation(TetVertex* v0, TetVertex* v1) const {
		for (DelTriangle* t : T)
			if (t->v0() == v0) return (t->v1() == v1);
			else if (t->v1() == v0) return (t->v2() == v1);
			else if (t->v2() == v0) return (t->v0() == v1);
		//ip_error("PLC_Face::hasCoherentOrientation: Should not happen!\n");
		std::cout<<"[tetmesh.h - hasCoherentOrientation()]  Should not happen!\n"; exit(1);
	}

	// Flip the face orientation
	void flip() {
		for (DelEdge* e : E) e->swap();
		for (DelTriangle* t : T) t->invert();
	}

	void iterativeAbsorb() {
		if (empty()) return;

		std::vector<PLC_Face*> coplanar;
		mark<7>();
		coplanar.push_back(this);

		for (size_t i = 0; i < coplanar.size(); i++) {
			PLC_Face* f = coplanar[i];
			assert(f->E.size() == 3);

			for (DelEdge* e : f->E) {
				PLC_Segment* s = (PLC_Segment*)e->getInfo();
				if (s!=NULL && s->isFlat()) {
					PLC_Face* opp = s->inc_PLCfaces[0];
					if (opp == f) opp = s->inc_PLCfaces[1];
					if (!opp->isMarked<7>()) {
						opp->mark<7>();
						coplanar.push_back(opp);
					}
				}
			}
		}
		for (PLC_Face* f : coplanar) f->unmark<7>();

		for (size_t i = 1; i < coplanar.size(); i++) absorb(coplanar[i]);

		if (coplanar.size() > 1) mergeEdges();


		// This is UGLY!!! Don't really understand why this inversion happens....
		for (DelEdge* e : E)
			if (vOrient2D(e->v0(), e->v1(), e->firstTriangle()->oppositeVertex(e)) < 0)
				e->invert();
	}

	bool containsPoint(const pointType* p) {
		for (DelTriangle* t : T) if (t->is_internal) {
			const pointType* p0 = t->v0()->getPoint();
			const pointType* p1 = t->v1()->getPoint();
			const pointType* p2 = t->v2()->getPoint();
			if (pointType::pointInTriangle(*p, *p0, *p1, *p2)) return true;
		}
		return false;
	}

	bool pointInBorder(const pointType* p) {
		for (DelEdge* e : E) {
			const DelTriangle* t1 = e->firstTriangle(), * t2 = e->secondTriangle();
			const pointType* p0 = e->v0()->getPoint();
			const pointType* p1 = e->v1()->getPoint();
			if ((t1->is_internal && (t2 == NULL || !t2->is_internal)) ||
				(t2 != NULL && t2->is_internal && !t1->is_internal))
				if (pointType::pointInSegment(*p, *p0, *p1)) return true;
		}
		return false;
	}

	bool hasEncroachedEdges() const {
		for (DelEdge* e : E) if (e->isEncroached()) return true;
		return false;
	}

	// Fill T with a local 2D Delaunay triangulation of the vertices
	// in the bounding segments and in current triangles.
	// The triangulation covers the entire convex hull.
	// Inner triangles are set as 'is_internal'
	// In 'nv' is not NULL, that vertex is also added to the triangulation
	void delaunize(TetVertex *nv =NULL) {
		if (segments.size() == 3 && T.size() == 1 && nv==NULL) {
			DelTriangle* t = T.back();
			t->is_internal = true; // Mark as inner
			return; // Nothing else to do. Already Delaunay
		}

		if (nv == NULL) {
			delaunizeInitialTriangles();

			// Mark DelEdges corresp to segments (dumb way...)
			for (PLC_Segment* s : segments) for (DelEdge* e : E) if (*e == *s) { e->mark<7>(); break; }

			// Mark inner triangles by starting from out and swapping across marked edges
			for (DelEdge* e : E) if (e->isOnBoundary()) { // Find a boundary edge to start from
				DelTriangle* o;
				DelTriangle* t = e->firstTriangle();
				// First is internal only if e is marked
				t->is_internal = (e->isMarked<7>());
				t->mark<7>();
				DelTriangles ts = { t };
				while (!ts.empty()) {
					t = ts.back(); ts.pop_back();

					e = t->e0();
					do {
						o = e->oppositeTriangle(t);
						if (o != NULL && !o->isMarked<7>()) {
							o->is_internal = (e->isMarked<7>()) ? (!t->is_internal) : (t->is_internal);
							o->mark<7>();
							ts.push_back(o);
						}

						e = t->nextEdge(e);
					} while (e != t->e0());
				}
				break;
			}

			for (DelEdge* e : E) e->unmark<7>();
			for (DelTriangle* t : T) t->unmark<7>();
		}
		else {
			insertExistingVertex(nv);
		}
		assert(T.size() > 1);
	}
};

bool PLC_Segment::isFlat() const {
	if (inc_PLCfaces.size() != 2) return false;
	PLC_Face* f1 = inc_PLCfaces[0];
	PLC_Face* f2 = inc_PLCfaces[1];
	if (f1->T.size() != 1 || f2->T.size() != 1) return false;

	DelTriangle* t1 = f1->T.front();
	DelTriangle* t2 = f2->T.front();

	const TetVertex* o1 = t1->v0();
	if (hasVertex(o1)) o1 = t1->v1();
	if (hasVertex(o1)) o1 = t1->v2();

	const TetVertex* o2 = t2->v0();
	if (hasVertex(o2)) o2 = t2->v1();
	if (hasVertex(o2)) o2 = t2->v2();

	//const TetVertex* o1 = NULL;
	//for (DelEdge* e : f1->E) if (!(*e == *this)) {
	//	o1 = e->v0();
	//	if (hasVertex(o1)) o1 = e->v1();
	//	assert(!hasVertex(o1));
	//	if (f1->vOrient2D(e->v0(), e->v1(), o1)) break;
	//	else o1 = NULL;
	//}
	//assert(o1 != NULL);
	//const TetVertex* o2 = NULL;
	//for (DelEdge* e : f2->E) if (!(*e == *this)) {
	//	o2 = e->v0();
	//	if (hasVertex(o2)) o2 = e->v1();
	//	assert(!hasVertex(o2));
	//	if (f2->vOrient2D(e->v0(), e->v1(), o2)) break;
	//	else o2 = NULL;
	//}
	//assert(o2 != NULL);

	const pointType* p1 = v0()->getPoint();
	const pointType* p2 = v1()->getPoint();
	const pointType* po1 = o1->getPoint();
	const pointType* po2 = o2->getPoint();

	return pointType::orient3D(*p1, *p2, *po1, *po2) == 0;
}

PLC_Segment* PLC_Segment::split(TetVertex* vm) {
	PLC_Segment* s2 = new PLC_Segment(*this);

	setV1(vm);
	s2->setV0(vm);
	for (PLC_Face* f : inc_PLCfaces)
		f->segments.push_back(s2);
	return s2;
}

typedef std::vector<PLC_Segment*> PLC_Segments;
typedef std::vector<PLC_Face*> PLC_Faces;


//////////////////////////////////////////////////////////////////////////
// 
// MAIN TETRAHEDRIZATION WHICH INCLUDES THE PLC
//
//////////////////////////////////////////////////////////////////////////

int vOrient3D(const TetVertex* v1, const TetVertex* v2, const TetVertex* v3, const TetVertex* v4) {
	return pointType::orient3D(*v1->getPoint(), *v2->getPoint(), *v3->getPoint(), *v4->getPoint());
}

class Tetrahedrization
{
protected:
	// Vertices, edges, faces and tets
	TetVertices V;
	TetEdges    E;
	TetFaces    F;
	Tetrahedra  T;

	// The PLC
	PLC_Segments S;
	PLC_Faces G;

	// Dirty PLC elements
	PLC_Segments dirty_Segments;
	DelTriangles dirty_Triangles;

	// This is true only if faces in G have a local 2D Delaunay triangulation
	bool facesAreDelaunized;

public:

	Tetrahedrization() : facesAreDelaunized(false) {}

	size_t num_vertices() { return V.size(); } const
	size_t num_vertices() const { return V.size(); } const
	size_t num_tetrahedra() { return T.size(); } const
	size_t num_faces() const { return F.size(); }

	void queueDirtySegment(PLC_Segment* s) { if (!s->isDirty()) { s->setDirty(); dirty_Segments.push_back(s); } }
	void queueDirtyTriangle(DelTriangle* s) { if (!s->isDirty()) { s->setDirty(); dirty_Triangles.push_back(s); } }

	const Tetrahedra& tets() const { return T; }
	const TetVertices& vrts() const { return V; }
	const TetFaces& faces() const { return F; }
	const PLC_Segments& get_PLCsegs() const { return S; }
	const PLC_Faces& get_PLCfaces() const { return G; }

	void initWithTetgen(int num_vertices, double vertices[], int num_triangles, uint32_t triangles[], bool quality, bool no_erosion);

	double initFromVerticesAndTets(const std::vector<pointType*>& vertices, const std::vector<uint32_t>& tet_nodes) {
		initVertices(vertices);

		for (uint32_t t = 0; t < (uint32_t)tet_nodes.size(); t += 4)
			if (tet_nodes[t+3]!=INFINITE_VERTEX)
				createTet(V[tet_nodes[t   ]], V[tet_nodes[t + 2]], V[tet_nodes[t + 1]], V[tet_nodes[t + 3]]);

		for (auto* v : V) deleteFullVE(v);
		for (auto* e : E) deleteFullEF(e);

		// closest points
		double ad, cp = DBL_MAX;
		for (auto* e : E) if ((ad = euclideanDistance(e->v0(), e->v1())) < cp) cp = ad;

		return cp / euclideanDistance(V.back(), V[V.size()-8]);
	}

	// replaceVertex : add by Lorenzo 15/11/2024
	void replaceVertex(TetVertex* old_vrt_ptr, TetVertex* new_vrt_ptr){
		std::replace(V.begin(), V.end(), old_vrt_ptr, new_vrt_ptr);
	}

	void delaunizePLCFaces() {
		for (PLC_Face* f : G) f->delaunize();
		facesAreDelaunized = true;
	}

	PLC_Segment* CreatePLC_Segment(TetVertex* v1, TetVertex* v2) {
		PLC_Segments* ve = (PLC_Segments*)v1->getInfo();
		for (PLC_Segment* s : *ve) if (s->oppositeVertex(v1) == v2) return s;
		PLC_Segment* s = new PLC_Segment(v1, v2);
		S.push_back(s);
		ve->push_back(s);
		((PLC_Segments*)v2->getInfo())->push_back(s);
		return s;
	}

	void CreatePLC_Face(uint32_t v1i, uint32_t v2i, uint32_t v3i) {
		TetVertex* v1 = V[v1i], * v2 = V[v2i], * v3 = V[v3i];

		PLC_Segment* e1 = CreatePLC_Segment(v2, v3);
		PLC_Segment* e2 = CreatePLC_Segment(v3, v1);
		PLC_Segment* e3 = CreatePLC_Segment(v1, v2);
		PLC_Face* f = new PLC_Face(e1, e2, e3);
		G.push_back(f);
		e1->addIncidentFace(f);
		e2->addIncidentFace(f);
		e3->addIncidentFace(f);
	}

	void CreatePLC_Face(std::vector<uint32_t>& v) {
		//std::reverse(v.begin(), v.end());
		TetVertex* first = V[v.front()];
		TetVertex* next = V[v.back()];
		PLC_Segment* e = CreatePLC_Segment(next, first);
		PLC_Segments s; 
		s.push_back(e);
		for (size_t i = 1; i < v.size(); i++) {
			next = V[v[i]];
			e = CreatePLC_Segment(first, next);
			s.push_back(e);
			first = next;
		}

		PLC_Face* f = new PLC_Face(s);
		G.push_back(f);
		for (PLC_Segment *t : s) t->addIncidentFace(f);
	}

	void mergePLCFaces() {
		// Add external box faces
		uint32_t fb = (uint32_t)V.size() - 8;
		const uint32_t btri[] = { 0, 1, 5, 0, 5, 4, 1, 3, 7, 1, 7, 5, 5, 7, 6, 6, 4, 5, 2, 7, 3, 7, 2, 6, 1, 2, 3, 1, 0, 2, 0, 6, 2, 6, 0, 4 };
		for (uint32_t t = 0; t < 12; t++) CreatePLC_Face(btri[t * 3] + fb, btri[t * 3 + 2] + fb, btri[t * 3 + 1] + fb);

		for (TetVertex* v : V) {
			delete ((PLC_Segments*)v->getInfo());
			v->setInfo(NULL);
		}

		// iterativeAbsorb moves by adjacency across flat segments and empties
		// incident faces of flat segments
		for (PLC_Face* f : G) f->iterativeAbsorb();

		// For each face, make segments
		for (PLC_Face* f : G) f->makeSegments();

		// Delete unused segments (with no incident faces) and remove them from S
		destroyUnusedElements< PLC_Segment* >(S, [](PLC_Segment* x) { return x->inc_PLCfaces.empty(); });

		// Delete unused faces (empty edges) and remove them from G
		destroyUnusedElements< PLC_Face* >(G, [](PLC_Face* x) { return x->numEdges() == 0; });

		// Initialize the vs relation at vertices
		initVSrelation();
	}

	void addPLCFaces(const uint32_t* tris, uint32_t num_tri) {
		for (TetVertex* v : V) v->setInfo(new PLC_Segments());
		for (uint32_t t = 0; t < num_tri; t++) CreatePLC_Face(tris[t * 3], tris[t * 3 + 1], tris[t * 3 + 2]);
		mergePLCFaces();
	}

	void addPLCFaces(std::vector<std::vector<uint32_t>>& f) {
		for (TetVertex* v : V) v->setInfo(new PLC_Segments());
		for (auto& v : f) CreatePLC_Face(v);
		mergePLCFaces();
	}

	bool isMissingSegment(PLC_Segment* s) const {
		return (s->v0()->getEdge(s->v1()) == NULL);
	}

	void recoverAllSegments(bool remove_unlinked =true) {
		recoverSegmentsInVector(S);

		if (remove_unlinked) removeUnlinkedElements();
	}

	bool isMissingTriangle(DelTriangle* t) {
		TetFaces vf; t->v0()->VF(vf);
		for (TetFace* f : vf) if (f->hasVertex(t->v1()) && t->hasVertex(t->v2())) return false;
		return true;
	}

	bool vInDiametralSphere(const TetVertex* v1, const TetVertex* v2, const TetVertex* v3, const TetVertex* v4) const {
		const pointType* p1 = v1->getPoint();
		const pointType* p2 = v2->getPoint();
		const pointType* p3 = v3->getPoint();
		const pointType* p4 = v4->getPoint();
		return pointType::inGabrielSphere(*p1, *p2, *p3, *p4) <= 0;
	}


	static TetEdge* getEdge(TetVertex* v0, TetVertex* v1) {

		const TetFace* f0 = v0->getIncidentEdge()->getIncidentFace();
		Tetrahedron* n, * t = f0->t1();

		TETMESH_STATIC Tetrahedra tets;
		tets.push_back(t); t->mark<7>();
		size_t i;
		TetEdge* e;
		for (i = 0; i < tets.size(); i++) {
			t = tets[i];
			if (t->hasVertex(v1)) {
				for (Tetrahedron* y : tets) y->unmark<7>();
				tets.clear();
				e = t->e0(); if (e->hasVertex(v0) && e->hasVertex(v1)) return e;
				e = t->e1(); if (e->hasVertex(v0) && e->hasVertex(v1)) return e;
				e = t->e2(); if (e->hasVertex(v0) && e->hasVertex(v1)) return e;
				e = t->e3(); if (e->hasVertex(v0) && e->hasVertex(v1)) return e;
				e = t->e4(); if (e->hasVertex(v0) && e->hasVertex(v1)) return e;
				return t->e5();
			}
			n = t->t0(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
			n = t->t1(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
			n = t->t2(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
			n = t->t3(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
		}

		for (Tetrahedron* t : tets) t->unmark<7>();
		tets.clear();
		return NULL;
	}

	TetFace* getTriangle(TetVertex* v0, TetVertex* v1, TetVertex* v2) const {

		const TetFace* f0 = v0->getIncidentEdge()->getIncidentFace();
		Tetrahedron* n, * t = f0->t1();

		TETMESH_STATIC Tetrahedra tets;
		tets.push_back(t); t->mark<7>();
		size_t i;
		for (i = 0; i < tets.size(); i++) {
			t = tets[i];
			if (t->hasVertex(v1) && t->hasVertex(v2)) {
				for (Tetrahedron* y : tets) y->unmark<7>();
				tets.clear();
				if (t->f0()->hasVertex(v0) && t->f0()->hasVertex(v1) && t->f0()->hasVertex(v2)) return t->f0();
				if (t->f1()->hasVertex(v0) && t->f1()->hasVertex(v1) && t->f1()->hasVertex(v2)) return t->f1();
				if (t->f2()->hasVertex(v0) && t->f2()->hasVertex(v1) && t->f2()->hasVertex(v2)) return t->f2();
				return t->f3();
			}
			n = t->t0(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
			n = t->t1(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
			n = t->t2(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
			n = t->t3(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(v0)) { tets.push_back(n); n->mark<7>(); }
		}

		for (Tetrahedron* t : tets) t->unmark<7>();
		tets.clear();
		return NULL;
	}

	bool isEncroachedTriangle(DelTriangle* t) {
		TetVertex* v0 = t->v0(), * v1 = t->v1(), * v2 = t->v2();
		TetFace* f = getTriangle(v0, v1, v2);
		if (f == NULL) return true;
		return f->isEncroached();
	}

	void initVSrelation() {
		for (TetVertex* v : V) v->setIncidentSegments(new PLC_Segments);
		for (PLC_Segment* s : S) {
			s->v0()->getIncidentSegments()->push_back(s);
			s->v1()->getIncidentSegments()->push_back(s);
		}
	}

	void deleteVSRelation() {
		for (TetVertex* v : V) {
			delete v->getIncidentSegments();
			v->setIncidentSegments(NULL);
		}
	}

	void getSegmentsInLink(TetVertex* v, PLC_Segments& segs) const {
		TETMESH_STATIC TetVertices vv;
		v->VV(vv);
		v->mark<6>(); for (TetVertex* w : vv) w->mark<6>();

		for (TetVertex* w : vv) {
			PLC_Segments* ws = ((TetVertex *)w)->getIncidentSegments();
			for (PLC_Segment* s : *ws) if (!s->isMarked<6>() && s->v0()->isMarked<6>() && s->v1()->isMarked<6>()) {
				s->mark<6>();
				segs.push_back(s);
			}
		}

		v->unmark<6>(); for (TetVertex* w : vv) w->unmark<6>();
		vv.clear();
	}


	bool recoverAllDirtySegments() {
		bool rec = false;

		while (!dirty_Segments.empty()) {
			//printf("\r\r\r\r\r\r\r\r\r\r%zu                    ", dirty_Segments.size()); fflush(stdout);
			PLC_Segment* s = dirty_Segments.back();
			dirty_Segments.pop_back();
			s->unsetDirty();
			TetEdge* e = getEdge(s->v0(), s->v1());// s->v0()->getEdge(s->v1());

			// Split if segment is not in mesh or if it is encroached
			if (e == NULL || e->isEncroached()) {				
				splitSegment(s); rec = true;
			}
			else e->segment = s;
		}
		//		printf("\n");
		return rec;
	}

	bool recoverSegmentsInVector(PLC_Segments& segs) {
		// Fill dirty list
		for (PLC_Segment* s : segs) {
			TetEdge* e = getEdge(s->v0(), s->v1()); // s->v0()->getEdge(s->v1());

			// Add to dirty list if segment is not in mesh or if it is encroached
			if (e == NULL || e->isEncroached()) queueDirtySegment(s);
			// Otherwise link the tet edge to its constraining segment
			else e->segment = s;
		}

		return recoverAllDirtySegments();
	}

	void delaunizeFace(PLC_Face* f, TetVertex* v) {
		f->insertExistingVertex(v);
		for (DelTriangle* t2 : f->getTriangles()) if (t2->is_internal && t2->hasVertex(v)) queueDirtyTriangle(t2);
	}

	void delaunizeFaceOnEdge(PLC_Face* f, TetVertex* v, TetVertex *sv0, TetVertex *sv1) {
		f->insertExistingVertexOnSegment(v, sv0, sv1);
		for (DelTriangle* t2 : f->getTriangles()) if (t2->is_internal && t2->hasVertex(v)) queueDirtyTriangle(t2);
	}

	bool recoverAllDirtyTriangles(bool verbose =false) {
		bool rec = recoverAllDirtySegments();

		while (!dirty_Triangles.empty()) {
			#ifdef DISP_PROGRESS
			if (verbose){ 
				printf("\r%zu tris missing..                     ", dirty_Triangles.size()); fflush(stdout);
			}
			#endif
			DelTriangle* t = dirty_Triangles.back();
			dirty_Triangles.pop_back();
			t->unsetDirty();
			TetFace* h = getTriangle(t->v0(), t->v1(), t->v2());
			if (h == NULL  || (h->isOnBoundary() && h->isEncroached())) {
				PLC_Face* f = (PLC_Face*)t->getInfo();
				recoverDelaunayTriangle(f, t);
				rec = true;
				recoverAllDirtySegments();
			}
			else h->deltri = t;
		}
		if (verbose) printf("\n");

		return rec;
	}

	bool recoverAllFaces() {
		// Fill dirty list
		for (PLC_Face* f : G) {
			for (DelTriangle* t : f->getTriangles()) if (t->is_internal) {
				TetFace* h = getTriangle(t->v0(), t->v1(), t->v2());
				if (h == NULL  || (h->isOnBoundary() && h->isEncroached())) queueDirtyTriangle(t);
				else h->deltri = t;
			}
		}

		bool ret = recoverAllDirtyTriangles(true);
		removeUnlinkedElements();
		return ret;
	}

	int countEncroachedDeltriangles() {
		// Fill dirty list
		int et = 0;
		for (PLC_Face* f : G) {
			for (DelTriangle* t : f->getTriangles()) if (t->is_internal) {
				TetFace* h = getTriangle(t->v0(), t->v1(), t->v2());
				assert(h != NULL);
				if (h->isEncroached()) et++;
			}
		}

		return et;
	}

	bool insertNewVertex(TetVertex* vm, Tetrahedron* t0) {
		Tetrahedron* st = searchTet(vm->getPoint(), t0);
		if (st == NULL) return false;
		insertExistingVertex(vm, st);
		vm->setIncidentSegments(new PLC_Segments);
		return true;
	}

	void splitFace(PLC_Face* f, pointType* cc, TetVertex *start) {
		assert(f->containsPoint(cc));
		assert(!f->pointInBorder(cc));

		TetVertex* vm = new TetVertex(cc);
		V.push_back(vm);

		insertNewVertex(vm, start->getIncidentEdge()->getIncidentFace()->t1());

		delaunizeFace(f, vm);
	}

	void recoverDelaunayTriangle(PLC_Face* f, DelTriangle* t) {
		const pointType* v0p = t->v0()->getPoint(), * v1p = t->v1()->getPoint(), * v2p = t->v2()->getPoint();
		explicitPoint3D ccv;
		computeCircumCenter(v0p, v1p, v2p, ccv);

		for (PLC_Segment* s : f->segments) {
			// If cc encroaches upon one of the segments of f, split segment
			if (s->isEncroachedBy(&ccv)) {
				splitSegment(s);
				// The triangle might not be deleted in this case. Need to check
				const DelTriangles& T = f->getTriangles();
				if (std::find(T.begin(), T.end(), t)!=T.end()) queueDirtyTriangle(t);
				return;
			}
		}

		implicitPoint3D_BPT* cc = convertToBPT(&ccv, f->ref_t[0], f->ref_t[1], f->ref_t[2]);
		splitFace(f, cc, t->v0());
	}

	void markInternalTets() {
		Tetrahedron* t, * y;
		TetFace* f;
		for (TetFace* g : F) if (g->isOnBoundary()) { f = g; break; }
		t = f->t1();
		t->is_internal = false;
		Tetrahedra todo;
		todo.push_back(t); t->mark<7>();

		while (!todo.empty()) {
			t = todo.back(); todo.pop_back();
			f = t->f0(); y = f->oppositeTet(t); if (y != NULL && !y->isMarked<7>()) { y->mark<7>(); todo.push_back(y); y->is_internal = (f->deltri == NULL) ? (t->is_internal) : (!t->is_internal); }
			f = t->f1(); y = f->oppositeTet(t); if (y != NULL && !y->isMarked<7>()) { y->mark<7>(); todo.push_back(y); y->is_internal = (f->deltri == NULL) ? (t->is_internal) : (!t->is_internal); }
			f = t->f2(); y = f->oppositeTet(t); if (y != NULL && !y->isMarked<7>()) { y->mark<7>(); todo.push_back(y); y->is_internal = (f->deltri == NULL) ? (t->is_internal) : (!t->is_internal); }
			f = t->f3(); y = f->oppositeTet(t); if (y != NULL && !y->isMarked<7>()) { y->mark<7>(); todo.push_back(y); y->is_internal = (f->deltri == NULL) ? (t->is_internal) : (!t->is_internal); }
		}
		for (Tetrahedron* y : T) y->unmark<7>();
	}

	void checkSegmentEncroachments() const {
			for (PLC_Segment* s : S)
				if (s->isEncroached()) {
					if (s->v0()->getEdge(s->v1()) != NULL) printf("ENCROACHED SEGMENT!\n");
					else printf("MISSING SEGMENT!\n");
					assert(0);
				}
	}


	DelTriangle* getEncroachedTriangleInFace(PLC_Face* f, const pointType* p) {
		for (DelTriangle* t : f->T) if (t->is_internal && t->isEncroachedBy(p)) return t;
		return NULL;
	}

	bool encroachesPLCSegment(const pointType* p, const Tetrahedra& cavity, PLC_Segment** es) {
		TetEdges cedges;
		TetEdge* e;
		for (Tetrahedron* t : cavity) {
			e = t->e0(); if (e->segment != NULL && !e->isMarked<7>()) { e->mark<7>(); cedges.push_back(e); }
			e = t->e1(); if (e->segment != NULL && !e->isMarked<7>()) { e->mark<7>(); cedges.push_back(e); }
			e = t->e2(); if (e->segment != NULL && !e->isMarked<7>()) { e->mark<7>(); cedges.push_back(e); }
			e = t->e3(); if (e->segment != NULL && !e->isMarked<7>()) { e->mark<7>(); cedges.push_back(e); }
			e = t->e4(); if (e->segment != NULL && !e->isMarked<7>()) { e->mark<7>(); cedges.push_back(e); }
			e = t->e5(); if (e->segment != NULL && !e->isMarked<7>()) { e->mark<7>(); cedges.push_back(e); }
		}

		*es = NULL;
		for (TetEdge* e : cedges) {
			if (e->segment->isEncroachedBy(p)) {
				*es = e->segment;
				break;
			}
		}
		for (TetEdge* e : cedges) e->unmark<7>();

		// Always check for encroachment with the outer bounding box
		if ((*es)==NULL) {
			for (auto f = G.end() - 6; f != G.end(); f++) {
				for (PLC_Segment* s : (*f)->segments) if (s->inc_PLCfaces[0] == (*f) && s->isEncroachedBy(p)) { *es = s; return true; }
			}
			return false;
		}

		return ((*es) != NULL);
	}

	bool encroachesPLCTriangle(const pointType* p, const Tetrahedra& cavity, PLC_Face** ef, DelTriangle** et) {

		TetFaces cfaces;
		TetFace* f;
		for (Tetrahedron* t : cavity) {
			f = t->f0(); if (f->deltri != NULL && !f->isMarked<7>()) { f->mark<7>(); cfaces.push_back(f); }
			f = t->f1(); if (f->deltri != NULL && !f->isMarked<7>()) { f->mark<7>(); cfaces.push_back(f); }
			f = t->f2(); if (f->deltri != NULL && !f->isMarked<7>()) { f->mark<7>(); cfaces.push_back(f); }
			f = t->f3(); if (f->deltri != NULL && !f->isMarked<7>()) { f->mark<7>(); cfaces.push_back(f); }
		}

		*et = NULL;
		for (TetFace* f : cfaces) {
			if (pointType::inGabrielSphere(*p, *f->v0()->getPoint(), *f->v1()->getPoint(), *f->v2()->getPoint()) <= 0) {
				*et = f->deltri;
				*ef = (PLC_Face*)f->deltri->getInfo();
				break;
			}
		}
		for (TetFace* f : cfaces) f->unmark<7>();

		// Always check for encroachment with the outer bounding box
		if ((*et)==NULL) {
			for (auto f = G.end() - 6; f != G.end(); f++) {
				DelTriangle* t = getEncroachedTriangleInFace(*f, p);
				if (t != NULL) { *ef = *f; *et = t; return true; }
			}
			return false;
		}

		return ((*et) != NULL);
	}

	// Calculate the cost and init tet's info field with an explicitPoint of its circumcenter
	static double computeTetCost(Tetrahedron* t, bool remove_slivers) {
		double ccc[3];

		const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
		const double shortest_edge_sqlen = getTetShortestEdgeSqLength(v);
		const double sqrad = circumsphere_ludecomp(v[0], v[1], v[2], v[3], ccc);
		const double cost = sqrad / shortest_edge_sqlen;

		// The following two lines try to account for slivers.
		// If set there is no longer a guarantee of termination, but it usually works well in practice
		if (remove_slivers) {
			const double shortest_sqheight = getTetShortestHeightSqLength(v);
			if (cost < 4.0 && (52 * shortest_sqheight) < shortest_edge_sqlen) return DBL_MAX;
		}

		return cost;
	}

	static void setTetCost(Tetrahedron* t, bool remove_slivers) {
		const double ratio = computeTetCost(t, remove_slivers);
		const uint64_t tcui = *((uint64_t*)&ratio);
		t->setInfo((void*)tcui);
	}

	static double setAndGetTetCost(Tetrahedron* t, bool remove_slivers) {
		const double ratio = computeTetCost(t, remove_slivers);
		const uint64_t tcui = *((uint64_t*)&ratio);
		t->setInfo((void*)tcui);
		return ratio;
	}

	static double getTetCost(Tetrahedron* t) {
		const uint64_t tcui = (uint64_t)t->getInfo();
		const double cost = *((double*)&tcui);
		return cost;
	}

	bool optimizeOneTet(Tetrahedron* t, double off_center_tr) {
		PLC_Face* encroached_face;
		DelTriangle* encroached_triangle;
		PLC_Segment* encroached_segment;
		bool split = false;

		double ccc[3];
		const pointType* v[4] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
		getTetCircumSpherePrecise(v, ccc, off_center_tr);
		explicitPoint3D cc(ccc[0], ccc[1], ccc[2]);

		TETMESH_STATIC Tetrahedra cavity;
		TetVertex vm(&cc, true);
		getCavity(&vm, cavity, t);

		//   Check whether cc encroaches upon any segment
		if (encroachesPLCSegment(&cc, cavity, &encroached_segment)) {
			splitSegment(encroached_segment);
			recoverAllDirtyTriangles();
			split = true;	
		}
		//   Check whether cc encroaches upon any face triangle
		else if (encroachesPLCTriangle(&cc, cavity, &encroached_face, &encroached_triangle)) {
			recoverDelaunayTriangle(encroached_face, encroached_triangle);
			split = true;
			recoverAllDirtyTriangles();
		}
		//   If cc does not encroach upon any triangle, just insert it in the mesh
		else if (!cavity.empty()) {
			EXIT_ON_THRESHOLD_NUMVERTICES;
			TetVertex* vm = new TetVertex(new explicitPoint3D(cc));
			V.push_back(vm);
			retetrahedrizeCavity(vm, cavity);
			vm->setIncidentSegments(new PLC_Segments);
			split = true;
		}

		cavity.clear();
		return split;
	}

	// Add by Lorenzo 15/11/2024
	// Calculate the cost of a virtual tetrahedron
	static double computeVirtTetCost(const pointType* v0, const pointType* v1, const pointType* v2, const pointType* v3) {
		double ccc[3];
		const pointType* v[4] = { v0, v1, v2, v3 };
		const double shortest_edge_sqlen = getTetShortestEdgeSqLength(v);
		const double sqrad = circumsphere_ludecomp(v[0], v[1], v[2], v[3], ccc);
		return sqrad / shortest_edge_sqlen;
	}

	//////////////////////////
	//
	// Main optimization function
	// Assumes that the tetmesh has no encroachment.
	// 
	// Parameters:
	//	threshold_ratio: max ratio between circumsphere radius and shortest edge of a tet
	//	remove_slivers: a heuristic is used to remove slivers during the optimization
	//	use_offcenters: use off-centers instead of circumcenters (A. Ungor)
	//
	// The algorithm is guaranteed to terminate if 'remove_slivers' is not set. However,
	// in practice it most often terminates even if 'remove_slivers' is true.
	//
	//////////////////////////
	
	void optimizeTets(double threshold_ratio = 2.0, bool remove_slivers =false, bool use_offcenters =true, uint32_t max_num_vertices =UINT32_MAX) {
	
		if (num_vertices() >= max_num_vertices) {
			deleteVSRelation();
			return;
		}

		threshold_ratio *= threshold_ratio; // Make it squared to account for squares everywhere
		bool split;
		const double offcenter_tr = (use_offcenters) ? (threshold_ratio) : (0);

		// Create a priority queue for tets
		auto tcmp = [](Tetrahedron* left, Tetrahedron* right) { return getTetCost(left) < getTetCost(right); };
		std::priority_queue < Tetrahedron*, Tetrahedra, decltype(tcmp)> queue;

		// Pre-calculate a 'cost' for each tet
		for (Tetrahedron* t : T) setTetCost(t, remove_slivers);

		do {
			split = false;

			for (Tetrahedron* t : T) if (getTetCost(t) > threshold_ratio) queue.push(t);
			int cnt = 1000;
			while (!queue.empty()) {
				Tetrahedron* t = queue.top();
				if (--cnt == 0) {
					#ifdef DISP_PROGRESS
					printf("\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\rTo be checked: %zu               ", queue.size());
					fflush(stdout);
					#endif
					cnt = 1000;
				}
				queue.pop();
				if (!t->isLinked()) continue;
				size_t tets_before = T.size();

				if (optimizeOneTet(t, offcenter_tr)) {
					if (num_vertices() >= max_num_vertices) {
						split = false;
						break;
					}
					split = true;
				}
				
				for (size_t i = tets_before; i < T.size(); i++) if (T[i]->isLinked()) {
					if (setAndGetTetCost(T[i], remove_slivers) > threshold_ratio) queue.push(T[i]);
				}
			}
			removeUnlinkedElements();
			printf("\n");
		} while (split);

		for (Tetrahedron* t : T) t->setInfo(NULL);

		deleteVSRelation();

		removeUnlinkedElements();
	}

	// REPORT AND STATS

	double maxEnergy() const {
		double ae, max_energy = 0.0;
		for (Tetrahedron* t : T) if ((ae = t->tetEnergy()) > max_energy) max_energy = ae;
		return max_energy;
	}

	void minMaxDihedral(double& min, double& max) const {
		min = DBL_MAX, max = 0.0;
		for (Tetrahedron* t : T) t->minMaxDihedral(min, max);
	}

	// double minFaceAngle() const {
	// 	double m, ma = DBL_MAX;
	// 	for (TetFace* f : F)
	// 		if ((m = ::minFaceAngle(f->v0()->getPoint(), f->v1()->getPoint(), f->v2()->getPoint())) < ma) ma = m;
	// 	return ma;
	// }

	// never called
	// double minTetAngle() const {
	// 	double m, ma = DBL_MAX;
	// 	for (Tetrahedron* t : T)
	// 		if ((m = ::minTetAngle(t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint())) < ma) ma = m;
	// 	return ma;
	// }

	double minEdgeLength() const {
		double min = DBL_MAX, t_min;
		for (Tetrahedron* t : T) {
			const pointType* v[] = { t->v0()->getPoint(), t->v1()->getPoint(), t->v2()->getPoint(), t->v3()->getPoint() };
			if ((t_min = getTetShortestEdgeSqLength(v)) < min) min = t_min;
		}
		return sqrt(min);
	}

	void averageAngles(double& av_min_faceEE, double& av_max_faceEE, double& av_min_tetFF, double& av_max_tetFF) {
		av_min_faceEE = 0.0; av_max_faceEE = 0.0;
		av_min_tetFF = 0.0; av_max_tetFF = 0.0;
		double min=DBL_MAX-1.0, max=0.0, tmp;
		for (Tetrahedron* tet : T) {
			min=DBL_MAX; max=0.0;
			const pointType* v[] = { tet->v0()->getPoint(), tet->v1()->getPoint(), tet->v2()->getPoint(), tet->v3()->getPoint() };
			tmp = getDihedralAngle(v[0],v[1],v[2],v[3]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getDihedralAngle(v[0],v[2],v[3],v[1]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getDihedralAngle(v[0],v[3],v[1],v[2]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getDihedralAngle(v[1],v[2],v[0],v[3]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getDihedralAngle(v[1],v[3],v[2],v[0]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getDihedralAngle(v[2],v[3],v[0],v[1]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			if(min<0.0 || min>360.0) min = 0.0; // CORRECTION for near degenerate tetrahedra
			if(max<0.0 || max>360.0) max = 180.0; // CORRECTION for near degenerate tetrahedra
			av_min_tetFF += min;
			av_max_tetFF += max;
		}
		av_min_tetFF /= num_tetrahedra();
		av_max_tetFF /= num_tetrahedra();

		for (TetFace* tri : F) {
			min=DBL_MAX; max=0.0;
			const pointType* v[] = { tri->v0()->getPoint(), tri->v1()->getPoint(), tri->v2()->getPoint()};
			tmp = getAngle(v[0],v[1],v[2]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getAngle(v[1],v[2],v[0]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			tmp = getAngle(v[2],v[0],v[1]); if(tmp<min) min=tmp; else if(tmp>max) max=tmp;
			if(min<0.0 || min>180.0) min = 0.0; // CORRECTION for near degenerate triangles
			if(max<0.0 || max>180.0) max = 180.0; // CORRECTION for near degenerate triangles
			av_min_faceEE += min;
			av_max_faceEE += max;
		}
		av_min_faceEE /= num_faces();
		av_max_faceEE /= num_faces();
	}

	void minMaxTetAngle(double& minEE_IN, double& maxEE_IN, double& minEE_EX, double& maxEE_EX, double& minFF_IN, double& maxFF_IN, double& minFF_EX, double& maxFF_EX) const {
		minEE_IN = DBL_MAX; minEE_EX = DBL_MAX;
		maxEE_IN = 0.0; maxEE_EX = 0.0;

		minFF_IN = DBL_MAX; minFF_EX = DBL_MAX;
		maxFF_IN = 0.0; maxFF_EX = 0.0;

		for (Tetrahedron* tet : T) {
			const pointType* v[] = { tet->v0()->getPoint(), tet->v1()->getPoint(), tet->v2()->getPoint(), tet->v3()->getPoint() };

			if (tet->is_internal) {
				getMinMaxTetFaceAngles(minEE_IN, maxEE_IN, v);
				getMinMaxTetDihedralAngles(minFF_IN, maxFF_IN, v);
			}
			else {
				getMinMaxTetFaceAngles(minEE_EX, maxEE_EX, v);
				getMinMaxTetDihedralAngles(minFF_EX, maxFF_EX, v);
			}
		}

	}

	void maxTetEnergy(double& maxIN, double& maxEX) const {
		double m;
		maxIN = 0.0; maxEX = 0.0;
		for (Tetrahedron* t : T) {
			m = t->tetEnergy();
			if (t->is_internal) { if (m > maxIN) maxIN = m; }
			else { if (m > maxEX) maxEX = m; }
		}
	}

	void printReport(bool input_encloses_vol, const char *mesh_name) {
		printf("\n-- %s --\n", mesh_name);
		printf("Mesh has %zu vertices and %zu tetrahedra\n", V.size(), T.size());

		double mai, mae;
		maxTetEnergy(mai, mae);
		printf("\nMax tet energy:\n\tInternal tets: %f\n", mai);
		if(input_encloses_vol) printf("\tExternal tets: %f\n", mae);

		double minPlanIN, maxPlanIN, minPlanEX, maxPlanEX;
		double minDihedIN, maxDihedIN, minDihedEX, maxDihedEX;
		minMaxTetAngle(minPlanIN, maxPlanIN, minPlanEX, maxPlanEX, minDihedIN, maxDihedIN, minDihedEX, maxDihedEX);
		printf("\nExtremal angles (DEGs)\n");
		printf("Min dihedral angle:\n\tInternal mesh: %f\n", minDihedIN);
		if(input_encloses_vol) printf("\tExternal mesh: %f\n", minDihedEX);
		printf("Max dihedral angle:\n\tInternal mesh: %f\n", maxDihedIN);
		if(input_encloses_vol) printf("\tExternal mesh: %f\n", maxDihedEX);
		printf("Min face angle:\n\tInternal mesh: %f\n", minPlanIN);
		if(input_encloses_vol) printf("\tExternal mesh: %f\n", minPlanEX);
		printf("Max face angle:\n\tInternal mesh: %f\n", maxPlanIN);
		if(input_encloses_vol) printf("\tExternal mesh: %f\n", maxPlanEX);
		double avMinPlan, avMaxPlan, avMinDihed, avMaxDihed;
		averageAngles(avMinPlan, avMaxPlan, avMinDihed, avMaxDihed);
		printf("\nAverage angles (DEGs)\n");
		printf("Average min per-face face angle: %f\n", avMinPlan);
		printf("Average max per-face face angle: %f\n", avMaxPlan);
		printf("Average min per-tet dihedral angle: %f\n", avMinDihed);
		printf("Average max per-tet dihedral angle: %f\n", avMaxDihed);

	}
	// -----------------------------------
	void insertExistingVertex(TetVertex* vm, Tetrahedron *st) {

		TETMESH_STATIC Tetrahedra cavity;
		Tetrahedron *t0 = searchTet(vm->getPoint(), st);
		getCavity(vm, cavity, t0);

		// DEBUG start
		// if(cavity.size() > 10000){
		// 	std::cout << "\nPROGRAM ABORTED: too large cavity (" << cavity.size() << " tetrahedra).\n\n\n";
		// 	exit(13);
		// }
		// DEBUG end

		retetrahedrizeCavity(vm, cavity);
		cavity.clear();
	}

	void splitSegment(PLC_Segment* s, pointType* vmp) {
		TetVertex* vm = new TetVertex(vmp);
		V.push_back(vm);

		insertNewVertex(vm, s->v0()->getIncidentEdge()->getIncidentFace()->t1());

		PLC_Segment* ns = s->split(vm);
		S.push_back(ns);
		queueDirtySegment(s);
		ns->setDirty(); dirty_Segments.push_back(ns);

		PLC_Segments* ss = vm->getIncidentSegments();
		ss->push_back(s);
		ss->push_back(ns);

		PLC_Segments* ovs = ns->v1()->getIncidentSegments();
		std::replace(ovs->begin(), ovs->end(), s, ns);

		if (facesAreDelaunized)
			for (PLC_Face* fi : s->inc_PLCfaces) delaunizeFaceOnEdge(fi, vm, s->v0(), ns->v1());
	}

	void splitSegment(PLC_Segment* s) {
		pointType* vmp;

		const pointType* v0p = s->v0()->getPoint(), * v1p = s->v1()->getPoint();

		vmp = createMidPoint(v0p, v1p);

		if (vmp->isLNC() && (vmp->toLNC().T() <= 0.0 || vmp->toLNC().T() >= 1.0)) {

			// DEBUG
			double x0, y0, z0; v0p->getApproxXYZCoordinates(x0, y0, z0);
			double x1, y1, z1; v1p->getApproxXYZCoordinates(x1, y1, z1);
			std::cout << "v0 = " << vector3d(x0, y0, z0) << " , v1 = " << vector3d(x1, y1, z1) << "\n";
			std::cout << "t=" << vmp->toLNC().T() << "\n";
			std::cout << "of len = " << sqrt(vector3d(x0, y0, z0).dist_sq(vector3d(x1, y1, z1))) << "\n";
			std::cout << "on original < " << vmp->toLNC().P() << " , " << vmp->toLNC().Q() << " > \n";
			double xp, yp, zp; (vmp->toLNC().P()).getApproxXYZCoordinates(xp, yp, zp);
			double xq, yq, zq; (vmp->toLNC().Q()).getApproxXYZCoordinates(xq, yq, zq);
			std::cout << "original len = " << sqrt(vector3d(xp, yp, zp).dist_sq(vector3d(xq, yq, zq))) << "\n";

			// ip_error("Could not split further! Did you verify whether your PLC is non-acute?\n");
			std::cout<<"[tetmesh.h - splitSegment()]  Could not split further! Verify whether your PLC is non-acute\n"; exit(1);
		}

		if (vmp->isBPT() && (vmp->toBPT().U() <= 0.0 || vmp->toBPT().U() >= 1.0 || vmp->toBPT().V() <= 0.0 || vmp->toBPT().V() >= 1.0)) {

			// DEBUG
			double x0, y0, z0; v0p->getApproxXYZCoordinates(x0, y0, z0);
			double x1, y1, z1; v1p->getApproxXYZCoordinates(x1, y1, z1);
			std::cout << "v0 = " << vector3d(x0, y0, z0) << " , v1 = " << vector3d(x1, y1, z1) << "\n";
			std::cout << "u=" << vmp->toBPT().U() << "\n";
			std::cout << "v=" << vmp->toBPT().V() << "\n";
			std::cout << "of len = " << sqrt(vector3d(x0, y0, z0).dist_sq(vector3d(x1, y1, z1))) << "\n";
			std::cout << "on original < " << vmp->toBPT().P() << " , " << vmp->toBPT().Q() << " , " << vmp->toBPT().R() << " > \n";
			double xp, yp, zp; (vmp->toBPT().P()).getApproxXYZCoordinates(xp, yp, zp);
			double xq, yq, zq; (vmp->toBPT().Q()).getApproxXYZCoordinates(xq, yq, zq);
			double xr, yr, zr; (vmp->toBPT().R()).getApproxXYZCoordinates(xr, yr, zr);
			std::cout << "original len PQ= " << sqrt(vector3d(xp, yp, zp).dist_sq(vector3d(xq, yq, zq))) << "\n";
			std::cout << "original len QR= " << sqrt(vector3d(xr, yr, zr).dist_sq(vector3d(xq, yq, zq))) << "\n";
			std::cout << "original len RP= " << sqrt(vector3d(xp, yp, zp).dist_sq(vector3d(xr, yr, zr))) << "\n";

			// ip_error("Could not split further! Did you verify whether your PLC is non-acute?\n");
			std::cout<<"[tetmesh.h - splitSegment()]  Could not split further! Verify whether your PLC is non-acute\n"; exit(1);
		}

		assert(!pointType::coincident(*v0p, *vmp));
		assert(!pointType::coincident(*v1p, *vmp));

		splitSegment(s, vmp);
	}

	bool isTetPositive(Tetrahedron* t) {
			const pointType* v1 = t->v0()->getPoint();
			const pointType* v2 = t->v1()->getPoint();
			const pointType* v3 = t->v2()->getPoint();
			const pointType* v4 = t->v3()->getPoint();

			return (pointType::orient3D(*v1, *v2, *v3, *v4) > 0);
	}

	bool saveTET(const char* filename)
	{
		uint64_t idx = 0;
		for (TetVertex* v : V) v->setInfo((void *)idx++);
		FILE* fp;
		if ((fp = fopen(filename, "w")) == NULL) return false;

		fprintf(fp, "%zu vertices\n%zu tets\n", V.size(), T.size());
		double x, y, z;
		for (TetVertex* v : V) {
			const pointType* p = v->getPoint();
			p->getApproxXYZCoordinates(x, y, z);
			fprintf(fp, "%f %f %f\n", x, y, z);
		}
		for (Tetrahedron* t : T) {
			size_t i1 = (size_t)t->v0()->getInfo();
			size_t i2 = (size_t)t->v1()->getInfo();
			size_t i3 = (size_t)t->v2()->getInfo();
			size_t i4 = (size_t)t->v3()->getInfo();
			fprintf(fp, "4 %zu %zu %zu %zu\n", i1, i2, i3, i4);
		}

		for (TetVertex* v : V) v->setInfo(NULL);
		fclose(fp);
		return true;
	}

	bool saveMEDIT(const char* filename)
	{
		uint64_t idx = 0;
		for (TetVertex* v : V) v->setInfo((void*)idx++);
		FILE* fp;
		if ((fp = fopen(filename, "w")) == NULL) return false;

		fprintf(fp, "MeshVersionFormatted 2\nDimension\n3\n");
		fprintf(fp, "Vertices %zu\n", V.size());

		double x, y, z;
		for (TetVertex* v : V) {
			const pointType* p = v->getPoint();
			p->getApproxXYZCoordinates(x, y, z);
			fprintf(fp, "%f %f %f 1\n", x, y, z);
		}

		fprintf(fp, "Tetrahedra %zu\n", T.size());

		bool use_inner = true;
		for (int i = 0; i < 2; i++) {
			for (Tetrahedron* t : T) if (t->is_internal == use_inner) {
				size_t i1 = (size_t)t->v0()->getInfo();
				size_t i2 = (size_t)t->v1()->getInfo();
				size_t i3 = (size_t)t->v2()->getInfo();
				size_t i4 = (size_t)t->v3()->getInfo();
				fprintf(fp, "%zu %zu %zu %zu %d\n", ++i1, ++i2, ++i3, ++i4, i+1);
			}
			use_inner = false;
		}

		for (TetVertex* v : V) v->setInfo(NULL);
		fclose(fp);
		return true;
	}

	bool saveOFFBoundary(const char* filename)
	{
		size_t num_v = 0, num_t = 0;

		for (TetVertex* v : V) v->setInfo(0);
			
		for (TetFace* f : F) if (f->isOnBoundary()) {
			f->v0()->setInfo((void *)1);
			f->v1()->setInfo((void *)1);
			f->v2()->setInfo((void *)1);
			num_t++;
		}

		for (TetVertex* v : V) if (v->getInfo()) {
			num_v++;
			v->setInfo((void*)(num_v));
		}

		FILE* fp;
		if ((fp = fopen(filename, "w")) == NULL) return false;

		fprintf(fp, "OFF\n");
		fprintf(fp, "%zu %zu 0\n", num_v, num_t);

		double x, y, z;
		for (TetVertex* v : V) if (v->getInfo()) {
			const pointType* p = v->getPoint();
			p->getApproxXYZCoordinates(x, y, z);
			fprintf(fp, "%f %f %f\n", x, y, z);
		}

		size_t i1, i2, i3;
		for (TetFace* f : F) if (f->isOnBoundary()) {
			i1 = (size_t)f->v0()->getInfo();
			i2 = (size_t)f->v1()->getInfo();
			i3 = (size_t)f->v2()->getInfo();
			fprintf(fp, "3 %zu %zu %zu\n", i1-1, i2-1, i3-1);
		}

		for (TetVertex* v : V) v->setInfo(0);

		fclose(fp);
		return true;
	}

	bool saveOFFInterface(const char* filename)
	{
		size_t num_v = 0, num_t = 0;

		for (TetVertex* v : V) v->setInfo(0);

		for (TetFace* f : F) if (f->isInterface()) {
			f->v0()->setInfo((void*)1);
			f->v1()->setInfo((void*)1);
			f->v2()->setInfo((void*)1);
			num_t++;
		}

		for (TetVertex* v : V) if (v->getInfo()) {
			num_v++;
			v->setInfo((void*)(num_v));
		}

		FILE* fp;
		if ((fp = fopen(filename, "w")) == NULL) return false;

		fprintf(fp, "OFF\n");
		fprintf(fp, "%zu %zu 0\n", num_v, num_t);

		double x, y, z;
		for (TetVertex* v : V) if (v->getInfo()) {
			const pointType* p = v->getPoint();
			p->getApproxXYZCoordinates(x, y, z);
			fprintf(fp, "%f %f %f\n", x, y, z);
		}

		size_t i1, i2, i3;
		for (TetFace* f : F) if (f->isInterface()) {
			i1 = (size_t)f->v0()->getInfo();
			i2 = (size_t)f->v1()->getInfo();
			i3 = (size_t)f->v2()->getInfo();
			if (f->t1()->is_internal) fprintf(fp, "3 %zu %zu %zu\n", i1 - 1, i2 - 1, i3 - 1);
			else fprintf(fp, "3 %zu %zu %zu\n", i3 - 1, i2 - 1, i1 - 1);
		}

		for (TetVertex* v : V) v->setInfo(0);

		fclose(fp);
		return true;
	}


	void removeUnlinkedElements()
	{
		destroyUnusedElements<TetVertex*>(V, [](TetVertex* x) { return !x->isLinked(); });
		destroyUnusedElements<TetEdge *>(E, [](TetEdge* x) { return !x->isLinked(); });
		destroyUnusedElements<TetFace *>(F, [](TetFace* x) { return !x->isLinked(); });
		destroyUnusedElements<Tetrahedron *>(T, [](Tetrahedron* x) { return !x->isLinked(); });
	}

	Tetrahedron* searchTet(const pointType* p, Tetrahedron* t = NULL) {
		if (t == NULL) t = T.back();

		do {
			assert(t != NULL);
			const pointType* v1 = t->v0()->getPoint();
			const pointType* v2 = t->v1()->getPoint();
			const pointType* v3 = t->v2()->getPoint();
			const pointType* v4 = t->v3()->getPoint();
			if (pointType::orient3D(*p, *v1, *v4, *v3) < 0) t = t->t1();
			else if (pointType::orient3D(*p, *v1, *v2, *v4) < 0) t = t->t2();
			else if (pointType::orient3D(*p, *v1, *v3, *v2) < 0) t = t->t3();
			else if (pointType::orient3D(*p, *v2, *v3, *v4) < 0) t = t->t0();
			else return t;
		} while (t != NULL);

		return NULL;
	}

	int symbolicPerturbation(const TetVertex *indices[5]) const {
		int swaps = 0;
		int n = 5;
		int count;
		do {
			count = 0;
			n--;
			for (int i = 0; i < n; i++) {
				if (indices[i]->getIndex() > indices[i + 1]->getIndex()) {
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

	bool inTetCircumsphere(TetVertex* v, const Tetrahedron* t) const {
		const TetVertex* v1 = t->v0();
		const TetVertex* v2 = t->v1();
		const TetVertex* v3 = t->v2();
		const TetVertex* v4 = t->v3();
		const pointType* p1 = v1->getPoint();
		const pointType* p2 = v2->getPoint();
		const pointType* p3 = v3->getPoint();
		const pointType* p4 = v4->getPoint();
		int det = pointType::inSphere(*v->getPoint(), *p1, *p2, *p3, *p4);
		if (det == 0) {
			const TetVertex* idx[5] = { v, v1, v2, v3, v4 };
			det = symbolicPerturbation(idx);
		}
		return (det > 0);
	}

	void getCavity(TetVertex *v, Tetrahedra& cavity, Tetrahedron* t0 = NULL) {
		if (t0 == NULL) t0 = searchTet(v->getPoint());
		assert(t0 != NULL);

		Tetrahedron* s, * t;
		cavity.push_back(t0); t0->mark<7>();
		for (size_t i = 0; i < cavity.size(); i++) {
			t = cavity[i];
			s = t->t0(); if (s != NULL && !s->isMarked<7>() && inTetCircumsphere(v, s)) { cavity.push_back(s); s->mark<7>(); }
			s = t->t1(); if (s != NULL && !s->isMarked<7>() && inTetCircumsphere(v, s)) { cavity.push_back(s); s->mark<7>(); }
			s = t->t2(); if (s != NULL && !s->isMarked<7>() && inTetCircumsphere(v, s)) { cavity.push_back(s); s->mark<7>(); }
			s = t->t3(); if (s != NULL && !s->isMarked<7>() && inTetCircumsphere(v, s)) { cavity.push_back(s); s->mark<7>(); }
		}

		for (Tetrahedron* t : cavity) t->unmark<7>();
	}

	void checkFaceOrientations() const {
		for (TetFace* f : F) if (f->isLinked()) {
			assert(vOrient3D(f->t1()->oppositeVertex(f), f->v0(), f->v1(), f->v2()) > 0);
			if (f->t2()!=NULL)
			assert(vOrient3D(f->t2()->oppositeVertex(f), f->v0(), f->v2(), f->v1()) > 0);
		}
		std::cout << "checkFaceOrientations passed!\n";
	}

	void retetrahedrizeCavity(TetVertex* v, Tetrahedra& cavity) {
		TETMESH_STATIC TetFaces boundary, innerFaces;
		TETMESH_STATIC TetEdges bdEdges, innerEdges;
		TetFace* f;
		TetEdge* e;
		Tetrahedron* s;

		// Set dirty segments and triangles
		for (Tetrahedron* t : cavity) {
			PLC_Segment* q;
			DelTriangle* r;
			r = t->f0()->deltri; if (r != NULL) queueDirtyTriangle(r);
			r = t->f1()->deltri; if (r != NULL) queueDirtyTriangle(r);
			r = t->f2()->deltri; if (r != NULL) queueDirtyTriangle(r);
			r = t->f3()->deltri; if (r != NULL) queueDirtyTriangle(r);
			q = t->e0()->segment; if (q != NULL) queueDirtySegment(q);
			q = t->e1()->segment; if (q != NULL) queueDirtySegment(q);
			q = t->e2()->segment; if (q != NULL) queueDirtySegment(q);
			q = t->e3()->segment; if (q != NULL) queueDirtySegment(q);
			q = t->e4()->segment; if (q != NULL) queueDirtySegment(q);
			q = t->e5()->segment; if (q != NULL) queueDirtySegment(q);
		}

		// Mark cavity tets
		for (Tetrahedron* t : cavity) t->mark<7>();

		// Collect boundary and inner faces
		for (Tetrahedron* t : cavity) {
			f = t->f0(); s = f->oppositeTet(t);
			if (s == NULL || !s->isMarked<7>()) { boundary.push_back(f); }
			else if (!f->isMarked<7>()) { innerFaces.push_back(f); f->mark<7>(); }
			f = t->f1(); s = f->oppositeTet(t);
			if (s == NULL || !s->isMarked<7>()) { boundary.push_back(f); }
			else if (!f->isMarked<7>()) {	innerFaces.push_back(f); f->mark<7>(); }
			f = t->f2(); s = f->oppositeTet(t); 
			if (s == NULL || !s->isMarked<7>()) { boundary.push_back(f); }
			else if (!f->isMarked<7>()) { innerFaces.push_back(f); f->mark<7>(); }
			f = t->f3(); s = f->oppositeTet(t); 
			if (s == NULL || !s->isMarked<7>()) { boundary.push_back(f); }
			else if (!f->isMarked<7>()) { innerFaces.push_back(f); f->mark<7>(); }
		}

		// For each boundary face, if coplanar with v remove from boundary, move to inner and mark
		for (TetFace* fb : boundary)
			if (vOrient3D(v, fb->v0(), fb->v1(), fb->v2()) == 0) {
				if (!fb->isOnBoundary()) {
					assert(!pointType::coincident(*v->getPoint(), *fb->v0()->getPoint()));
					assert(!pointType::coincident(*v->getPoint(), *fb->v1()->getPoint()));
					assert(!pointType::coincident(*v->getPoint(), *fb->v2()->getPoint()));
					assert(fb->isOnBoundary());
				} 
				fb->mark<7>();
				innerFaces.push_back(fb);
			}
		std::erase_if(boundary, [](TetFace* fb) { return fb->isMarked<7>(); });

		// Collect boundary edges and init their full EF relation
		for (TetFace* f : boundary) {
			e = f->e0(); if (!e->isMarked<7>()) { bdEdges.push_back(e); e->mark<7>(); initFullEF(e, f); e->setIncidentFace(f); }
			else { ((TetFaces*)e->getInfo())->push_back(f); }
			e = f->e1(); if (!e->isMarked<7>()) { bdEdges.push_back(e); e->mark<7>(); initFullEF(e, f); e->setIncidentFace(f); }
			else { ((TetFaces*)e->getInfo())->push_back(f); }
			e = f->e2(); if (!e->isMarked<7>()) { bdEdges.push_back(e); e->mark<7>(); initFullEF(e, f); e->setIncidentFace(f); }
			else { ((TetFaces*)e->getInfo())->push_back(f); }
		}


		// Init full VE relation for boundary vertices
		for (TetEdge* e : bdEdges) {
			if (e->v0()->getInfo() == NULL) { initFullVE(e->v0(), e); e->v0()->setIncidentEdge(e); }
			else ((TetEdges*)e->v0()->getInfo())->push_back(e);
			if (e->v1()->getInfo() == NULL) { initFullVE(e->v1(), e); e->v1()->setIncidentEdge(e); }
			else ((TetEdges*)e->v1()->getInfo())->push_back(e);
		}

		// Unlink inner faces and edges
		for (TetFace* f : innerFaces) {
			e = f->e0(); if (!e->isMarked<7>()) e->unlink();
			e = f->e1(); if (!e->isMarked<7>()) e->unlink();
			e = f->e2(); if (!e->isMarked<7>()) e->unlink();
			f->unlink();
		}

		// Unlink old tets
		for (Tetrahedron* t : cavity) t->unlink();

		// Remove unlinked tets from FT relation
		for (TetFace* f : boundary) if (f->t1()->isMarked<7>()) f->setFirstTet(NULL); else f->setSecondTet(NULL);

		// Init full VE relation for v
		initFullVE(v);

		// Make the new tets
		for (TetFace* f : boundary) {
			assert(vOrient3D(v, f->v0(), f->v1(), f->v2()) != 0);
			if (f->t1() == NULL) 
				assert(vOrient3D(v, f->v0(), f->v1(), f->v2()) > 0);
			else
				assert(vOrient3D(v, f->v0(), f->v2(), f->v1()) > 0);

			CreateTet(v, f);
		}

		// Clear the temporary full relations
		for (TetEdge* e : bdEdges) {
			if (e->v1()->getInfo() != NULL) { deleteFullVE(e->v1()); }
			if (e->v0()->getInfo() != NULL) { deleteFullVE(e->v0()); }
		}
		for (TetEdge* e : bdEdges) { deleteFullEF(e); }

		// Clear full EF for each edge in VE(v)
		for (TetEdge* e : *((TetEdges*)v->getInfo())) { deleteFullEF(e); }

		// Clear full VE of v
		deleteFullVE(v);

		// Unmark boundary edges
		for (TetEdge* e : bdEdges) e->unmark<7>();

		boundary.clear(); innerFaces.clear();
		bdEdges.clear(); innerEdges.clear();
	}

protected:
	void initFullVE(TetVertex* v) { v->setInfo(new TetEdges()); }
	void initFullVE(TetVertex* v, TetEdge* e) { v->setInfo(new TetEdges({ e })); }
	void deleteFullVE(TetVertex* v) { delete ((TetEdges*)v->getInfo()); v->setInfo(NULL); }

	void initFullEF(TetEdge* e) { e->setInfo(new TetFaces()); }
	void initFullEF(TetEdge* e, TetFace *f) { e->setInfo(new TetFaces({f})); }
	void deleteFullEF(TetEdge* e) { delete ((TetFaces*)e->getInfo()); e->setInfo(NULL); }

	void initVertices(const std::vector<pointType*>& vertices) {
		V.reserve(vertices.size());

		// Symbolic perturbation uses vertex pointers to disambiguate, therefore it is
		// crucial that these pointers are in the same order as in the vector. That is
		// why we allocate them all here as a single vector.
		TetVertex* vptrs = new TetVertex[vertices.size()];

		uint32_t i = 0;
		for (pointType* p : vertices) {
			vptrs[i].setPoint(p);
			V.push_back(vptrs + i);
			initFullVE(V.back());
			i++;
		}
	}

	TetEdge* CreateEdge(TetVertex* v1, TetVertex* v2)
	{
		TetEdges* ve = (TetEdges*)v1->getInfo();
		for (TetEdge* e : *ve) if (e->hasVertex(v2)) return e;

		TetEdge* e = new TetEdge(v1, v2); e->setInfo(new TetFaces);
		((TetEdges*)v1->getInfo())->push_back(e); v1->setIncidentEdge(e);
		((TetEdges*)v2->getInfo())->push_back(e); v2->setIncidentEdge(e);

		E.push_back(e);
		return e;
	}

	TetFace* CreateFace(TetEdge* e1, TetEdge* e2, TetEdge* e3)
	{
		TetFaces* ef = (TetFaces*)e1->getInfo();
		for (TetFace* f : *ef) if (f->hasEdges(e2, e3)) return f;

		TetFace* f = new TetFace(e1, e2, e3);
		e1->setIncidentFace(f); ((TetFaces*)e1->getInfo())->push_back(f);
		e2->setIncidentFace(f); ((TetFaces*)e2->getInfo())->push_back(f);
		e3->setIncidentFace(f); ((TetFaces*)e3->getInfo())->push_back(f);

		F.push_back(f);
		return f;
	}

	Tetrahedron* createTet(TetVertex* v1, TetVertex* v2, TetVertex* v3, TetVertex* v4) {
		TetEdge* e1 = CreateEdge(v3, v4);
		TetEdge* e2 = CreateEdge(v2, v4);
		TetEdge* e3 = CreateEdge(v2, v3);
		TetEdge* e4 = CreateEdge(v1, v4);
		TetEdge* e5 = CreateEdge(v1, v3);
		TetEdge* e6 = CreateEdge(v1, v2);
		TetFace* f1 = CreateFace(e2, e3, e1);
		TetFace* f2 = CreateFace(e5, e4, e1);
		TetFace* f3 = CreateFace(e4, e6, e2);
		TetFace* f4 = CreateFace(e6, e5, e3);
		Tetrahedron* t = new Tetrahedron(f1, f2, f3, f4); T.push_back(t);
		if (f1->t1() == NULL) f1->setFirstTet(t); else f1->setSecondTet(t);
		if (f2->t1() == NULL) f2->setFirstTet(t); else f2->setSecondTet(t);
		if (f3->t1() == NULL) f3->setFirstTet(t); else f3->setSecondTet(t);
		if (f4->t1() == NULL) f4->setFirstTet(t); else f4->setSecondTet(t);
		
		assert(isTetPositive(t));
		return t;
	}

	Tetrahedron* CreateTet(TetVertex* v1, TetFace *f1) {
		TetVertex* v2 = f1->v0(), * v3 = f1->v1(), * v4 = f1->v2();

		//if (vOrient3D(v1, v2, v3, v4) > 0) return createTet(v1, v2, v3, v4);
		//else return createTet(v1, v2, v4, v3);

		if (f1->t1() == NULL) 
			return createTet(v1, v2, v3, v4);
		else 
			return createTet(v1, v2, v4, v3);
	}

	public:

	int checkConnectivity(bool allow_unlinked_elements =false)
	{
		TetEdge* c;
		TetEdges ve;
		TetFaces ef;
		unsigned int j;
		bool nm;

		for (TetVertex *v : V)
		{
			if (v == NULL) assert(0 && "checkConnectivity: detected NULL element in V list!\n");
			if (allow_unlinked_elements && !v->isLinked()) continue;
			if (v->getIncidentEdge() == NULL) assert(0 && "checkConnectivity: detected NULL e0 pointer for a vertex!\n");
			if (!v->getIncidentEdge()->hasVertex(v)) 
				assert(0 && "checkConnectivity: detected wrong e0 pointer for a vertex!\n");
		}

		for (TetEdge *e : E)
		{
			if (e == NULL) assert(0 && "checkConnectivity: detected NULL element in E list!\n");
			if (allow_unlinked_elements && !e->isLinked()) continue;
			if (e->v1() == NULL || e->v0() == NULL) assert(0 && "checkConnectivity: detected edge with one or two NULL end-points!\n");
			if (e->v1() == e->v0()) assert(0 && "checkConnectivity: detected edge with two coincident end-points!\n");
			if (e->getIncidentFace() == NULL) assert(0 && "checkConnectivity: detected NULL f0 pointer for an edge!\n");
			if (!e->getIncidentFace()->hasEdge(e)) assert(0 && "checkConnectivity: detected wrong f0 pointer for an edge!\n");
		}

		for (TetFace *f : F)
		{
			if (f == NULL) assert(0 && "checkConnectivity: detected NULL element in F list!\n");
			if (allow_unlinked_elements && !f->isLinked()) continue;
			if (f->e0() == NULL || f->e1() == NULL || f->e2() == NULL)
				assert(0 && "checkConnectivity: detected face with NULL bounding edges!\n");
			if (f->e0() == f->e1() || f->e0() == f->e2() || f->e1() == f->e2())
				assert(0 && "checkConnectivity: detected face with coincident edges!\n");
			if (f->t1() == NULL && f->t2() == NULL) assert(0 && "checkConnectivity: detected face with no incident tets!\n");

			if (f->t1() != NULL)
			{
				if (!f->t1()->hasFace(f)) assert(0 && "checkConnectivity: detected wrong t1 tet at a face!\n");
				if (!f->t1()->sameOrientation(f)) assert(0 && "checkConnectivity: wrong orientation at t1!\n");
			}
			if (f->t2() != NULL)
			{
				if (!f->t2()->hasFace(f)) assert(0 && "checkConnectivity: detected wrong t2 tet at a face!\n");
				if (f->t2()->sameOrientation(f)) assert(0 && "checkConnectivity: wrong orientation at t2!\n");
			}
		}

		for (Tetrahedron *t : T)
		{
			if (t == NULL) assert(0 && "checkConnectivity: detected NULL element in T list!\n");
			if (allow_unlinked_elements && !t->isLinked()) continue;
			if (t->f1() == NULL || t->f2() == NULL || t->f3() == NULL || t->f0() == NULL)
				assert(0 && "checkConnectivity: detected NULL as a tet face!\n");
			if (t->f1() == t->f2() || t->f1() == t->f3() || t->f1() == t->f0() ||
				t->f2() == t->f3() || t->f2() == t->f0() || t->f3() == t->f0()
				) assert(0 && "checkConnectivity: detected coincident faces for a tet!\n");
			if (t->e1() == NULL || t->e2() == NULL || t->e3() == NULL ||
				t->e4() == NULL || t->e5() == NULL || t->e0() == NULL
				) assert(0 && "checkConnectivity: tet faces do not share edges!\n");
			if (t->f1()->t1() != t && t->f1()->t2() != t)
				assert(0 && "checkConnectivity: detected tet with 1st face not pointing to the tet itself!\n");
			if (t->f2()->t1() != t && t->f2()->t2() != t)
				assert(0 && "checkConnectivity: detected tet with 2nd face not pointing to the tet itself!\n");
			if (t->f3()->t1() != t && t->f3()->t2() != t)
				assert(0 && "checkConnectivity: detected tet with 3rd face not pointing to the tet itself!\n");
			if (t->f0()->t1() != t && t->f0()->t2() != t)
				assert(0 && "checkConnectivity: detected tet with 4th face not pointing to the tet itself!\n");
		}


		for (TetFace *f : F) if (f->isLinked())
		{
			for (c = f->e0(); c != NULL; c = (c == f->e2()) ? (NULL) : (f->nextEdge(c)))
			{
				c->EF(ef);
				for (nm = true, j = 0; j < ef.size(); j++)
				{
					if (ef[j] != f && ef[j]->hasEdges(f->nextEdge(c), f->prevEdge(c)))
						assert(0 && "checkConnectivity: detected duplicate face!\n");
					if (ef[j] == f) nm = false;
				}
				ef.clear();
				if (nm) assert(0 && "checkConnectivity: detected non manifold edge (missing face)!\n");
			}
		}

		for (TetEdge *e : E) if (e->isLinked())
		{
			e->v1()->VE(ve);
			for (nm = true, j = 0; j < ve.size(); j++)
			{
				if (ve[j] != e && ve[j]->oppositeVertex(e->v1()) == e->v0()) assert(0 && "checkConnectivity: detected duplicate edge!\n");
				if (ve[j] == e) nm = false;
			}
			ve.clear();
			if (nm) assert(0 && "checkConnectivity: detected non manifold vertex!\n");
			e->v0()->VE(ve);
			for (nm = true, j = 0; j < ve.size(); j++)
			{
				if (ve[j] != e && ve[j]->oppositeVertex(e->v0()) == e->v1()) assert(0 && "checkConnectivity: detected duplicate edge!\n");
				if (ve[j] == e) nm = false;
			}
			ve.clear();
			if (nm) assert(0 && "checkConnectivity: detected non manifold vertex!\n");
		}

		printf("checkConnectivity Passed.\n");

		return 0;
	}

	void checkDelaunayCondition() {
		Tetrahedron* t;
		const TetVertex* v;
		for (TetFace* f : F) if (f->isLinked() && !f->isOnBoundary()) {
			t = f->t1();
			v = f->t2()->oppositeVertex(f);
			if (inTetCircumsphere((TetVertex*)v, t)) assert(0 && "t1 circumsphere contains opposite vertex!\n");
			t = f->t2();
			v = f->t1()->oppositeVertex(f);
			if (inTetCircumsphere((TetVertex*)v, t)) assert(0 && "t2 circumsphere contains opposite vertex!\n");
		}
		printf("checkDelaunayCondition Passed.\n");
	}


	void savePLCFaces(const char* filename) {
		size_t num_v = 0, num_t = 0;

		for (TetVertex* v : V) v->setInfo(0);

		for (auto g = G.begin(); g < G.end() - 6; g++)
			for (DelTriangle* f : (*g)->getTriangles()) if (f->is_internal) {
				f->v0()->setInfo((void*)1);
				f->v1()->setInfo((void*)1);
				f->v2()->setInfo((void*)1);
				num_t++;
			}

		for (TetVertex* v : V) if (v->getInfo()) {
			num_v++;
			v->setInfo((void*)(num_v));
		}

		FILE* fp;
		if ((fp = fopen(filename, "w")) == NULL) return;

		fprintf(fp, "OFF\n");
		fprintf(fp, "%zu %zu 0\n", num_v, num_t);

		double x, y, z;
		for (TetVertex* v : V) if (v->getInfo()) {
			const pointType* p = v->getPoint();
			p->getApproxXYZCoordinates(x, y, z);
			fprintf(fp, "%f %f %f\n", x, y, z);
		}

		size_t i1, i2, i3;
		for (auto g = G.begin(); g < G.end() - 6; g++)
			for (DelTriangle* f : (*g)->getTriangles()) if (f->is_internal) {
				i1 = (size_t)f->v0()->getInfo();
				i2 = (size_t)f->v1()->getInfo();
				i3 = (size_t)f->v2()->getInfo();
				fprintf(fp, "3 %zu %zu %zu\n", i1 - 1, i2 - 1, i3 - 1);
			}

		for (TetVertex* v : V) v->setInfo(0);

		fclose(fp);
	}

	void checkAllFaces() {
		for (PLC_Face* f : G) f->check();
	}

	// Added by Lorenzo 12/12/2024
	void export_DelTris_asTriVrtsInds(std::vector<uint32_t>& del_tri, bool exclude_bbtris) const {
		del_tri.reserve(F.size());
		for(const TetFace* tri : F) if(tri->deltri != NULL){
			if(exclude_bbtris && (tri->t1()==NULL || tri->t2()==NULL) ) continue; // do not use bounding box boundary faces
			const DelTriangle* t = tri->deltri;
			del_tri.insert(del_tri.end(), 
					{ (uint32_t) t->v0()->getIndex(), 
					  (uint32_t) t->v1()->getIndex(), 
					  (uint32_t) t->v2()->getIndex() } );
		}
	}
	size_t count_DelTris(bool exclude_bbtris) const {
		size_t counter = 0;
		for(const TetFace* tri : F) if(tri->deltri != NULL){
			if(exclude_bbtris && (tri->t1()==NULL || tri->t2()==NULL) ) continue; // do not count bounding box boundary faces
			counter++;
		}
		return counter;
	}
};


void TetVertex::VT(Tetrahedra& tets) const
{
	assert(tets.empty());
	const TetFace* f0 = e0->getIncidentFace();
	Tetrahedron* n, *t = f0->t1();
	assert(t != NULL);

	tets.push_back(t); t->mark<7>();
	for (size_t i = 0; i < tets.size(); i++) {
		t = tets[i];
		n = t->t0(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(this)) { tets.push_back(n); n->mark<7>(); }
		n = t->t1(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(this)) { tets.push_back(n); n->mark<7>(); }
		n = t->t2(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(this)) { tets.push_back(n); n->mark<7>(); }
		n = t->t3(); if (n != NULL && !n->isMarked<7>() && n->hasVertex(this)) { tets.push_back(n); n->mark<7>(); }
	}

	for (Tetrahedron* t : tets) t->unmark<7>();
}

void TetVertex::VF(TetFaces& faces) const
{
	assert(faces.empty());
	TETMESH_STATIC Tetrahedra tets;
	VT(tets);
	TetFace* f;
	for (Tetrahedron* t : tets) {
		f = t->f1(); if (!f->isMarked<7>() && f->hasVertex(this)) { faces.push_back(f); f->mark<7>(); }
		f = t->f2(); if (!f->isMarked<7>() && f->hasVertex(this)) { faces.push_back(f); f->mark<7>(); }
		f = t->f3(); if (!f->isMarked<7>() && f->hasVertex(this)) { faces.push_back(f); f->mark<7>(); }
		f = t->f0(); if (!f->isMarked<7>() && f->hasVertex(this)) { faces.push_back(f); f->mark<7>(); }
	}
	tets.clear();

	for (TetFace* f : faces) f->unmark<7>();
}

void TetVertex::VE(TetEdges& edges) const
{
	assert(edges.empty());
	TETMESH_STATIC TetFaces faces;
	VF(faces);
	TetEdge* e;
	for (TetFace* f : faces) {
		e = f->e0(); if (!e->isMarked<7>() && e->hasVertex((TetVertex*)this)) { edges.push_back(e); e->mark<7>(); }
		e = f->e1(); if (!e->isMarked<7>() && e->hasVertex((TetVertex*)this)) { edges.push_back(e); e->mark<7>(); }
		e = f->e2(); if (!e->isMarked<7>() && e->hasVertex((TetVertex*)this)) { edges.push_back(e); e->mark<7>(); }
	}
	faces.clear();

	for (TetEdge* e : edges) e->unmark<7>();
}

void TetVertex::VV(TetVertices& vs) const
{
	assert(vs.empty());
	TETMESH_STATIC TetEdges edges;
	VE(edges);
	for (TetEdge* e : edges) vs.push_back((TetVertex *)e->oppositeVertex(this));
	edges.clear();
}


TetEdge* TetVertex::getEdge(const TetVertex* vo) const {
	TETMESH_STATIC TetEdges edges;
	VE(edges);
	TetEdge* ret = NULL;
	for (TetEdge* e : edges) if (e->hasVertex((TetVertex*)vo)) { ret = e; break; }

	edges.clear();
	return ret;
}

void TetEdge::ET(Tetrahedra& lt) const {
	assert(lt.empty());

	TetFace* f = f0;
	Tetrahedron* t = f->t1();
	while (t != NULL)
	{
		lt.push_back(t);
		f = t->oppositeEdgeFace(this, f);
		if (f == f0) return;
		t = f->oppositeTet(t);
	}
	std::reverse(lt.begin(), lt.end());

	f = f0;
	t = f->t2();
	while (t != NULL)
	{
		lt.push_back(t);
		f = t->oppositeEdgeFace(this, f);
		t = f->oppositeTet(t);
	}
}

void TetEdge::EF(TetFaces& lf) const {
	assert(lf.empty());

	Tetrahedron* t;
	TetFace* f;

	lf.push_back(f0);

	f = f0;
	if ((t = f->t1()) != NULL)
	{
		f = t->oppositeEdgeFace(this, f);
		while (f != f0)
		{
			lf.push_back(f);
			if ((t = f->oppositeTet(t)) == NULL) break;
			f = t->oppositeEdgeFace(this, f);
		}
	}

	if (f == f0) return;

	for (unsigned int i = 0; i < (lf.size() / 2); i++)
	{
		f = lf[lf.size() - 1 - i]; lf[lf.size() - 1 - i] = lf[i]; lf[i] = f;
	}

	f = f0;
	if ((t = f->t2()) != NULL)
	{
		f = t->oppositeEdgeFace(this, f);
		while (f != f0)
		{
			lf.push_back(f);
			if ((t = f->oppositeTet(t)) == NULL) break;
			f = t->oppositeEdgeFace(this, f);
		}
	}
}

TetFace* TetEdge::getFace(TetVertex* vo) const {
	TETMESH_STATIC TetFaces faces;
	EF(faces);
	TetFace* ret = NULL;
	for (TetFace* f : faces) if (f->hasVertex(vo)) { ret = f; break; }

	faces.clear();
	return ret;
}

bool TetFace::isInterface() const {
	return (it2 == NULL && it1->is_internal) || (it2 != NULL && it1->is_internal != it2->is_internal);
}

bool PLC_Segment::isEncroached() const {
	TETMESH_STATIC TetFaces faces;
	TetEdge* e = Tetrahedrization::getEdge((TetVertex *)v0(), (TetVertex*)v1()); // v0()->getEdge(v1());
	if (e == NULL) return true; // If missing -> encroached

	e->EF(faces);
	for (TetFace* f : faces) {
		const pointType* pq = f->oppositeVertex(e)->getPoint();
		if (isEncroachedBy(pq)) { faces.clear(); return true; }
	}

	faces.clear();
	return false;
}

bool TetEdge::isEncroached() const {
	TETMESH_STATIC TetFaces faces;
	EF(faces);
	for (TetFace* f : faces) {
		const pointType* pq = f->oppositeVertex((TetEdge *)this)->getPoint();
		if (isEncroachedBy(pq)) { faces.clear(); return true; }
	}

	faces.clear();
	return false;
}

bool TetFace::isEncroached() const {
	const pointType* v0p = v0()->getPoint();
	const pointType* v1p = v1()->getPoint();
	const pointType* v2p = v2()->getPoint();
	const TetVertex* o1 = t1()->oppositeVertex(this);
	if (pointType::inGabrielSphere(*o1->getPoint(), *v0p, *v1p, *v2p) <= 0) return true;
	if (t2() == NULL) return false;
	const TetVertex* o2 = t2()->oppositeVertex(this);
	if (pointType::inGabrielSphere(*o2->getPoint(), *v0p, *v1p, *v2p) <= 0) return true;
	return false;
}

bool TetFace::isWeaklyEncroached() const {
	const pointType* v0p = v0()->getPoint();
	const pointType* v1p = v1()->getPoint();
	const pointType* v2p = v2()->getPoint();
	const TetVertex* o1 = t1()->oppositeVertex(this);
	if (pointType::inGabrielSphere(*o1->getPoint(), *v0p, *v1p, *v2p) < 0) return true;
	if (t2() == NULL) return false;
	const TetVertex* o2 = t2()->oppositeVertex(this);
	if (pointType::inGabrielSphere(*o2->getPoint(), *v0p, *v1p, *v2p) < 0) return true;
	return false;
}





//////////////////////////////////////////////////////////////////////////
// 
// LOCAL 2D DELAUNAY TRIANGULATION FOR PLC FACES
//
//////////////////////////////////////////////////////////////////////////

int PLC_Face::orient2D(const pointType* p1, const pointType* p2, const pointType* p3) const {
	switch (max_normal_component) {
	case 3: return pointType::orient2Dxy(*p1, *p2, *p3);
	case 1: return pointType::orient2Dyz(*p1, *p2, *p3);
	case 2: return pointType::orient2Dzx(*p1, *p2, *p3);
	case -3: return pointType::orient2Dxy(*p1, *p3, *p2);
	case -1: return pointType::orient2Dyz(*p1, *p3, *p2);
	default:
		assert(max_normal_component == -2);
		return pointType::orient2Dzx(*p1, *p3, *p2);
	}
}

int PLC_Face::vOrient2D(const TetVertex* v1, const TetVertex* v2, const TetVertex* v3) const {
	return orient2D(v1->getPoint(), v2->getPoint(), v3->getPoint());
}

int PLC_Face::symbolicPerturbation(const TetVertex* indices[4]) const {
	int swaps = 0;
	int n = 4;
	int count;
	do {
		count = 0;
		n--;
		for (int i = 0; i < n; i++) {
			if (indices[i]->getIndex() > indices[i + 1]->getIndex()) {
				std::swap(indices[i], indices[i + 1]);
				count++;
			}
		}
		swaps += count;
	} while (count);

	n = vOrient2D(indices[1], indices[2], indices[3]);
	if (n) return (swaps % 2) ? (n) : -n;

	n = vOrient2D(indices[0], indices[2], indices[3]);
	return (swaps % 2) ? (-n) : (n);
}


bool PLC_Face::vInDiametralSphere(const TetVertex* v1, const TetVertex* v2, const TetVertex* v3, const TetVertex* v4) const {
	const pointType* p1 = v1->getPoint();
	const pointType* p2 = v2->getPoint();
	const pointType* p3 = v3->getPoint();
	const pointType* p4 = v4->getPoint();

	int det = pointType::inGabrielSphere(*p1, *p2, *p3, *p4);

	if (det == 0) {
		const TetVertex* idx[4] = { v1, v2, v3, v4 };
		det = symbolicPerturbation(idx);
	}

	return (det < 0);
}

DelTriangle* PLC_Face::searchTriangle(TetVertex* v) const {
	const pointType* p = v->getPoint();
	DelTriangle* t = T.back();

	do {
		assert(t != NULL);
		const pointType* v0 = t->v0()->getPoint();
		const pointType* v1 = t->v1()->getPoint();
		const pointType* v2 = t->v2()->getPoint();
		if (orient2D(p, v0, v1) < 0) t = t->e2()->oppositeTriangle(t);
		else if (orient2D(p, v1, v2) < 0) t = t->e0()->oppositeTriangle(t);
		else if (orient2D(p, v2, v0) < 0) t = t->e1()->oppositeTriangle(t);
		else return t;
	} while (t != NULL);

	return NULL;
}

void PLC_Face::splitAndDelaunizeEdge(TetVertex* v, DelEdge* e) {
	TetVertex* v0 = e->v0(), * v1 = e->v1();
	DelTriangle* t0 = e->firstTriangle(), * t1 = e->secondTriangle();

	//assert(vOrient2D(v, v0, v1) == 0);

	DelEdge* ne = new DelEdge(v, v1); E.push_back(ne);
	e->replaceVertex(v1, v);

	DelEdge* ne0 = t0->nextEdge(e), *pe0 = t0->prevEdge(e);
	TetVertex* ov0 = (TetVertex *)t0->oppositeVertex(e);
	DelEdge* no0 = new DelEdge(v, ov0); E.push_back(no0);
	DelTriangle* nt0 = new DelTriangle(ne, ne0, no0, this); T.push_back(nt0);
	nt0->is_internal = t0->is_internal;
	no0->setFirstTriangle(t0); no0->setSecondTriangle(nt0);
	ne->setFirstTriangle(nt0);
	ne0->replaceTriangle(t0, nt0);
	t0->replaceEdge(ne0, no0);

	DelEdges toswap;

	// Corner case: border edge
	if (t1 != NULL) {
		DelEdge* pe1 = t1->prevEdge(e), * ne1 = t1->nextEdge(e);
		TetVertex* ov1 = (TetVertex*)t1->oppositeVertex(e);
		DelEdge* no1 = new DelEdge(v, ov1); E.push_back(no1);
		DelTriangle* nt1 = new DelTriangle(ne, no1, pe1, this); T.push_back(nt1);
		nt1->is_internal = t1->is_internal;
		no1->setFirstTriangle(nt1); no1->setSecondTriangle(t1);
		ne->setSecondTriangle(nt1);
		pe1->replaceTriangle(t1, nt1);
		t1->replaceEdge(pe1, no1);
		toswap = { ne0, pe0, ne1, pe1, no0, no1 };
	}
	else toswap = { ne0, pe0 };

	iterativeSwap(toswap);
}

bool PLC_Face::nearlyCollinear(const pointType* ap, const pointType* bp, const pointType* cp) {
	const vector3d a(ap), b(bp), c(cp);

	vector3d ab = a - b;
	vector3d cb = c - b;
	ab *= (1.0 / sqrt(ab.sq_length()));
	cb *= (1.0 / sqrt(cb.sq_length()));
	double cross = (ab & cb).sq_length();
	return cross < 1.0e-20;
}

void PLC_Face::splitAndDelaunizeTriangle(TetVertex* v, DelTriangle* t) {
	TetVertex* v0 = t->v0(), * v1 = t->v1(), * v2 = t->v2();
	DelEdge* e0 = t->e0(), * e1 = t->e1(), * e2 = t->e2();

	const pointType* vp = v->getPoint();
	const pointType* v0p = v0->getPoint();
	const pointType* v1p = v1->getPoint();
	const pointType* v2p = v2->getPoint();

	// Corner case: vertex is on one of the edges -> swap
	if (orient2D(v0p, vp, v1p) == 0) splitAndDelaunizeEdge(v, e2);
	else if (orient2D(v1p, vp, v2p) == 0) splitAndDelaunizeEdge(v, e0);
	else if (orient2D(v2p, vp, v0p) == 0) splitAndDelaunizeEdge(v, e1);
	else {
		// Split triangle
		DelEdge* ne0 = new DelEdge(v, v0); E.push_back(ne0);
		DelEdge* ne1 = new DelEdge(v, v1); E.push_back(ne1);
		DelEdge* ne2 = new DelEdge(v, v2); E.push_back(ne2);
		DelTriangle* nt1 = new DelTriangle(ne0, ne2, e1, this); T.push_back(nt1);
		DelTriangle* nt2 = new DelTriangle(ne1, ne0, e2, this); T.push_back(nt2);
		nt1->is_internal = nt2->is_internal = t->is_internal;
		ne0->setFirstTriangle(nt2); ne0->setSecondTriangle(nt1);
		ne1->setFirstTriangle(t); ne1->setSecondTriangle(nt2);
		ne2->setFirstTriangle(nt1); ne2->setSecondTriangle(t);
		e1->replaceTriangle(t, nt1);
		e2->replaceTriangle(t, nt2);
		t->replaceEdge(e1, ne2);
		t->replaceEdge(e2, ne1);

		DelEdges toswap = { e0, e1, e2 };
		iterativeSwap(toswap);
	}
}

void PLC_Face::iterativeSwap(DelEdges& toswap) {
	DelEdge* diamond[4];
	while (!toswap.empty()) {
		auto t = toswap.back(); toswap.pop_back();
		if (swapIfNotDelaunay(t, diamond)) {
			for (int i = 0; i < 4; i++)
				toswap.push_back(diamond[i]);
		}
	}
}

bool PLC_Face::swapIfNotDelaunay(DelEdge* e, DelEdge* diamond[4]) {
	const DelTriangle* et1 = e->firstTriangle(), * et2 = e->secondTriangle();
	const TetVertex* vb = et1->oppositeVertex(e);
	if (et2 == NULL || (et1->is_internal != et2->is_internal)) return false;
	const TetVertex* vc = e->v0(), * va = e->v1();
	const TetVertex* vq = et2->oppositeVertex(e);

	if (vInDiametralSphere(vq, va, vb, vc)) {
		e->swap();
		diamond[0] = et1->nextEdge(e);
		diamond[1] = et1->prevEdge(e);
		diamond[2] = et2->nextEdge(e);
		diamond[3] = et2->prevEdge(e);
		return true;
	}

	return false;
}

void PLC_Face::addExternalPoint(TetVertex* v) {
	DelEdges toswap;
	DelEdge* prev_bd=NULL, * next_bd;
	DelEdge* prev_ne, * next_ne;

	// Find one boundary edge visible from v
	for (DelEdge* e : E) 
		if (e->isOnBoundary() && vOrient2D(v, e->v1(), e->v0()) > 0) {
			prev_bd = e->prevOnBoundary();
			next_bd = e->nextOnBoundary();
			next_ne = new DelEdge(v, e->v1()); E.push_back(next_ne);
			prev_ne = new DelEdge(e->v0(), v); E.push_back(prev_ne);
			DelTriangle* t = new DelTriangle(next_ne, e, prev_ne, this); T.push_back(t);
			t->is_internal = false;
			next_ne->setFirstTriangle(t);
			prev_ne->setFirstTriangle(t);
			e->setSecondTriangle(t);
			toswap.push_back(e);
			break;
		}
	
	// Move along bondary and keep making triangles while visibility holds
	while (vOrient2D(v, prev_bd->v1(), prev_bd->v0()) > 0) {
		DelEdge* e = prev_bd;
		prev_bd = e->prevOnBoundary();
		DelEdge* nne = E.back();
		DelEdge* pne = new DelEdge(e->v0(), v); E.push_back(pne);
		DelTriangle* t = new DelTriangle(nne, e, pne, this); T.push_back(t);
		t->is_internal = false;
		nne->setSecondTriangle(t);
		pne->setFirstTriangle(t);
		e->setSecondTriangle(t);
		toswap.push_back(e);
	}

	// Same as above but in opposite direction
	DelEdge* pne = next_ne;
	while (vOrient2D(v, next_bd->v1(), next_bd->v0()) > 0) {
		DelEdge* e = next_bd;
		next_bd = e->nextOnBoundary();
		DelEdge* nne = new DelEdge(v, e->v1()); E.push_back(nne);
		DelTriangle* t = new DelTriangle(nne, e, pne, this); T.push_back(t);
		t->is_internal = false;
		pne->setSecondTriangle(t);
		nne->setFirstTriangle(t);
		e->setSecondTriangle(t);
		toswap.push_back(e);
		pne = nne;
	}

	iterativeSwap(toswap);
}


void PLC_Face::insertExistingVertex(TetVertex *v) {
	DelTriangle* t = searchTriangle(v);
	if (t != NULL) splitAndDelaunizeTriangle(v, t);
	else addExternalPoint(v);
}

void PLC_Face::insertExistingVertexOnSegment(TetVertex* v, TetVertex* sv0, TetVertex* sv1) {
	for (DelEdge* e : E) if (e->hasVertex(sv0) && e->hasVertex(sv1)) {
		splitAndDelaunizeEdge(v, e);
		return;
	}
	//ip_error("PLC_Face::insertExistingVertexOnSegment: Could not find edge\n");
	std::cout<<"[tetmesh.h - insertExistingVertexOnSegment()]  Could not find edge\n"; exit(1);
}

void DelEdge::swap() {
	TetVertex* ev0 = v0(), *ev1 = v1();
	assert(t2 != NULL);
	TetVertex* ov0 = (TetVertex*)t1->oppositeVertex(this);
	TetVertex* ov1 = (TetVertex*)t2->oppositeVertex(this);
	DelEdge* ne0 = t1->nextEdge(this), * pe0 = t1->prevEdge(this);
	DelEdge* pe1 = t2->prevEdge(this), * ne1 = t2->nextEdge(this);
	_v[0] = ov1;
	_v[1] = ov0;
	DelEdge** t1e = t1->_e, ** t2e = t2->_e;
	t1e[0] = pe0; t1e[1] = ne1; t1e[2] = this;
	t2e[0] = pe1; t2e[1] = ne0; t2e[2] = this;
	ne0->replaceTriangle(t1, t2);
	ne1->replaceTriangle(t2, t1);
}

void PLC_Face::delaunizeInitialTriangles() {
	// Take three vertices to start from
	DelTriangle* t = T.back();
	TetVertex* v0 = t->v0(), * v1 = t->v1(), * v2 = t->v2();
	v0->mark<7>();
	v1->mark<7>();
	v2->mark<7>();

	// Collect all other vertices
	TetVertex* v;
	TetVertices vs;
	for (PLC_Segment* s : segments) {
		v = s->v0(); if (!v->isMarked<7>()) { v->mark<7>(); vs.push_back(v); }
		v = s->v1(); if (!v->isMarked<7>()) { v->mark<7>(); vs.push_back(v); }
	}
	for (DelTriangle* t : T) {
		v = t->v0(); if (!v->isMarked<7>()) { v->mark<7>(); vs.push_back(v); }
		v = t->v1(); if (!v->isMarked<7>()) { v->mark<7>(); vs.push_back(v); }
		v = t->v2(); if (!v->isMarked<7>()) { v->mark<7>(); vs.push_back(v); }
	}
	for (TetVertex* v : vs) v->unmark<7>();
	v0->unmark<7>();
	v1->unmark<7>();
	v2->unmark<7>();

	// Destroy current triangles
	destroyUnusedElements<DelEdge*>(E, [](DelEdge* x) { return true; });
	destroyUnusedElements<DelTriangle*>(T, [](DelTriangle* x) { return true; });

	DelEdge* e0 = new DelEdge(v0, v1); E.push_back(e0);
	DelEdge* e1 = new DelEdge(v1, v2); E.push_back(e1);
	DelEdge* e2 = new DelEdge(v2, v0); E.push_back(e2);
	DelTriangle* t0 = new DelTriangle(e0, e1, e2, this); T.push_back(t0);
	e0->setFirstTriangle(t0);
	e1->setFirstTriangle(t0);
	e2->setFirstTriangle(t0);
	t0->is_internal = false;

	assert(vOrient2D(v0, v1, v2) > 0);

	for (TetVertex* v : vs)
		insertExistingVertex(v);
}

void PLC_Face::check() {
	for (DelEdge* e : E)
		if (!e->firstTriangle()->hasEdge(e)) 
			assert(0 && "missing edge in incident tri\n");
		else if (e->secondTriangle() != NULL && !e->secondTriangle()->hasEdge(e))
			assert(0 && "missing edge in incident tri 2\n");
	for (DelTriangle* t : T)
		if (!t->e0()->commonVertex(t->e1()) || !t->e1()->commonVertex(t->e2()) || !t->e2()->commonVertex(t->e0())) 
			assert(0 && "Wrong edges\n");
	for (DelTriangle* t : T)
		if (!t->e0()->hasTriangle(t) || !t->e1()->hasTriangle(t) || !t->e2()->hasTriangle(t))
			assert(0 && "Wrong edges\n");
	for (DelTriangle* t : T)
			assert(vOrient2D(t->v0(), t->v1(), t->v2()) > 0);
	for (DelEdge* e : E) {
		assert(vOrient2D(e->v0(), e->v1(), e->firstTriangle()->oppositeVertex(e)) > 0);
		if (e->secondTriangle() != NULL)
			assert(vOrient2D(e->v1(), e->v0(), e->secondTriangle()->oppositeVertex(e)) > 0);
	}
}




