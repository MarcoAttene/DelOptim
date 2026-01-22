#define _USE_MATH_DEFINES
#include <math.h>

// --------------------- //
// DIHEDRAL ACUTE EXACT  //
// --------------------- //
// This section should be moved in another source file, it is EXATCT

inline int dihedralAcuteness_filtered(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz)
{
	const double ux = cx - dx;
	const double uy = cy - dy;
	const double uz = cz - dz;
	const double vx = bx - dx;
	const double vy = by - dy;
	const double vz = bz - dz;
	const double wx = ax - dx;
	const double wy = ay - dy;
	const double wz = az - dz;
	const double uxsq = ux * ux;
	const double uysq = uy * uy;
	const double uzsq = uz * uz;
	const double qx = uysq + uzsq;
	const double qy = uxsq + uzsq;
	const double qz = uxsq + uysq;
	const double tx = vx * wx;
	const double ty = vy * wy;
	const double tz = vz * wz;
	const double mx = qx * tx;
	const double my = qy * ty;
	const double mz = qz * tz;
	const double dl1 = mx + my;
	const double dl = dl1 + mz;
	const double rx = uy * uz;
	const double ry = ux * uz;
	const double rz = ux * uy;
	const double sx1 = vy * wz;
	const double sx2 = vz * wy;
	const double sx = sx1 + sx2;
	const double sy1 = vx * wz;
	const double sy2 = vz * wx;
	const double sy = sy1 + sy2;
	const double sz1 = vx * wy;
	const double sz2 = vy * wx;
	const double sz = sz1 + sz2;
	const double nx = rx * sx;
	const double ny = ry * sy;
	const double nz = rz * sz;
	const double dr1 = nx + ny;
	const double dr = dr1 + nz;
	const double d = dl - dr;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(ux)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(uy)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(uz)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(vx)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(vy)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(vz)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(wx)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(wy)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(wz)) > max_var) max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= epsilon;
	epsilon *= 1.332267629550189e-14;
	if (d > epsilon) return IP_Sign::POSITIVE;
	if (-d > epsilon) return IP_Sign::NEGATIVE;
	return Filtered_Sign::UNCERTAIN;
}

inline int dihedralAcuteness_interval(interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz, interval_number dx, interval_number dy, interval_number dz)
{
	setFPUModeToRoundUP();
	const interval_number ux(cx - dx);
	const interval_number uy(cy - dy);
	const interval_number uz(cz - dz);
	const interval_number vx(bx - dx);
	const interval_number vy(by - dy);
	const interval_number vz(bz - dz);
	const interval_number wx(ax - dx);
	const interval_number wy(ay - dy);
	const interval_number wz(az - dz);
	const interval_number uxsq(ux * ux);
	const interval_number uysq(uy * uy);
	const interval_number uzsq(uz * uz);
	const interval_number qx(uysq + uzsq);
	const interval_number qy(uxsq + uzsq);
	const interval_number qz(uxsq + uysq);
	const interval_number tx(vx * wx);
	const interval_number ty(vy * wy);
	const interval_number tz(vz * wz);
	const interval_number mx(qx * tx);
	const interval_number my(qy * ty);
	const interval_number mz(qz * tz);
	const interval_number dl1(mx + my);
	const interval_number dl(dl1 + mz);
	const interval_number rx(uy * uz);
	const interval_number ry(ux * uz);
	const interval_number rz(ux * uy);
	const interval_number sx1(vy * wz);
	const interval_number sx2(vz * wy);
	const interval_number sx(sx1 + sx2);
	const interval_number sy1(vx * wz);
	const interval_number sy2(vz * wx);
	const interval_number sy(sy1 + sy2);
	const interval_number sz1(vx * wy);
	const interval_number sz2(vy * wx);
	const interval_number sz(sz1 + sz2);
	const interval_number nx(rx * sx);
	const interval_number ny(ry * sy);
	const interval_number nz(rz * sz);
	const interval_number dr1(nx + ny);
	const interval_number dr(dr1 + nz);
	const interval_number d(dl - dr);
	setFPUModeToRoundNEAR();

	if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
	return d.sign();
}

inline int dihedralAcuteness_bigfloat(bigfloat ax, bigfloat ay, bigfloat az, bigfloat bx, bigfloat by, bigfloat bz, bigfloat cx, bigfloat cy, bigfloat cz, bigfloat dx, bigfloat dy, bigfloat dz)
{
	const bigfloat ux(cx - dx);
	const bigfloat uy(cy - dy);
	const bigfloat uz(cz - dz);
	const bigfloat vx(bx - dx);
	const bigfloat vy(by - dy);
	const bigfloat vz(bz - dz);
	const bigfloat wx(ax - dx);
	const bigfloat wy(ay - dy);
	const bigfloat wz(az - dz);
	const bigfloat uxsq(ux * ux);
	const bigfloat uysq(uy * uy);
	const bigfloat uzsq(uz * uz);
	const bigfloat qx(uysq + uzsq);
	const bigfloat qy(uxsq + uzsq);
	const bigfloat qz(uxsq + uysq);
	const bigfloat tx(vx * wx);
	const bigfloat ty(vy * wy);
	const bigfloat tz(vz * wz);
	const bigfloat mx(qx * tx);
	const bigfloat my(qy * ty);
	const bigfloat mz(qz * tz);
	const bigfloat dl1(mx + my);
	const bigfloat dl(dl1 + mz);
	const bigfloat rx(uy * uz);
	const bigfloat ry(ux * uz);
	const bigfloat rz(ux * uy);
	const bigfloat sx1(vy * wz);
	const bigfloat sx2(vz * wy);
	const bigfloat sx(sx1 + sx2);
	const bigfloat sy1(vx * wz);
	const bigfloat sy2(vz * wx);
	const bigfloat sy(sy1 + sy2);
	const bigfloat sz1(vx * wy);
	const bigfloat sz2(vy * wx);
	const bigfloat sz(sz1 + sz2);
	const bigfloat nx(rx * sx);
	const bigfloat ny(ry * sy);
	const bigfloat nz(rz * sz);
	const bigfloat dr1(nx + ny);
	const bigfloat dr(dr1 + nz);
	const bigfloat d(dl - dr);
	return sgn(d);
}

inline int dihedralAcuteness_exact(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz)
{
	double ux[2];
	expansionObject::two_Diff(cx, dx, ux);
	double uy[2];
	expansionObject::two_Diff(cy, dy, uy);
	double uz[2];
	expansionObject::two_Diff(cz, dz, uz);
	double vx[2];
	expansionObject::two_Diff(bx, dx, vx);
	double vy[2];
	expansionObject::two_Diff(by, dy, vy);
	double vz[2];
	expansionObject::two_Diff(bz, dz, vz);
	double wx[2];
	expansionObject::two_Diff(ax, dx, wx);
	double wy[2];
	expansionObject::two_Diff(ay, dy, wy);
	double wz[2];
	expansionObject::two_Diff(az, dz, wz);
	double uxsq[8];
	int uxsq_len = expansionObject::Gen_Product(2, ux, 2, ux, uxsq);
	double uysq[8];
	int uysq_len = expansionObject::Gen_Product(2, uy, 2, uy, uysq);
	double uzsq[8];
	int uzsq_len = expansionObject::Gen_Product(2, uz, 2, uz, uzsq);
	double qx[16];
	int qx_len = expansionObject::Gen_Sum(uysq_len, uysq, uzsq_len, uzsq, qx);
	double qy[16];
	int qy_len = expansionObject::Gen_Sum(uxsq_len, uxsq, uzsq_len, uzsq, qy);
	double qz[16];
	int qz_len = expansionObject::Gen_Sum(uxsq_len, uxsq, uysq_len, uysq, qz);
	double tx[8];
	int tx_len = expansionObject::Gen_Product(2, vx, 2, wx, tx);
	double ty[8];
	int ty_len = expansionObject::Gen_Product(2, vy, 2, wy, ty);
	double tz[8];
	int tz_len = expansionObject::Gen_Product(2, vz, 2, wz, tz);
	double mx_p[128], * mx = mx_p;
	int mx_len = expansionObject::Gen_Product_With_PreAlloc(qx_len, qx, tx_len, tx, &mx, 128);
	double my_p[128], * my = my_p;
	int my_len = expansionObject::Gen_Product_With_PreAlloc(qy_len, qy, ty_len, ty, &my, 128);
	double mz_p[128], * mz = mz_p;
	int mz_len = expansionObject::Gen_Product_With_PreAlloc(qz_len, qz, tz_len, tz, &mz, 128);
	double dl1_p[128], * dl1 = dl1_p;
	int dl1_len = expansionObject::Gen_Sum_With_PreAlloc(mx_len, mx, my_len, my, &dl1, 128);
	double dl_p[128], * dl = dl_p;
	int dl_len = expansionObject::Gen_Sum_With_PreAlloc(dl1_len, dl1, mz_len, mz, &dl, 128);
	double rx[8];
	int rx_len = expansionObject::Gen_Product(2, uy, 2, uz, rx);
	double ry[8];
	int ry_len = expansionObject::Gen_Product(2, ux, 2, uz, ry);
	double rz[8];
	int rz_len = expansionObject::Gen_Product(2, ux, 2, uy, rz);
	double sx1[8];
	int sx1_len = expansionObject::Gen_Product(2, vy, 2, wz, sx1);
	double sx2[8];
	int sx2_len = expansionObject::Gen_Product(2, vz, 2, wy, sx2);
	double sx[16];
	int sx_len = expansionObject::Gen_Sum(sx1_len, sx1, sx2_len, sx2, sx);
	double sy1[8];
	int sy1_len = expansionObject::Gen_Product(2, vx, 2, wz, sy1);
	double sy2[8];
	int sy2_len = expansionObject::Gen_Product(2, vz, 2, wx, sy2);
	double sy[16];
	int sy_len = expansionObject::Gen_Sum(sy1_len, sy1, sy2_len, sy2, sy);
	double sz1[8];
	int sz1_len = expansionObject::Gen_Product(2, vx, 2, wy, sz1);
	double sz2[8];
	int sz2_len = expansionObject::Gen_Product(2, vy, 2, wx, sz2);
	double sz[16];
	int sz_len = expansionObject::Gen_Sum(sz1_len, sz1, sz2_len, sz2, sz);
	double nx_p[128], * nx = nx_p;
	int nx_len = expansionObject::Gen_Product_With_PreAlloc(rx_len, rx, sx_len, sx, &nx, 128);
	double ny_p[128], * ny = ny_p;
	int ny_len = expansionObject::Gen_Product_With_PreAlloc(ry_len, ry, sy_len, sy, &ny, 128);
	double nz_p[128], * nz = nz_p;
	int nz_len = expansionObject::Gen_Product_With_PreAlloc(rz_len, rz, sz_len, sz, &nz, 128);
	double dr1_p[128], * dr1 = dr1_p;
	int dr1_len = expansionObject::Gen_Sum_With_PreAlloc(nx_len, nx, ny_len, ny, &dr1, 128);
	double dr_p[128], * dr = dr_p;
	int dr_len = expansionObject::Gen_Sum_With_PreAlloc(dr1_len, dr1, nz_len, nz, &dr, 128);
	double d_p[128], * d = d_p;
	int d_len = expansionObject::Gen_Diff_With_PreAlloc(dl_len, dl, dr_len, dr, &d, 128);

	double return_value = d[d_len - 1];
	if (d_p != d) FreeDoubles(d);
	if (dr_p != dr) FreeDoubles(dr);
	if (dr1_p != dr1) FreeDoubles(dr1);
	if (nz_p != nz) FreeDoubles(nz);
	if (ny_p != ny) FreeDoubles(ny);
	if (nx_p != nx) FreeDoubles(nx);
	if (dl_p != dl) FreeDoubles(dl);
	if (dl1_p != dl1) FreeDoubles(dl1);
	if (mz_p != mz) FreeDoubles(mz);
	if (my_p != my) FreeDoubles(my);
	if (mx_p != mx) FreeDoubles(mx);

	if (return_value > 0) return IP_Sign::POSITIVE;
	if (return_value < 0) return IP_Sign::NEGATIVE;
	if (return_value == 0) return IP_Sign::ZERO;
	return IP_Sign::UNDEFINED;
}

// +1 if there is an acute dihedral angle at edge c,d of the tetrahedron a,b,c,d
inline int dihedralAcuteness(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz)
{
	int ret;
	ret = dihedralAcuteness_filtered(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
	if (ret != Filtered_Sign::UNCERTAIN) return ret;
	ret = dihedralAcuteness_interval(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
	if (ret != Filtered_Sign::UNCERTAIN) return ret;
	return dihedralAcuteness_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
}


// ANGLES

// TRUE if triangles <e0,e1,u> and <e0,e1,v> ( where <e0,e1> has the same orientation in both triangles)
// forms an acute dihedral angle at <e0,e1>.
bool isAcuteDihedral_exact(const pointType* e0, const pointType* e1, const pointType* u, const pointType* v) {
	vector3d Oe0(e0); // D
	vector3d Oe1(e1); // C
	vector3d Ou(u); // B
	vector3d Ov(v); // A
	return dihedralAcuteness(Ov.c[0], Ov.c[1], Ov.c[2], Ou.c[0], Ou.c[1], Ou.c[2], Oe1.c[0], Oe1.c[1], Oe1.c[2], Oe0.c[0], Oe0.c[1], Oe0.c[2]) > 0;
}

// TRUE if the two triangles <v0,v1,o0> and <v0,v1,o1> form an acute dihedral angle
inline bool isAcuteDihedral_withNormalVects(const pointType* v0, const pointType* v1, const pointType* o0, const pointType* o1) {
	const vector3d v0v(v0), v1v(v1);
	const vector3d vav1 = (vector3d(o0) - v1v) & (v1v - v0v);
	const vector3d vbv1 = (vector3d(o1) - v1v) & (v1v - v0v);
	return (vbv1.dot(vav1) > 0.0);
}

class vec3d_bf {
public:
	bigfloat c[3]; // 3 coordinates

	inline vec3d_bf() { }
	inline vec3d_bf(const bigfloat x, const bigfloat y, const bigfloat z) { c[0] = x; c[1] = y; c[2] = z; }
	inline vec3d_bf(const pointType* p) {
		if (p->isExplicit3D()) {
			const explicitPoint& e = p->toExplicit3D();
			c[0] = e.X();
			c[1] = e.Y();
			c[2] = e.Z();
		}
		else {
			bigfloat d;
			p->getBigfloatLambda(c[0], c[1], c[2], d);
			assert(d == bigfloat(1));
		}
	}

	inline vec3d_bf operator+(const vec3d_bf& v) const { return vec3d_bf(c[0] + v.c[0], c[1] + v.c[1], c[2] + v.c[2]); }
	inline vec3d_bf operator-(const vec3d_bf& v) const { return vec3d_bf(c[0] - v.c[0], c[1] - v.c[1], c[2] - v.c[2]); }
	inline vec3d_bf operator*(const bigfloat d) const { return vec3d_bf(c[0] * d, c[1] * d, c[2] * d); }

	inline bigfloat dot(const vec3d_bf& p) const { return (c[0] * p.c[0] + c[1] * p.c[1] + c[2] * p.c[2]); }
	inline vec3d_bf cross(const vec3d_bf& p) const { return vec3d_bf(c[1] * p.c[2] - c[2] * p.c[1], c[2] * p.c[0] - c[0] * p.c[2], c[0] * p.c[1] - c[1] * p.c[0]); }
	inline bigfloat tripleProd(const vec3d_bf& v2, const vec3d_bf& v3) const {
		return ((v2.c[0] * v3.c[1] * c[2]) - (v3.c[0] * v2.c[1] * c[2])) +
			((v3.c[0] * c[1] * v2.c[2]) - (c[0] * v3.c[1] * v2.c[2])) +
			((c[0] * v2.c[1] * v3.c[2]) - (v2.c[0] * c[1] * v3.c[2]));
	}

	inline bigfloat operator*(const vec3d_bf& d) const { return dot(d); }
	inline vec3d_bf operator&(const vec3d_bf& d) const { return cross(d); }

	// Squared length
	inline bigfloat sq_length() const { return dot(*this); }

	// Squared distance
	inline bigfloat dist_sq(const vec3d_bf& v) const { return ((*this) - v).sq_length(); }
};


inline double sqEuclideanDistance(const pointType* a, const pointType* b) {
	return vector3d(a).dist_sq(b);
}

inline double euclideanDistance(const pointType* a, const pointType* b) {
	return sqrt(sqEuclideanDistance(a,b));
}

inline double squaredDistanceFromLine(const vector3d vv, const vector3d x1v, const vector3d x2v)
{
	vector3d x21((x2v)-(x1v));
	vector3d x10((x1v)-(vv));
	x10 = x21 & x10;
	return (x10 * x10) / (x21 * x21);
}

// The circumcenter is not exact. It should be guaranteed into the face, but rounding might make it slightly out.
// In such a case, it is virtually guaranteed that the circumcenter encroaches upon one of the bounding segments,
// meaning that it would not be used.

inline implicitPoint3D_BPT* convertToBPT(const vector3d& p, const explicitPoint* rv0, const explicitPoint* rv1, const explicitPoint* rv2) {
	// Translation into a BPT
	vector3d r0(rv0), r1(rv1), r2(rv2);
	vector3d r_cross = (r1 - r0).cross(r2 - r0);
	vector3d v_cross = (r1 - p).cross(r2 - p);
	vector3d u_cross = (r2 - p).cross(r0 - p);
	double ra = r_cross.sq_length();
	double v = sqrt(v_cross.sq_length() / ra);
	double u = sqrt(u_cross.sq_length() / ra);
	if (v_cross.dot(r_cross) < 0) v = -v;
	if (u_cross.dot(r_cross) < 0) u = -u;

	assert(fabs(v) < 1000000 && fabs(u) < 1000000);
	return new implicitPoint3D_BPT(*rv0, *rv1, *rv2, v, u);
}

inline double snapToPowerOfTwo(double d) {
	uint64_t n = *((uint64_t*)(&d));
	n &= (UINT64_MAX << 52);
	return *((double*)(&n));
}

inline pointType* createMidPoint(const pointType* v0, const pointType* v1)
{

	// return new implicitPoint3D_LNC(*v0, *v1, 0.5);

	// v0 and v1 can be Explicit, LNC or BPT (E, L, B)
	if (v0->getType() > v1->getType()) std::swap(v0, v1); // E < L < B

	// The following is the old implementation based on old predicates
	// which must take explicit point as input. 
	// QUESTION: some of the following contructors may optimize the process??
	// Cases:
	// EE -> easy, create LNC with t=0.5
	// EL -> E must be one endpoint of L (assert), create LNC (see code below)
	// EB -> should be impossible, assert
	// LL -> endpoints must coincide (assert), create LNC with avg t
	// LB -> endpoints of L must be a subset of vertices of B (assert). Convert L to B and create BPT with avg u and v
	// BB -> vertices must coincide (assert), create BPT with avg u and v

	// EE
	if (v0->isExplicit3D() && v1->isExplicit3D()) {
		// return new implicitPoint_LNC(v0->toExplicit3D(), v1->toExplicit3D(), 0.5);
		return new implicitPoint_LNC(*v0, *v1, 0.5);
	}

	// EL
	if (v0->isExplicit3D() && v1->isLNC()) {
		const implicitPoint3D_LNC& l = v1->toLNC();

		if(l.P().isExplicit3D() && l.Q().isExplicit3D()){
			// if (&l.P() == v0) return new implicitPoint_LNC(v0->toExplicit3D(), l.Q().toExplicit3D(), l.T() / 2);
			if (&l.P() == v0) return new implicitPoint_LNC(*v0, l.Q(), l.T() * 0.5);
			// return new implicitPoint_LNC(l.P().toExplicit3D(), v0->toExplicit3D(), (l.T() + 1.0) / 2);
			else if (&l.Q() == v0) return new implicitPoint_LNC(l.P(), *v0, (l.T() + 1.0) * 0.5);
		}

		else return new implicitPoint3D_LNC(*v0, *v1, 0.5);
	}

	// EB
	// assert(!(v0->isExplicit3D() && v1->isBPT()));

	// LL
	if (v0->isLNC() && v1->isLNC()) {

		const implicitPoint3D_LNC& l0 = v0->toLNC();
		const implicitPoint3D_LNC& l1 = v1->toLNC();

		if(l0.P().isExplicit3D() && l0.Q().isExplicit3D() &&
			l1.P().isExplicit3D() && l1.Q().isExplicit3D() ){

				if (&l0.P() == &l1.P() && &l0.Q() == &l1.Q()) // same orientations
					//return new implicitPoint_LNC(l0.P().toExplicit3D(), l0.Q().toExplicit3D(), (l0.T() + l1.T()) / 2);
					return new implicitPoint_LNC(l0.P(), l0.Q(), (l0.T() + l1.T()) * 0.5);
				else if (&l0.P() == &l1.Q() && &l0.Q() == &l1.P()) // opposite orientations
					// return new implicitPoint_LNC(l0.P().toExplicit3D(), l0.Q().toExplicit3D(), (l0.T() + 1.0 - l1.T()) / 2);
					return new implicitPoint_LNC(l0.P(), l0.Q(), (l0.T() + 1.0 - l1.T()) * 0.5);
				else {

					double u, v;
					const explicitPoint3D* A = NULL;
					const explicitPoint3D* B = NULL;
					const explicitPoint3D* C = NULL;
					if (&l0.P() == &l1.P()) {
						A = &l0.Q().toExplicit3D();  u = l0.T();
						B = &l1.Q().toExplicit3D();  v = l1.T();
						C = &l0.P().toExplicit3D();
					}
					else if (&l0.P() == &l1.Q()) {
						A = &l0.Q().toExplicit3D();  u = 1.0 - l0.T();
						B = &l1.P().toExplicit3D();  v = l1.T();
						C = &l0.P().toExplicit3D();
					}
					else if (&l0.Q() == &l1.P()) {
						A = &l0.P().toExplicit3D();  u = l0.T();
						B = &l1.Q().toExplicit3D();  v = 1.0 - l1.T();
						C = &l0.Q().toExplicit3D();
					}
					else if(&l0.Q() == &l1.Q()){
						A = &l0.P().toExplicit3D();  u = 1.0 - l0.T();
						B = &l1.P().toExplicit3D();  v = 1.0 - l1.T();
						C = &l0.Q().toExplicit3D();
					}
					// else assert(false);
					else return new implicitPoint3D_LNC(*v0, *v1, 0.5);

					return new implicitPoint3D_BPT(*A, *B, *C, u / 2, v / 2);
				}

			return new implicitPoint3D_LNC(*v0, *v1, 0.5);
		}

	}

	return new implicitPoint3D_LNC(*v0, *v1, 0.5);

	// // LB
	// if (v0->isLNC() && v1->isBPT()) {

	// 	std::cout<<std::endl<<"midpoint of   v0 LNC - v1 BPT"<<std::endl;

	// 	const implicitPoint3D_LNC& l0 = v0->toLNC();
	// 	const implicitPoint3D_BPT& b1 = v1->toBPT();
	// 	assert(&l0.P() == &b1.P() || &l0.P() == &b1.Q() || &l0.P() == &b1.R());
	// 	assert(&l0.Q() == &b1.P() || &l0.Q() == &b1.Q() || &l0.Q() == &b1.R());

	// 	// Convert LNC's t to BPT's (u,v)
	// 	double u = b1.U(), v = b1.V(), t = l0.T();
	// 	if (&l0.P() == &b1.Q() && &l0.Q() == &b1.P()) { u += (1 - t); v += t; }
	// 	else if (&l0.P() == &b1.P() && &l0.Q() == &b1.Q()) { u += t; v += (1 - t); }
	// 	else if (&l0.P() == &b1.Q() && &l0.Q() == &b1.R()) { u += (1 - t); }
	// 	else if (&l0.P() == &b1.R() && &l0.Q() == &b1.Q()) { u += t; }
	// 	else if (&l0.P() == &b1.R() && &l0.Q() == &b1.P()) { v += t; }
	// 	else {
	// 		assert(&l0.P() == &b1.P() && &l0.Q() == &b1.R());
	// 		v += (1-t);
	// 	}

	// 	u /= 2;
	// 	v /= 2;

	// 	return new implicitPoint3D_BPT(b1.P(), b1.Q(), b1.R(), v, u);
	// }

	// // BB
	// if (v0->isBPT() && v1->isBPT()) {

	// 	std::cout<<std::endl<<"midpoint of   v0 BPT - v1 BPT"<<std::endl;

	// 	const implicitPoint3D_BPT& b0 = v0->toBPT();
	// 	const implicitPoint3D_BPT& b1 = v1->toBPT();
	// 	const explicitPoint3D* b1p = &b1.P(), *b1q = &b1.Q(), *b1r = &b1.R();
	// 	double b1ou = b1.U(), b1ov = b1.V();
	// 	double b1v, b1u;
	// 	if (&b0.P() == b1p && &b0.Q() == b1q && &b0.R() == b1r) { b1v = b1ov; b1u = b1ou; }
	// 	else if (&b0.P() == b1q && &b0.Q() == b1r && &b0.R() == b1p) { b1v = b1ou; b1u = 1.0 - (b1ov + b1ou); }
	// 	else if (&b0.P() == b1r && &b0.Q() == b1p && &b0.R() == b1q) { b1v = 1.0 - (b1ov + b1ou); b1u = b1ov; }
	// 	else {
	// 		std::swap(b1p, b1q);
	// 		std::swap(b1ou, b1ov);
	// 		if (&b0.P() == b1p && &b0.Q() == b1q && &b0.R() == b1r) { b1v = b1ov; b1u = b1ou; }
	// 		else if (&b0.P() == b1q && &b0.Q() == b1r && &b0.R() == b1p) { b1v = b1ou; b1u = 1.0 - (b1ov + b1ou); }
	// 		else if (&b0.P() == b1r && &b0.Q() == b1p && &b0.R() == b1q) { b1v = 1.0 - (b1ov + b1ou); b1u = b1ov; }
	// 		else{ 
	// 			// ip_error("apparently two BPTs have different vertices...\n");
	// 			std::cout<<"[fp_geometry.h - createMidPoint()] apparently two BPTs have different vertices...\n"; exit(1);
	// 		}

	// 		//implicitPoint3D_BPT pr(b0.P(), b0.Q(), b0.R(), (b0.V() + b1v) / 2, (b0.U() + b1u) / 2);
	// 		//std::cout << *v0 << "\n";
	// 		//std::cout << *v1 << "\n";
	// 		//std::cout << pr << "\n";
	// 		//getchar();
	// 	}
	// 	return new implicitPoint3D_BPT(b0.P(), b0.Q(), b0.R(), (b0.V() + b1v) / 2, (b0.U() + b1u) / 2);
	// }

	// No other possibility is supported
	// ip_error("createMidPoint - Should never reach this point!\n");
	std::cout<<"[fp_geometry.h - createMidPoint()] Should never reach this point!\n"; exit(1);
	return NULL;
}

inline implicitPoint3D_BPT* createCircumCenter(
	const pointType* v0, const pointType* v1, const pointType* v2,
	const explicitPoint* rv0, const explicitPoint* rv1, const explicitPoint* rv2) {
	// Approximated circumcenter
	vector3d a(v0), b(v1), c(v2);
	vector3d ac = c - a;
	vector3d ab = b - a;
	vector3d abXac = ab.cross(ac);
	const double cpl = abXac.sq_length();
	vector3d cc0 = (abXac.cross(ab) * ac.sq_length() + ac.cross(abXac) * ab.sq_length()) * (0.5 / cpl);
	vector3d cc = a + cc0;

	return convertToBPT(cc, rv0, rv1, rv2);
}

inline void computeCircumCenter(
	const pointType* v0, const pointType* v1, const pointType* v2, explicitPoint3D& cc) {
	// Approximated circumcenter
	vector3d a(v0), b(v1), c(v2);
	vector3d ac = c - a;
	vector3d ab = b - a;
	vector3d abXac = ab.cross(ac);
	const double cpl = abXac.sq_length();
	vector3d cc0 = (abXac.cross(ab) * ac.sq_length() + ac.cross(abXac) * ab.sq_length()) * (0.5 / cpl);
	vector3d ccv = a + cc0;
	cc = explicitPoint3D(ccv.c[0], ccv.c[1], ccv.c[2]);
}

inline double createTetCircumcenterPrecise(const pointType* v0, const pointType* v1, const pointType* v2, const pointType* v3, double* cc, double offcenter_tr) {
	const vec3d_bf a(v0), b(v1), c(v2), d(v3); // Exact because we are using LNCs and BPTs only

	// Use coordinates relative to point 'a' of the tetrahedron.
	const vec3d_bf ba = b - a;
	const vec3d_bf ca = c - a;
	const vec3d_bf da = d - a;

	// Squares of lengths of the edges incident to 'a'.
	bigfloat len_ba = ba.sq_length();
	bigfloat len_ca = ca.sq_length();
	bigfloat len_da = da.sq_length();

	// Cross products of these edges.
	vec3d_bf cross_cd = ca & da;
	vec3d_bf cross_db = da & ba;
	vec3d_bf cross_bc = ba & ca;

	// Calculate the denominator of the formula.
	bigfloat denominator = 2.0 * (ba * cross_cd).get_d();

	// Calculate offset (from 'a') of circumcenter. Everything is scaled by 'denominator' to avoid divisions
	vec3d_bf cc0;
	cc0.c[0] = (len_ba * cross_cd.c[0] + len_ca * cross_db.c[0] + len_da * cross_bc.c[0]);
	cc0.c[1] = (len_ba * cross_cd.c[1] + len_ca * cross_db.c[1] + len_da * cross_bc.c[1]);
	cc0.c[2] = (len_ba * cross_cd.c[2] + len_ca * cross_db.c[2] + len_da * cross_bc.c[2]);
	bigfloat radius = cc0.sq_length();
	cc0 = cc0 + (a*denominator);

	// Results must be scaled back by dividing by demoninator
	// Rounding comes into play only here
	double dden = denominator.get_d();
	cc[0] = cc0.c[0].get_d() / dden;
	cc[1] = cc0.c[1].get_d() / dden;
	cc[2] = cc0.c[2].get_d() / dden;
	const double sq_radius = radius.get_d() / (dden * dden);

	// Change to use off-centers
	if (offcenter_tr != 0) {
		const vec3d_bf pp[4] = { a, b, c, d };
		vec3d_bf p, q;
		double al, ml = DBL_MAX;
		for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if ((al = pp[i].dist_sq(pp[j]).get_d()) < ml) {
			ml = al;
			p = pp[i];
			q = pp[j];
		}

		if ((sq_radius / ml) > offcenter_tr) { // Verify quality - might be a good tet if remove_slivers is used
			const double mult = 0.6 * (sqrt(offcenter_tr - 0.25) + sqrt(offcenter_tr)); 

			const double l = sqrt(ml) * mult;

			vec3d_bf cbf(cc[0], cc[1], cc[2]);
			vec3d_bf m = (p + q) * 0.5;
			vec3d_bf dn = (cbf - m) * (1.0 / sqrt((cbf - m).sq_length().get_d()));
			cc0 = m + (dn * l);

			if (((m - cbf).sq_length() - l).sgn() > 0) {
				cc[0] = cc0.c[0].get_d();
				cc[1] = cc0.c[1].get_d();
				cc[2] = cc0.c[2].get_d();
			}
		}

	}
	// End change

	return sq_radius;
}

inline double createTetCircumcenterPrecise(const pointType* v0, const pointType* v1, const pointType* v2, const pointType* v3, vector3d& cc, double offcenter_tr)
{
	return createTetCircumcenterPrecise(v0, v1, v2, v3, cc.c, offcenter_tr);
}

// The tet circumcenter is not exact. It should be guaranteed into the mesh, but rounding might make it slightly out.
// In such a case, it is virtually guaranteed that the circumcenter encroaches upon one of the PLC faces,
// meaning that it would not be used.

inline double createTetCircumcenter(const pointType* v0, const pointType* v1, const pointType* v2, const pointType* v3, vector3d& cc)
{
	const vector3d a(v0), b(v1), c(v2), d(v3);

	// Use coordinates relative to point 'a' of the tetrahedron.

	const vector3d ba = b - a;
	const vector3d ca = c - a;
	const vector3d da = d - a;

	// Squares of lengths of the edges incident to 'a'.
	double len_ba = ba.sq_length();
	double len_ca = ca.sq_length();
	double len_da = da.sq_length();

	// Cross products of these edges.
	vector3d cross_cd = ca & da;
	vector3d cross_db = da & ba;
	vector3d cross_bc = ba & ca;

	// Calculate the denominator of the formula.
	double sden = (ba * cross_cd);
	double denominator = 0.5 / sden;

	// Calculate offset (from 'a') of circumcenter.
	cc.c[0] = (len_ba * cross_cd.c[0] + len_ca * cross_db.c[0] + len_da * cross_bc.c[0]) * denominator;
	cc.c[1] = (len_ba * cross_cd.c[1] + len_ca * cross_db.c[1] + len_da * cross_bc.c[1]) * denominator;
	cc.c[2] = (len_ba * cross_cd.c[2] + len_ca * cross_db.c[2] + len_da * cross_bc.c[2]) * denominator;
	double radius = cc.sq_length();
	cc += a;

	if (pointType::inSphere(explicitPoint(cc.c[0], cc.c[1], cc.c[2]), *v0, *v1, *v2, *v3) < 0) 
		return createTetCircumcenterPrecise(v0, v1, v2, v3, cc, 0);

	return radius;
}


inline double tetrahedronEnergy(const vector3d& w1, const vector3d& w2, const vector3d& w3, const vector3d& w4) {
	vector3d e1, e2, e3;

	// Always use smallest point as origin for the three vectors
	// to guarantee coherent results across different vertex permutations
	if (w1 < w2 && w1 < w3 && w1 < w4) { e1 = w2 - w1; e2 = w3 - w1; e3 = w4 - w1; }
	else if (w2 < w1 && w2 < w3 && w2 < w4) { e2 = w1 - w2; e1 = w3 - w2; e3 = w4 - w2; }
	else if (w3 < w1 && w3 < w2 && w3 < w4) { e1 = w1 - w3; e2 = w2 - w3; e3 = w4 - w3; }
	else { e2 = w1 - w4; e1 = w2 - w4; e3 = w3 - w4; }

	//const vector3d e1 = v2 - v1, e2 = v3 - v1, e3 = v4 - v1;
	const double* t1 = e1.c, * t2 = e2.c, * t3 = e3.c;
	const vector3d j1(-t1[0] + t1[1] + t1[2], t1[0] - t1[1] + t1[2], t1[0] + t1[1] - t1[2]);
	const vector3d j2(-t2[0] + t2[1] + t2[2], t2[0] - t2[1] + t2[2], t2[0] + t2[1] - t2[2]);
	const vector3d j3(-t3[0] + t3[1] + t3[2], t3[0] - t3[1] + t3[2], t3[0] + t3[1] - t3[2]);

	const double num = (j1 * j1) + (j2 * j2) + (j3 * j3);
	const double det = j1.tripleProd(j2, j3);
	if (det <= 0) return DBL_MAX;

	return num / cbrt(det * det); // pow(det, (2.0 / 3.0));
}

inline double distanceLineLine(const vector3d& t, const vector3d& A, const vector3d& A1, const vector3d& B1)
{
	vector3d uu1 = (t - A) & (A1 - B1);
	double nom = (A - A1) * (uu1);
	return (nom * nom) / uu1.sq_length();
}

inline double getTetShortestEdgeSqLength(const pointType* v[4]) {
	vector3d pp[4] = { v[0], v[1], v[2], v[3] };
	double al, ml = DBL_MAX;
	for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if ((al = pp[i].dist_sq(pp[j])) < ml) ml = al;
	return ml;
}

inline double getTetShortestHeightSqLength(const pointType* v[4]) {
	// vector3d pp[4] = { v[0], v[1], v[2], v[3] };
	double al, ml = DBL_MAX;
	if ((al = distanceLineLine(v[0], v[1], v[2], v[3])) < ml) ml = al;
	if ((al = distanceLineLine(v[0], v[2], v[1], v[3])) < ml) ml = al;
	if ((al = distanceLineLine(v[0], v[3], v[2], v[1])) < ml) ml = al;
	return ml;
}

inline double getTetCircumSphere(const pointType* v[4], double ccc[3]) {
	vector3d cc;
	double radius = createTetCircumcenter(v[0], v[1], v[2], v[3], cc);
	ccc[0] = cc.c[0];
	ccc[1] = cc.c[1];
	ccc[2] = cc.c[2];
	return radius;
}

inline void getTetBarycenter(const pointType* v[4], double ccc[3]) {
	const vector3d v0(v[0]), v1(v[1]), v2(v[2]), v3(v[3]);
	const double ew = 4.0, iw = 1.0;
	const double w[4] = { v[0]->isExplicit3D() ? ew : iw, v[1]->isExplicit3D() ? ew : iw, v[2]->isExplicit3D() ? ew : iw, v[3]->isExplicit3D() ? ew : iw };
	const vector3d cc = (v0 * w[0] + v1 * w[1] + v2 * w[2] + v3 * w[3]) * (1.0/(w[0] + w[1] + w[2] + w[3]));
	ccc[0] = cc.c[0];
	ccc[1] = cc.c[1];
	ccc[2] = cc.c[2];
}

inline double getTriangleCircumSphere(const pointType* v[3], double ccc[3]) {
	// Approximated circumcenter
	vector3d a(v[0]), b(v[1]), c(v[2]);
	vector3d ac = c - a;
	vector3d ab = b - a;
	vector3d abXac = ab.cross(ac);
	vector3d cc0 = (abXac.cross(ab) * ac.sq_length() + ac.cross(abXac) * ab.sq_length()) * (0.5 / abXac.sq_length());
	vector3d cc = a + cc0;
	ccc[0] = cc.c[0];
	ccc[1] = cc.c[1];
	ccc[2] = cc.c[2];
	return DBL_MAX;
}

inline double getTetCircumSpherePrecise(const pointType* v[4], double ccc[3], double offcenter_tr) {
	vector3d cc;
	double radius = createTetCircumcenterPrecise(v[0], v[1], v[2], v[3], cc, offcenter_tr);

	ccc[0] = cc.c[0];
	ccc[1] = cc.c[1];
	ccc[2] = cc.c[2];

	return radius;
}

// ANGLES

// dihedral angle cosine
inline double dihedralAngleCos(const pointType* v0, const pointType* v1, const pointType* o0, const pointType* o1) {
	const vector3d v0v(v0), v1v(v1);
	const vector3d vav1 = (vector3d(o0) - v1v) & (v1v - v0v); // Normal at t1
	const vector3d vbv1 = (vector3d(o1) - v1v) & (v1v - v0v); // Normal at t2
	return ((vav1 * vbv1) / (sqrt(vav1.sq_length() * vbv1.sq_length())));
}

inline double minTetDihedralAngleCos(const pointType* v0, const pointType* v1, const pointType* v2, const pointType* v3) {
	double m, ma = dihedralAngleCos(v0, v1, v2, v3);
	if ((m = dihedralAngleCos(v0, v2, v1, v3)) < ma) ma = m;
	if ((m = dihedralAngleCos(v0, v3, v1, v2)) < ma) ma = m;
	if ((m = dihedralAngleCos(v1, v2, v0, v3)) < ma) ma = m;
	if ((m = dihedralAngleCos(v1, v3, v0, v2)) < ma) ma = m;
	if ((m = dihedralAngleCos(v2, v3, v0, v1)) < ma) ma = m;
	return -ma;
}

inline double getAngle(double a, double b, double c) {

	assert(a>0.0 && b>0.0 && c>0.0 && "ERROR: negative side\n");
	// assert((a < b+c) && (b < a+c) && (c < a+b) && "ERROR: triangle inequality does not hold\n"); 
	// if((a > b+c) || (b > a+c) || (c > a+b)){
	// 	std::cout<<"ERROR: triangle inequality does not hold\n"; exit(1);
	// } 

	if (a < b) std::swap(a, b); // must be a > b

	double mu;
	if (c < b) mu = c - (a - b);
	else	  mu = b - (a - c);

	double num = ((a - b) + c) * mu;
	double den = (a + (b + c)) * ((a - c) + b);

	assert(((den != 0.0) || (num != 0.0)) && "ERROR: Ivalid angle (NaN, both num and den are 0)\n" );

	if (den == 0.0) return 180.0;
	if (num == 0.0) return 0.0;

	return 180.0 * ((atan(sqrt((num / den))) * 2.0) / M_PI);
}

inline double getAngle_accurate(bigfloat a, bigfloat b, bigfloat c) {

	assert(a.sgn() > 0.0 && b.sgn() > 0.0 && c.sgn() > 0.0 &&
													"ERROR: negative side\n");
	assert( ((a-(b+c)).sgn() < 0.0) && ((b-(a+c)).sgn() < 0.0) && 
			((c-(a+b)).sgn() < 0.0) && "ERROR: triangle inequality violated\n"); 

	if ((a-b).sgn() < 0.0) std::swap(a, b); // must be a > b

	bigfloat mu;
	if ((c-b).sgn() < 0.0) mu = c - (a - b);
	else	  			   mu = b - (a - c);

	bigfloat num = ((a - b) + c) * mu;
	bigfloat den = (a + (b + c)) * ((a - c) + b);

	assert(((den != 0.0) || (num != 0.0)) && 
				"ERROR: Ivalid angle (NaN, both num and den are 0)\n" );

	if (den == 0.0) return 180.0;
	if (num == 0.0) return 0.0;

	bigrational num_br = num, den_br = den; 
	return 180.0 * ((atan((num_br / den_br).get_bigfloat(256).sqrt(256).get_d()) * 2.0) / M_PI);
}

// Returns the angle at v0 between vectors v1-v0 and v2-v0
inline double getAngle_accurate(const pointType* v0, const pointType* v1, const pointType* v2) {
	bigfloat a = ((vec3d_bf(v1) - vec3d_bf(v0)).sq_length()).sqrt(256);
	bigfloat b = ((vec3d_bf(v2) - vec3d_bf(v0)).sq_length()).sqrt(256);
	bigfloat c = ((vec3d_bf(v2) - vec3d_bf(v1)).sq_length()).sqrt(256);
	return getAngle_accurate(a, b, c);
}

inline double getAngle(const pointType* v0, const pointType* v1, const pointType* v2) {
	double a = sqrt((vector3d(v1) - vector3d(v0)).sq_length());
	double b = sqrt((vector3d(v2) - vector3d(v0)).sq_length());
	double c = sqrt((vector3d(v2) - vector3d(v1)).sq_length());
	return getAngle(a, b, c);
}

// OLD VERSION
// inline double minFaceAngle(const pointType* v0, const pointType* v1, const pointType* v2) {
// 	double minang = 0.0;
// 	const vector3d v[3] = { v0, v1, v2 };
// 	for (int i = 0; i < 3; i++) {  // <- DO NOT CHECK ALL FACE ANGLES
// 		for (int j = i + 1; j < 3; j++) {
// 			for (int k = j + 1; k < 3; k++) {
// 				const vector3d e1 = v[k] - v[i];
// 				const vector3d e2 = v[j] - v[i];
// 				const double ang = ((e1 * e2) / (sqrt(e1.sq_length() * e2.sq_length())));
// 				if (ang > minang) minang = ang;
// 			}
// 		}
// 	}

// 	return 180.0 * (acos(minang) / M_PI);
// }
// inline double minFaceAngle(const pointType* v0, const pointType* v1, const pointType* v2) {
// 	const vector3d v[3] = { v0, v1, v2 };
// 	double l0 = (v[2] - v[1]).sq_length(); 
// 	double l1 = (v[2] - v[0]).sq_length();
// 	double l2 = (v[0] - v[1]).sq_length();
// 	if(l1>l0) std::swap(l0,l1);
// 	if(l2>l1) std::swap(l2,l1); // l2 is the shortest edge: the angle between l0 and l1 is the mimumm over the triangular face. 
// 	return getAngle(sqrt(l0),sqrt(l1),sqrt(l2));
// }

static inline void update_minMax(double val, double& min, double& max) {
	if (val < min) min = val; else if (val > max) max = val;
}

void getMinMaxTetFaceAngles(double& min, double& max, const pointType* v[4]) {
	double val0, val1, val2, l0, l1, l2;
	int i0, i1, i2;
	for (int fi = 0; fi < 4; fi++) {
		i0 = (fi + 1) & 3;  i1 = (fi + 2) & 3;  i2 = (fi + 3) & 3;

		l0 = sqrt((vector3d(v[i1]) - vector3d(v[i2])).sq_length());
		l1 = sqrt((vector3d(v[i2]) - vector3d(v[i0])).sq_length());
		l2 = sqrt((vector3d(v[i0]) - vector3d(v[i1])).sq_length());

		val0 = getAngle(l1, l2, l0); update_minMax(val0, min, max);
		val1 = getAngle(l2, l0, l1); update_minMax(val1, min, max);
		val2 = 180.0 - (val0 + val1); update_minMax(val2, min, max);
		//val2 = getAngle(v[i0], v[i1], v[i2]);
		//if(abs(val0+val1+val2-180.0)>0.000000000001) std::cout<<"angle sum = "<<val0+val1+val3<<std::endl;
	}
}

inline double orient2d_val(double ax, double ay, double bx, double by, double cx, double cy) {
	return (ax - cx) * (by - cy) - (bx - cx) * (ay - cy);
}

inline double orient3d_val(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz) {
	double m11 = ax - dx, m12 = ay - dy, m13 = az - dz;
	double m21 = bx - dx, m22 = by - dy, m23 = bz - dz;
	double m31 = cx - dx, m32 = cy - dy, m33 = cz - dz;
	return (m11 * m22 * m33 + m12 * m23 * m31 + m13 * m21 * m32) - (m13 * m22 * m31 + m12 * m21 * m33 + m11 * m23 * m32);
}

inline double orient3d_val(const vector3d& v0, const vector3d& v1, const vector3d& v2, const vector3d& v3) {
	return orient3d_val(v0.c[0], v0.c[1], v0.c[2],
		v1.c[0], v1.c[1], v1.c[2],
		v2.c[0], v2.c[1], v2.c[2],
		v3.c[0], v3.c[1], v3.c[2]);
}

// Returns the cross product (v1-v0) x (v2-v0)
// From: "Lecture Notes on Geometric Robustnes" by J. R. Shewchuk
inline vector3d accurate_cross_prod(const vector3d& v1, const vector3d& v2, const vector3d& v0) {
	double rx = orient2d_val(v1.c[1], v1.c[2], v2.c[1], v2.c[2], v0.c[1], v0.c[2]);
	double ry = orient2d_val(v1.c[2], v1.c[0], v2.c[2], v2.c[0], v0.c[2], v0.c[0]);
	double rz = orient2d_val(v1.c[0], v1.c[1], v2.c[0], v2.c[1], v0.c[0], v0.c[1]);
	return vector3d(rx, ry, rz);
}

inline double getDihedralAngle_withNormalVectors(const pointType* v0, const pointType* v1, const pointType* o0, const pointType* o1) {
	const vector3d v0v(v0), v1v(v1);
	const vector3d vav1 = (vector3d(o0) - v0v) & (v1v - v0v); // Normal at t1
	const vector3d vbv1 = (vector3d(o1) - v0v) & (v1v - v0v); // Normal at t2
	double l1 = sqrt(vav1.sq_length()), l2 = sqrt(vbv1.sq_length());
	double l3 = sqrt((vav1 - vbv1).sq_length());
	return getAngle(l1, l2, l3);
}

inline double getDihedralAngle_withNormalVectors_bf(const pointType* v0, const pointType* v1, const pointType* o0, const pointType* o1) {
	const vec3d_bf v0v(v0), v1v(v1);
	const vec3d_bf d = (v1v - v0v);
	const vec3d_bf vav1 = (vec3d_bf(o0) - v0v) & d; // Normal at t1
	const vec3d_bf vbv1 = (vec3d_bf(o1) - v0v) & d; // Normal at t2
	bigfloat l1 = (vav1.sq_length()), l2 = (vbv1.sq_length());
	bigfloat l3 = ((vav1 - vbv1).sq_length());
	return getAngle(sqrt(l1.get_d()), sqrt(l2.get_d()), sqrt(l3.get_d()));
}

// Given tetrahedron <v0,v1,v2,v3> (oriented such that orient3d_val<0), 
// computes dihedral angle at tet edge <v0,v1>.
// From: "Lecture Notes on Geometric Robustnes" by J. R. Shewchuk
inline double getDihedralAngle(const pointType* v0, const pointType* v1, const pointType* v2, const pointType* v3) {
	const double toll = 7.7715611723761027e-016;
	bool avoid_atan2_below_toll = true;

	const vector3d Ov0(v0), Ov1(v1), Ov2(v2), Ov3(v3);
	double l01 = sqrt((Ov1 - Ov0).sq_length());

	// if (l01 == 0.0) ip_error("[getDihedralAngle] edge is degenerate\n");
	assert( (l01 != 0.0) && "[getDihedralAngle] edge is degenerate\n" );

	double num = -orient3d_val(Ov0, Ov1, Ov2, Ov3) * l01; // 6 * tet volume
	if (num < 0 && abs(num) > toll) {
		int o3d_tet = pointType::orient3D(*v0, *v1, *v2, *v3);
		// if (o3d_tet < 0) ip_error("[detDihedralAngle] invertedTet\n");
		assert( (o3d_tet >= 0) && "[detDihedralAngle] invertedTet\n" );
		
		std::cout << "[detDihedralAngle] WARNING num ( ori3d = " << num << " * len = " << l01 << ") sign corrected according to exact_orient3D\n";
		num = abs(num);
	}

	const vector3d n2 = accurate_cross_prod(Ov1, Ov0, Ov3); // (Ov1-Ov3) x (Ov0-Ov3)
	const vector3d n3 = accurate_cross_prod(Ov0, Ov1, Ov2); // (Ov0-Ov2) x (Ov1-Ov2)
	double den = (n3 * n2) * (-1.0);

	// if (den == 0.0 && num == 0.0) ip_error("[getDihedralAngle] invalid dihedral angle (NaN)\n");
	assert( (den != 0.0 || num != 0.0) && "[getDihedralAngle] invalid dihedral angle (NaN)\n" );
	
	// Use bigfloats when numerics make things unstable
	if (avoid_atan2_below_toll && (abs(num) < toll || abs(den) < toll))
		return getDihedralAngle_withNormalVectors_bf(v0, v1, v2, v3);

	return atan2(num, den) * 180.0 / M_PI; // atan2 should return angles in [0,pi]
}

void getMinMaxTetDihedralAngles(double& min, double& max, const pointType* v[4], int type = 0) {
	// type = 0 -> uses getDihedralAngle()
	// type = 1 -> uses getDihedralAngle_withNormalVectors()
	if (type == 0) {
		update_minMax(getDihedralAngle(v[0], v[1], v[2], v[3]), min, max);
		update_minMax(getDihedralAngle(v[0], v[2], v[3], v[1]), min, max);
		update_minMax(getDihedralAngle(v[0], v[3], v[1], v[2]), min, max);
		update_minMax(getDihedralAngle(v[1], v[2], v[0], v[3]), min, max);
		update_minMax(getDihedralAngle(v[1], v[3], v[2], v[0]), min, max);
		update_minMax(getDihedralAngle(v[2], v[3], v[0], v[1]), min, max);
	}
	else {
		update_minMax(getDihedralAngle_withNormalVectors(v[0], v[1], v[2], v[3]), min, max);
		update_minMax(getDihedralAngle_withNormalVectors(v[0], v[2], v[3], v[1]), min, max);
		update_minMax(getDihedralAngle_withNormalVectors(v[0], v[3], v[1], v[2]), min, max);
		update_minMax(getDihedralAngle_withNormalVectors(v[1], v[2], v[0], v[3]), min, max);
		update_minMax(getDihedralAngle_withNormalVectors(v[1], v[3], v[2], v[0]), min, max);
		update_minMax(getDihedralAngle_withNormalVectors(v[2], v[3], v[0], v[1]), min, max);
	}

}

// PROJECTIONS

bool pointProjectionInInnerTriangle(const pointType* p, const pointType* t0, const pointType* t1, const pointType* t2) {
	const vector3d pv(p), t0v(t0), t1v(t1), t2v(t2);

	const vector3d e2 = t2v - t1v;
	const vector3d tn = (t1v - t0v).cross(e2);

	const vector3d tn2 = (t1v - pv).cross(e2);
	if (tn2 * tn <= 0) return false;
	const vector3d e0 = t0v - t2v;
	const vector3d tn0 = (t2v - pv).cross(e0);
	if (tn0 * tn <= 0) return false;
	const vector3d e1 = t1v - t0v;
	const vector3d tn1 = (t0v - pv).cross(e1);
	if (tn1 * tn <= 0) return false;
	return true;
}

bool pointProjectionInInnerSegment(const pointType* p, const pointType* s0, const pointType* s1) {
	const vector3d pv(p), s0v(s0), s1v(s1);
	const vector3d s = s1v - s0v;
	const vector3d e0 = pv - s0v;
	const vector3d e1 = pv - s1v;

	return (e0 * s) > 0 && (e1 * s) < 0;
}

// DISTANCES

double pointSqDistanceFromLine(const pointType* p, const pointType* s0, const pointType* s1)
{
	const vector3d pv(p), s0v(s0), s1v(s1);
	const vector3d x21 = s1v - s0v;
	const vector3d x10 = s0v - pv;
	const vector3d x = x21 & x10;

	return (x * x) / (x21 * x21);
}

double pointSqDistanceFromPlane(const pointType* p, const pointType* t0, const pointType* t1, const pointType* t2)
{
	const vector3d pv(p), t0v(t0), t1v(t1), t2v(t2);
	const vector3d dirver = (t2v - t1v) & (t1v - t0v);
	const double CA2 = dirver * dirver;
	const double d = (dirver * pv) - (dirver * t0v);

	return (d * d) / CA2;
}

bool lu_decomp_and_solve(double lu[4][4], double* b)
{
	double s[4], X[4], d, p, l, m, t;
	int pi = 0, i, j, k, ps[4] = { 0, 1, 2, 3 };

	for (i = 0; i < 3; i++) {
		l = 0.0;
		for (j = 0; j < 3; j++) if (l < (t = fabs(lu[i][j]))) l = t;
		if (l == 0) return false;
		s[i] = 1.0 / l;
	}

	for (k = 0; k < 2; k++) {
		l = 0.0;
		for (i = k; i < 3; i++) if (l < (t = fabs(lu[ps[i]][k]) * s[ps[i]])) { l = t; pi = i; }
		if (l == 0.0)  return false;

		if (pi != k) std::swap(ps[k], ps[pi]);

		p = lu[ps[k]][k];
		for (i = k + 1; i < 3; i++) if ((lu[ps[i]][k] = m = lu[ps[i]][k] / p) != 0.0)
			for (j = k + 1; j < 3; j++) lu[ps[i]][j] -= m * lu[ps[k]][j];
	}

	if (lu[ps[2]][2] == 0.0) return false;

	for (i = 0; i < 3; i++) X[i] = 0.0;

	for (i = 0; i < 3; i++) {
		d = 0.0;
		for (j = 0; j < i; j++) d += lu[ps[i]][j] * X[j];
		X[i] = b[ps[i]] - d;
	}

	for (i = 2; i >= 0; i--) {
		d = 0.0;
		for (j = i + 1; j < 3; j++)	d += lu[ps[i]][j] * X[j];
		X[i] = (X[i] - d) / lu[ps[i]][i];
	}

	for (i = 0; i < 3; i++) b[i] = X[i];

	return true;
}

// Compute the circumradius using LU decomposition
// Return the squared radius
double circumsphere_ludecomp(const pointType* pap, const pointType* pbp, const pointType* pcp, const pointType* pdp, double* cent)
{
	const vector3d pa(pap), pb(pbp), pc(pcp), pd(pdp);
	double A[4][4], b[4];

	A[0][0] = pb.c[0] - pa.c[0];
	A[0][1] = pb.c[1] - pa.c[1];
	A[0][2] = pb.c[2] - pa.c[2];
	A[1][0] = pc.c[0] - pa.c[0];
	A[1][1] = pc.c[1] - pa.c[1];
	A[1][2] = pc.c[2] - pa.c[2];
	A[2][0] = pd.c[0] - pa.c[0];
	A[2][1] = pd.c[1] - pa.c[1];
	A[2][2] = pd.c[2] - pa.c[2];

	b[0] = 0.5 * (A[0][0] * A[0][0] + A[0][1] * A[0][1] + A[0][2] * A[0][2]);
	b[1] = 0.5 * (A[1][0] * A[1][0] + A[1][1] * A[1][1] + A[1][2] * A[1][2]);
	b[2] = 0.5 * (A[2][0] * A[2][0] + A[2][1] * A[2][1] + A[2][2] * A[2][2]);

	if (!lu_decomp_and_solve(A, b)) return 0.0;

	cent[0] = pa.c[0] + b[0];
	cent[1] = pa.c[1] + b[1];
	cent[2] = pa.c[2] + b[2];

	return (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
}

void getTriangleAngles(const pointType* vi0, const pointType* vi1, const pointType* vi2, double& a0, double& a1, double& a2) {
	double l0, l1, l2;
	l0 = sqrt((vector3d(vi1) - vector3d(vi2)).sq_length());
	l1 = sqrt((vector3d(vi2) - vector3d(vi0)).sq_length());
	l2 = sqrt((vector3d(vi0) - vector3d(vi1)).sq_length());
	a0 = getAngle(l1, l2, l0);
	a1 = getAngle(l2, l0, l1);
	a2 = getAngle(l0, l1, l2);
}

// Computes the cosine of the angle at A of the triangle <A,B,C>
// using the cosine law.
double cosOfAngle_at(const vector3d& A, const vector3d& B, const vector3d& C){
    double a_sq = B.dist_sq( C );
    double b_sq = C.dist_sq( A );
    double c_sq = A.dist_sq( B );
    return (b_sq + c_sq - a_sq) / (2 * sqrt(b_sq * c_sq) );
}
