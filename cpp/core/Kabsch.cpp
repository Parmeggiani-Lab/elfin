// These Kabsch(mostly just SVD) functions taken from the Hybridization
// protocol of the Rosetta software suite, TMalign.cc

#include <cmath>

#include "Kabsch.hpp"
#include "util.h"
#include "MathUtils.hpp"
#include "../input/JSONParser.hpp"

namespace elfin
{

std::vector<double>
point3fToVector(Point3f const & pt)
{
	std::vector<double> v;

	v.resize(3);
	v.at(0) = pt.x;
	v.at(1) = pt.y;
	v.at(2) = pt.z;

	return v;
}

Matrix<double>
points3fToVectors(Points3f const & pts)
{
	Matrix<double> out;

	out.resize(pts.size());
	for (int i = 0; i < out.size(); i++)
		out.at(i) = point3fToVector(pts.at(i));

	return out;
}

void
resample(
    Points3f & ref,
    Points3f & pts)
{
	const uint N = ref.size();

	// Compute  shape total lengths
	float refTotLen = 0.0f;
	for (int i = 1; i < N; i++)
		refTotLen += ref.at(i).distTo(ref.at(i - 1));

	float ptsTotLen = 0.0f;
	for (int i = 1; i < pts.size(); i++)
		ptsTotLen += pts.at(i).distTo(pts.at(i - 1));

	// Upsample pts
	Points3f resampled;

	// First and last points are the same
	resampled.push_back(pts.at(0));

	float refProp = 0.0f, ptsProp = 0.0f;
	int mpi = 1;
	for (int i = 1; i < pts.size(); i++)
	{
		const Point3f & baseFpPoint = pts.at(i - 1);
		const Point3f & nextFpPoint = pts.at(i);
		const float baseFpProportion = ptsProp;
		const float fpSegment = nextFpPoint.distTo(baseFpPoint)
		                        / ptsTotLen;
		const Vector3f vec = nextFpPoint - baseFpPoint;

		ptsProp += fpSegment;
		while (refProp <= ptsProp && mpi < N)
		{
			const float mpSegment =
			    ref.at(mpi).distTo(ref.at(mpi - 1))
			    / refTotLen;

			if (refProp + mpSegment > ptsProp)
				break;
			refProp += mpSegment;

			const float s = (refProp - baseFpProportion)
			                / fpSegment;
			resampled.push_back(baseFpPoint + (vec * s));

			mpi++;
		}
	}

	// Sometimes the last node is automatically added
	if (resampled.size() < N)
		resampled.push_back(pts.back());

	pts = resampled;
}

// void
// upsample(
//     Points3f & ref,
//     Points3f & pts)
// {
// 	// Upsample the shape with fewer points
// 	const bool refIsLonger = ref.size() > pts.size();
// 	const size_t N = refIsLonger ? ref.size() : pts.size();

// 	// Use a proportion based algorithm because
// 	// we want to assume both shapes are roughly
// 	// the same length, but not exactly
// 	const Points3f & morePoints = refIsLonger ? ref : pts;
// 	const Points3f & fewerPoints = refIsLonger ? pts : ref;

// 	// N === ref.size()

// 	// Compute longer shape total length
// 	float mpTotalLength = 0.0f;
// 	for (int i = 1; i < N; i++)
// 		mpTotalLength += morePoints.at(i).distTo(morePoints.at(i - 1));

// 	float fpTotalLength = 0.0f;
// 	for (int i = 1; i < fewerPoints.size(); i++)
// 		fpTotalLength += fewerPoints.at(i).distTo(fewerPoints.at(i - 1));

// 	// Upsample fewerPoints
// 	Points3f upsampled;

// 	// First and last points are the same
// 	upsampled.push_back(fewerPoints.at(0));

// 	float mpProportion = 0.0f, fpProportion = 0.0f;
// 	int mpi = 1;
// 	for (int i = 1; i < fewerPoints.size(); i++)
// 	{
// 		const Point3f & baseFpPoint = fewerPoints.at(i - 1);
// 		const Point3f & nextFpPoint = fewerPoints.at(i);
// 		const float baseFpProportion = fpProportion;
// 		const float fpSegment = nextFpPoint.distTo(baseFpPoint)
// 		                        / fpTotalLength;
// 		const Vector3f vec = nextFpPoint - baseFpPoint;

// 		fpProportion += fpSegment;
// 		while (mpProportion <= fpProportion && mpi < N)
// 		{
// 			const float mpSegment =
// 			    morePoints.at(mpi).distTo(morePoints.at(mpi - 1))
// 			    / mpTotalLength;

// 			if (mpProportion + mpSegment > fpProportion)
// 				break;
// 			mpProportion += mpSegment;

// 			const float s = (mpProportion - baseFpProportion)
// 			                / fpSegment;
// 			upsampled.push_back(baseFpPoint + (vec * s));

// 			mpi++;
// 		}
// 	}

// 	// Sometimes the last node is automatically added
// 	if (upsampled.size() < N)
// 		upsampled.push_back(fewerPoints.back());

// 	refIsLonger ? (pts = upsampled) : (ref = upsampled);
// }

// Implemetation of Kabsch algoritm for finding the best rotation matrix
// ---------------------------------------------------------------------------
// x    - x(i,m) are coordinates of atom m in set x            (input)
// y    - y(i,m) are coordinates of atom m in set y            (input)
// n    - n is number of atom pairs                            (input)
// mode  - 0:calculate rms only                                (input)
// 1:calculate rms,u,t                                 (takes longer)
// rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
// u    - u(i,j) is   rotation  matrix for best superposition  (output)
// t    - t(i)   is translation vector for best superposition  (output)
bool RosettaKabsch(
    std::vector<std::vector<double>> const & x,
    std::vector<std::vector<double>> const & y,
    int const n,
    int const mode,
    double *rms,
    std::vector<double>& t,
    std::vector<std::vector<double>> & u )
{
	int i, j, m, m1, l, k;
	double e0, rms1, d, h, g;
	double cth, sth, sqrth, p, det, sigma;
	double xc[3], yc[3];
	double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	double sqrt3 = 1.73205080756888, tol = 0.01;
	int ip[] = {0, 1, 3, 1, 2, 4, 3, 4, 5};
	int ip2312[] = {1, 2, 0, 1};

	int a_failed = 0, b_failed = 0;
	double epsilon = 0.00000001;

	//initializtation
	*rms = 0;
	rms1 = 0;
	e0 = 0;
	for ( i = 0; i < 3; i++ ) {
		xc[i] = 0.0;
		yc[i] = 0.0;
		t[i] = 0.0;
		for ( j = 0; j < 3; j++ ) {
			u[i][j] = 0.0;
			r[i][j] = 0.0;
			a[i][j] = 0.0;
			if ( i == j ) {
				u[i][j] = 1.0;
				a[i][j] = 1.0;
			}
		}
	}

	if ( n < 1 ) return false;

	//compute centers for vector sets x, y
	for ( i = 0; i < n; i++ ) {
		xc[0] += x[i][0];
		xc[1] += x[i][1];
		xc[2] += x[i][2];

		yc[0] += y[i][0];
		yc[1] += y[i][1];
		yc[2] += y[i][2];
	}
	for ( i = 0; i < 3; i++ ) {
		xc[i] = xc[i] / n;
		yc[i] = yc[i] / n;
	}

	//compute e0 and matrix r
	for ( m = 0; m < n; m++ ) {
		for ( i = 0; i < 3; i++ ) {
			e0 += (x[m][i] - xc[i]) * (x[m][i] - xc[i]) + \
			      (y[m][i] - yc[i]) * (y[m][i] - yc[i]);
			d = y[m][i] - yc[i];
			for ( j = 0; j < 3; j++ ) {
				r[i][j] += d * (x[m][j] - xc[j]);
			}
		}
	}

	//compute determinat of matrix r
	det = r[0][0] * ( r[1][1] * r[2][2] - r[1][2] * r[2][1] )\
	      - r[0][1] * ( r[1][0] * r[2][2] - r[1][2] * r[2][0] )\
	      + r[0][2] * ( r[1][0] * r[2][1] - r[1][1] * r[2][0] );
	sigma = det;

	//compute tras(r)*r
	m = 0;
	for ( j = 0; j < 3; j++ ) {
		for ( i = 0; i <= j; i++ ) {
			rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
			m++;
		}
	}

	double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
	double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
	                - rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
	det = det * det;

	for ( i = 0; i < 3; i++ ) {
		e[i] = spur;
	}

	if ( spur > 0 ) {
		d = spur * spur;
		h = d - cof;
		g = (spur * cof - det) / 2.0 - spur * h;

		if ( h > 0 ) {
			sqrth = sqrt(h);
			d = h * h * h - g * g;
			if ( d < 0.0 ) d = 0.0;
			d = atan2( sqrt(d), -g ) / 3.0;
			cth = sqrth * cos(d);
			sth = sqrth * sqrt3 * sin(d);
			e[0] = (spur + cth) + cth;
			e[1] = (spur - cth) + sth;
			e[2] = (spur - cth) - sth;

			if ( mode != 0 ) {
				for ( l = 0; l < 3; l = l + 2 ) {
					d = e[l];
					ss[0] = (d - rr[2]) * (d - rr[5])  - rr[4] * rr[4];
					ss[1] = (d - rr[5]) * rr[1]      + rr[3] * rr[4];
					ss[2] = (d - rr[0]) * (d - rr[5])  - rr[3] * rr[3];
					ss[3] = (d - rr[2]) * rr[3]      + rr[1] * rr[4];
					ss[4] = (d - rr[0]) * rr[4]      + rr[1] * rr[3];
					ss[5] = (d - rr[0]) * (d - rr[2])  - rr[1] * rr[1];

					if ( fabs(ss[0]) <= epsilon ) ss[0] = 0.0;
					if ( fabs(ss[1]) <= epsilon ) ss[1] = 0.0;
					if ( fabs(ss[2]) <= epsilon ) ss[2] = 0.0;
					if ( fabs(ss[3]) <= epsilon ) ss[3] = 0.0;
					if ( fabs(ss[4]) <= epsilon ) ss[4] = 0.0;
					if ( fabs(ss[5]) <= epsilon ) ss[5] = 0.0;

					if ( fabs(ss[0]) >= fabs(ss[2]) ) {
						j = 0;
						if ( fabs(ss[0]) < fabs(ss[5]) ) {
							j = 2;
						}
					} else if ( fabs(ss[2]) >= fabs(ss[5]) ) {
						j = 1;
					} else {
						j = 2;
					}

					d = 0.0;
					j = 3 * j;
					for ( i = 0; i < 3; i++ ) {
						k = ip[i + j];
						a[i][l] = ss[k];
						d = d + ss[k] * ss[k];
					}

					if ( d > epsilon ) d = 1.0 / sqrt(d);
					else d = 0.0;
					for ( i = 0; i < 3; i++ ) {
						a[i][l] = a[i][l] * d;
					}
				} //for l

				d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
				if ( (e[0] - e[1]) > (e[1] - e[2]) ) {
					m1 = 2; m = 0;
				} else {
					m1 = 0; m = 2;
				}
				p = 0;
				for ( i = 0; i < 3; i++ ) {
					a[i][m1] = a[i][m1] - d * a[i][m];
					p = p + a[i][m1] * a[i][m1];
				}
				if ( p <= tol ) {
					p = 1.0;
					for ( i = 0; i < 3; i++ ) {
						if ( p < fabs(a[i][m]) ) {
							continue;
						}
						p = fabs( a[i][m] );
						j = i;
					}
					k = ip2312[j];
					l = ip2312[j + 1];
					p = sqrt( a[k][m] * a[k][m] + a[l][m] * a[l][m] );
					if ( p > tol ) {
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m] / p;
						a[l][m1] =  a[k][m] / p;
					} else {
						a_failed = 1;
					}
				} else {
					p = 1.0 / sqrt(p);
					for ( i = 0; i < 3; i++ ) {
						a[i][m1] = a[i][m1] * p;
					}
				}
				if ( a_failed != 1 ) {
					a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
					a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
					a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
				}
			}//if(mode!=0)
		}//h>0

		//compute b anyway
		if ( mode != 0 && a_failed != 1 ) { //a is computed correctly
			//compute b
			for ( l = 0; l < 2; l++ ) {
				d = 0.0;
				for ( i = 0; i < 3; i++ ) {
					b[i][l] = r[i][0] * a[0][l] + r[i][1] * a[1][l] + r[i][2] * a[2][l];
					d = d + b[i][l] * b[i][l];
				}
				if ( d > epsilon ) {
					d = 1.0 / sqrt(d);
				} else {
					d = 0.0;
				}
				for ( i = 0; i < 3; i++ ) {
					b[i][l] = b[i][l] * d;
				}
			}
			d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
			p = 0.0;

			for ( i = 0; i < 3; i++ ) {
				b[i][1] = b[i][1] - d * b[i][0];
				p += b[i][1] * b[i][1];
			}

			if ( p <= tol ) {
				p = 1.0;
				for ( i = 0; i < 3; i++ ) {
					if ( p < fabs(b[i][0]) ) {
						continue;
					}
					p = fabs( b[i][0] );
					j = i;
				}
				k = ip2312[j];
				l = ip2312[j + 1];
				p = sqrt( b[k][0] * b[k][0] + b[l][0] * b[l][0] );
				if ( p > tol ) {
					b[j][1] = 0.0;
					b[k][1] = -b[l][0] / p;
					b[l][1] =  b[k][0] / p;
				} else {
					b_failed = 1;
				}
			} else {
				p = 1.0 / sqrt(p);
				for ( i = 0; i < 3; i++ ) {
					b[i][1] = b[i][1] * p;
				}
			}
			if ( b_failed != 1 ) {
				b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
				b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
				b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
				//compute u
				for ( i = 0; i < 3; i++ ) {
					for ( j = 0; j < 3; j++ ) {
						u[i][j] = b[i][0] * a[j][0] + b[i][1] * a[j][1]\
						          + b[i][2] * a[j][2];
					}
				}
			}

			//compute t
			for ( i = 0; i < 3; i++ ) {
				t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1])\
				       - u[i][2] * xc[2];
			}
		}//if(mode!=0 && a_failed!=1)
	} else {
		//compute t
		for ( i = 0; i < 3; i++ ) {
			t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) - u[i][2] * xc[2];
		}
	}

	//compute rms
	for ( i = 0; i < 3; i++ ) {
		if ( e[i] < 0 ) e[i] = 0;
		e[i] = sqrt( e[i] );
	}
	d = e[2];
	if ( sigma < 0.0 ) {
		d = - d;
	}
	d = (d + e[1]) + e[0];
	rms1 = (e0 - d) - d;
	if ( rms1 < 0.0 ) rms1 = 0.0;

	*rms = rms1;
	return true;
}

// A Wrapper to call the a bit more complicated Rosetta version
bool Kabsch(
    const Points3f & mobile,
    const Points3f & ref,
    Matrix<double> & rot,
    Vector3f & tran,
    double & rms,
    int mode = 1)
{
	Matrix<double> xx = points3fToVectors(mobile);
	Matrix<double> yy = points3fToVectors(ref);

	std::vector<double> tt = point3fToVector(tran);

	if (rot.size() < 3)
		rot.resize(3);
	for (auto & row : rot)
		if (row.size() < 3)
			row.resize(3);

	const size_t n = mobile.size();

	const bool retVal = RosettaKabsch(xx, yy, n, mode, &rms, tt, rot);
	tran = Vector3f(std::vector<float>(tt.begin(), tt.end()));

	return retVal;
}

float
kabschScore(
    const Genes & genes,
    Points3f ref)
{
	// First make a copy of genes into points
	Points3f mobile;
	mobile.resize(genes.size());

	for (int i = 0; i < genes.size(); i++)
	{
		const Point3f & pt = genes.at(i).com();
		mobile.at(i) = pt;
	}

	return kabschScore(mobile, ref);
}

float
kabschScore(
    Points3f mobile,
    Points3f ref)
{
	if (ref.size() != mobile.size())
		resample(ref, mobile);

	// Run Kabsch to get RMS
	Matrix<double> rot;
	Vector3f tran;
	double rms;

	const bool retVal = Kabsch(mobile, ref, rot, tran, rms, 0);

	panic_if(!retVal, "Kabsch failed!\n");

	return rms;
}

int _testKabsch()
{
	using namespace elfin;

	msg("Testing Kabsch\n");
	const Point3f arrA[] = {
		Point3f(4.7008892286345, 42.938597096873, 14.4318130193692),
		Point3f(-20.3679194392227, 27.5712678608402, -12.1390617339732),
		Point3f(24.4692807074156, -1.32083675968276, 31.1580458282477),
		Point3f(-31.1044984967455, -6.41414114190809, 3.28255887994549),
		Point3f(18.6775433365315, -5.32162505701938, -14.9272896423117),
		Point3f(-31.648884426273, -19.3650527983443, 43.9001561999887),
		Point3f(-13.1515403509663, 0.850865538112699, 37.5942811492984),
		Point3f(12.561856072969, 1.07715641721097, 5.01563428984222),
		Point3f(28.0227435151377, 31.7627708322262, 12.2475086001227),
		Point3f(-41.8874231134215, 29.4831416883453, 8.70447045314168),
	};

	const Point3f arrB[] =
	{
		Point3f(-29.2257707266972, -18.8897713349587, 9.48960740086143),
		Point3f(-19.8753669720509, 42.3379642103244, -23.7788252219155),
		Point3f(-2.90766514824093, -6.9792608670416, 10.2843089382083),
		Point3f(-26.9511839788441, -31.5183679875864, 21.1215780433683),
		Point3f(34.4308792695389, 40.4880968679893, -27.825326598276),
		Point3f(-30.5235710432951, 47.9748378356085, -38.2582349144194),
		Point3f(-27.4078219027601, -6.11300268738968, -20.3324126781673),
		Point3f(-32.9291952852141, -38.8880776559401, -18.1221698074118),
		Point3f(-27.2335702183446, -24.1935304087933, -7.58332402861928),
		Point3f(-6.43013158961009, -9.12801538874479, 0.785828466111815),
	};

	const Point3f actualR[] =
	{
		Point3f( 0.523673403299203, -0.276948392922051, -0.805646171923458),
		Point3f(-0.793788382691122, -0.501965361762521, -0.343410511043611),
		Point3f(-0.309299482996081, 0.819347522879342, -0.482704326238996),
	};

	const Vector3f actualTran(-1.08234396236629,
	                          5.08395199432057,
	                          -13.0170407784248);

	Points3f A(arrA, arrA + sizeof(arrA) / sizeof(arrA[0]));
	Points3f B(arrB, arrB + sizeof(arrB) / sizeof(arrB[0]));

	Matrix<double> rot;
	Vector3f tran;
	double rms;

	unsigned int failCount = 0;

	// Test Kabsch rotation and translation
	const bool retVal = Kabsch(A, B, rot, tran, rms);
	msg("Kabsch call retVal: %s\n", retVal ? "ok" : "failed");

	msg("Rot:\n");
	for (int i = 0; i < rot.size(); i++)
	{
		auto const & row = rot.at(i);

		raw("%16.6f %16.6f %16.6f\n",
		    row.at(0), row.at(1), row.at(2));
		if (!Vector3f(row.at(0), row.at(1), row.at(2)).approximates(actualR[i]))
		{
			failCount++;
			err("Rotation test failed: row does not approximate actual rotation row\n");
		}
	}

	msg("Tran: %s\n", tran.toString().c_str());
	if (!tran.approximates(actualTran))
	{
		failCount++;
		err("Translation test failed: does not approximate actual translation\n");
	}

	// Test upsampling
	Points3f Afewer = A;
	Afewer.erase(Afewer.begin() + (Afewer.size() / 2),
	             Afewer.begin() + (Afewer.size() / 2) + 1);

	if (Afewer.size() == B.size())
		die("Afewer and B sizes have not been made different!\n");

	resample(Afewer, B);

	if (Afewer.size() != B.size())
	{
		failCount++;
		err("Upsampling failed: Lengths: Afewer=%d B=%d\n",
		    Afewer.size(), B.size());
	}

	// Load necessary data to setup Gene
	RelaMat relaMat;
	NameIdMap nameIdMap;
	IdNameMap idNameMap;
	RadiiList radiiList;
	JSONParser().parseDB("../../resources/xDB.json", nameIdMap, idNameMap, relaMat, radiiList);

	Gene::setup(&idNameMap);

	// Test Kabsch scoring
	Genes G;
	for (int i = 0; i < A.size(); i++)
		G.push_back(Gene(i, A.at(i)));

	Point3f B0Copy = B.at(0);

	float score = kabschScore(G, B);
	msg("A-B Score: %.10f\n", score);
	if (!float_approximates(score, 7796.9331054688))
	{
		failCount++;
		err("A-B Score test failed\n");
	}

	if (!B.at(0).approximates(B0Copy))
	{
		failCount++;
		err("Scoring const-ness failed: shape B was modified during scoring\n");
	}


	// Test scoring identical shapes
	G.clear();
	for (int i = 0; i < B.size(); i++)
		G.push_back(Gene(i, B.at(i)));

	score = kabschScore(G, B);

	msg("B-B Score: %.10f\n", score);
	if (!float_approximates(score, 0.0))
	{
		failCount++;
		err("B-B self score failed\n");
	}

	// Test shifted shapes
	for (Gene & g : G)
		g.com() += Vector3f(-10, 20, 30);

	score = kabschScore(G, B);

	msg("B-B Shifted Score: %.10f\n", score);
	if (!float_approximates(score, 0.0))
	{
		failCount++;
		err("B-B shifted self score failed\n");
	}

	// Test scoring different sized (sub)shapes
	G.clear();
	for (int i = 0; i < B.size(); i++)
	{
		if (i == B.size() / 2)
			continue;
		G.push_back(Gene(i, B.at(i)));
	}

	score = kabschScore(G, B);

	msg("B[1:]-B score: %.10f\n", score);
	
	if (!float_approximates(score, 650.2928466797))
	{
		err("Resampled score differs\n");
		failCount++;
	}

	// Test verdict
	if (failCount == 0)
		msg("Passed!\n");
	else
		err("Failed! failCount=%d\n", failCount);

	return failCount;
}

} // namespace elfin
