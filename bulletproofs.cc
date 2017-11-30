/*
g++ -O2 -g -fPIC -DPIC -Wall -W -Wno-unused -I src -I contrib/epee/include/ -I external/easylogging++/ bulletproofs.cc build/debug/src/ringct/libringct.so build/debug/src/crypto/libcncrypto.so build/debug/src/common/libcommon.so build/debug/contrib/epee/src/libepee*.a build/debug/external/easylogging++/libeasylogging.so -lboost_filesystem -lboost_system -lboost_thread -lgmp  -lstdc++
LD_LIBRARY_PATH=build/debug/src/ringct/:build/debug/src/crypto/:build/debug/src/common:build/debug/external/easylogging++/ ./a.out [<ntrials>]
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "misc_log_ex.h"
#include "common/perf_timer.h"
extern "C" {
#include "crypto/crypto-ops.h"
}
#include "ringct/rctOps.h"
#include "ringct/rctSigs.h" // for borromean test

//#define DEBUG_BP
#define trace() printf("trace: %u\n", __LINE__)
#define PERF_TIMER_START(name) tools::PerformanceTimer *pt_##name = new tools::PerformanceTimer(#name, 1000000, el::Level::Info)
#define PERF_TIMER_END(name) do { delete pt_##name; pt_##name = NULL; } while(0)

static constexpr size_t N = 64;
static constexpr size_t logN = 6;
static rct::key Hi[N], Gi[N];
static ge_dsmp Gprecomp[64], Hprecomp[64];
static rct::keyV p2;

static rct::key two()
{
  static const rct::key TWO = { {0x02, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00  } };
  return TWO;
}
     
struct ProofTuple
{
  rct::key V, A, S, T1, T2;
  rct::key taux, mu;
  rct::keyV L, R;
  rct::key a, b, t;

  ProofTuple() {}
  ProofTuple(const rct::key &V, const rct::key &A, const rct::key &S, const rct::key &T1, const rct::key &T2, const rct::key &taux, const rct::key &mu, const rct::keyV &L, const rct::keyV &R, const rct::key &a, const rct::key &b, const rct::key &t):
    V(V), A(A), S(S), T1(T1), T2(T2), taux(taux), mu(mu), L(L), R(R), a(a), b(b), t(t) {}

  BEGIN_SERIALIZE_OBJECT()
    FIELD(V)
    FIELD(A)
    FIELD(S)
    FIELD(T1)
    FIELD(T2)
    FIELD(taux)
    FIELD(mu)
    FIELD(L)
    FIELD(R)
    FIELD(a)
    FIELD(b)
    FIELD(t)
  END_SERIALIZE()
};
     
/* Given two scalar arrays, construct a vector commitment */
static rct::key vector_exponent(const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  CHECK_AND_ASSERT_THROW_MES(a.size() <= N, "Incompatible sizes of a and N");
  rct::key res = rct::identity();
  for (size_t i = 0; i < a.size(); ++i)
  {
    rct::key term;
    rct::addKeys3(term, a[i], Gprecomp[i], b[i], Hprecomp[i]);
    rct::addKeys(res, res, term);
  }
  return res;
}

/* Compute a custom vector-scalar commitment */
static rct::key vector_exponent_custom(const rct::keyV &A, const rct::keyV &B, const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(A.size() == B.size(), "Incompatible sizes of A and B");
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  CHECK_AND_ASSERT_THROW_MES(a.size() == A.size(), "Incompatible sizes of a and A");
  CHECK_AND_ASSERT_THROW_MES(a.size() <= N, "Incompatible sizes of a and N");
  rct::key res = rct::identity();
  for (size_t i = 0; i < a.size(); ++i)
  {
    rct::key term;
#if 0
    // we happen to know where A and B might fall, so don't bother checking the rest
    ge_dsmp *Acache = NULL, *Bcache = NULL;
    ge_dsmp Acache_custom[1], Bcache_custom[1];
    if (Gi[i] == A[i])
      Acache = Gprecomp + i;
    else if (i<32 && Gi[i+32] == A[i])
      Acache = Gprecomp + i + 32;
    else
    {
      rct::precomp(Acache_custom[0], A[i]);
      Acache = Acache_custom;
    }
    if (i == 0 && B[i] == Hi[0])
      Bcache = Hprecomp;
    else
    {
      rct::precomp(Bcache_custom[0], B[i]);
      Bcache = Bcache_custom;
    }
    rct::addKeys3(term, a[i], *Acache, b[i], *Bcache);
#else
    ge_dsmp Acache, Bcache;
    rct::precomp(Bcache, B[i]);
    rct::addKeys3(term, a[i], A[i], b[i], Bcache);
#endif
    rct::addKeys(res, res, term);
  }
  return res;
}

/* Given a scalar, construct a vector of powers */
static rct::keyV vector_powers(rct::key x, size_t N)
{
  rct::keyV res(N);
  if (N == 0)
    return res;
  res[0] = rct::identity();
  if (N == 1)
    return res;
  res[1] = x;
  for (size_t i = 2; i < N; ++i)
  {
    sc_mul(res[i].bytes, res[i-1].bytes, x.bytes);
  }
  return res;
}
     
/* Given two scalar arrays, construct the inner product */
static rct::key inner_product(const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  rct::key res = rct::zero();
  for (size_t i = 0; i < a.size(); ++i)
  {
    rct::key m;
    sc_mul(m.bytes, a[i].bytes, b[i].bytes);
    sc_add(res.bytes, res.bytes, m.bytes);
  }
  return res;
}
     
/* Given two scalar arrays, construct the Hadamard product */
static rct::keyV hadamard(const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    sc_mul(res[i].bytes, a[i].bytes, b[i].bytes);
  }
  return res;
}

/* Given two curvepoint arrays, construct the Hadamard product */
static rct::keyV hadamard2(const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    rct::addKeys(res[i], a[i], b[i]);
  }
  return res;
}

/* Add two vectors */
static rct::keyV vector_add(const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    sc_add(res[i].bytes, a[i].bytes, b[i].bytes);
  }
  return res;
}
     
/* Subtract two vectors */
static rct::keyV vector_subtract(const rct::keyV &a, const rct::keyV &b)
{
  CHECK_AND_ASSERT_THROW_MES(a.size() == b.size(), "Incompatible sizes of a and b");
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    sc_sub(res[i].bytes, a[i].bytes, b[i].bytes);
  }
  return res;
}
     
/* Multiply a scalar and a vector */
static rct::keyV vector_scalar(const rct::keyV &a, const rct::key &x)
{
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    sc_mul(res[i].bytes, a[i].bytes, x.bytes);
  }
  return res;
}
     
/* Exponentiate a curve vector by a scalar */
static rct::keyV vector_scalar2(const rct::keyV &a, const rct::key &x)
{
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    rct::scalarmultKey(res[i], a[i], x);
  }
  return res;
}

/* Compute the inverse of a scalar, the stupid way */
static rct::key invert(const rct::key &x)
{
  rct::key inv;
#warning TODO: find a better invert func ? Though it seems pretty fast anyway

  mpz_t X;
  mpz_init(X);
  mpz_import(X, sizeof(rct::key)/sizeof(mp_limb_t), -1, sizeof(mp_limb_t), -1, 0, x.bytes);

  mpz_t L;
  mpz_init(L);
  mpz_import(L, sizeof(rct::key)/sizeof(mp_limb_t), -1, sizeof(mp_limb_t), -1, 0, rct::curveOrder().bytes);

  mpz_t minv;
  mpz_init2(minv, 256+sizeof(mp_limb_t)*8);
  CHECK_AND_ASSERT_THROW_MES(mpz_invert(minv, X, L), "Failed to invert");

  CHECK_AND_ASSERT_THROW_MES(mpz_size(minv) * sizeof(mp_limb_t) == sizeof(inv), "unexpected size of inverse");
  mpz_export(inv.bytes, NULL, -1, sizeof(mp_limb_t), -1, 0, minv);

  mpz_clear(X);
  mpz_clear(L);
  mpz_clear(minv);

#ifdef DEBUG_BP
  rct::key tmp;
  sc_mul(tmp.bytes, inv.bytes, x.bytes);
  CHECK_AND_ASSERT_THROW_MES(tmp == rct::identity(), "invert failed");
#endif
  return inv;
}
     
/* Compute the slice of a vector */
static rct::keyV slice(const rct::keyV &a, size_t start, size_t stop)
{
  CHECK_AND_ASSERT_THROW_MES(start < a.size(), "Invalid start index");
  CHECK_AND_ASSERT_THROW_MES(stop <= a.size(), "Invalid stop index");
  CHECK_AND_ASSERT_THROW_MES(start < stop, "Invalid start/stop indices");
  rct::keyV res(stop - start);
  for (size_t i = start; i < stop; ++i)
  {
    res[i - start] = a[i];
  }
  return res;
}

/* Given a value v (0..2^N-1) and a mask gamma, construct a range proof */
static ProofTuple PROVE(uint64_t v, const rct::key &gamma)
{
  PERF_TIMER_UNIT(PROVE, 1000000);

  rct::key V;
  rct::keyV aL(N), aR(N);

  // vG + gammaH
  PERF_TIMER_START(PROVE_v);
  rct::key sv = rct::zero();
  sv.bytes[0] = v & 255;
  sv.bytes[1] = (v >> 8) & 255;
  sv.bytes[2] = (v >> 16) & 255;
  sv.bytes[3] = (v >> 24) & 255;
  sv.bytes[4] = (v >> 32) & 255;
  sv.bytes[5] = (v >> 40) & 255;
  sv.bytes[6] = (v >> 48) & 255;
  sv.bytes[7] = (v >> 56) & 255;
  rct::addKeys2(V, sv, gamma, rct::H);
  PERF_TIMER_END(PROVE_v);
     
  PERF_TIMER_START(PROVE_aLaR);
  for (size_t i = N; i-- > 0; )
  {
    if (v & (((uint64_t)1)<<i))
    {
      aL[i] = rct::identity();
    }
    else
    {
      aL[i] = rct::zero();
    }
    sc_sub(aR[i].bytes, aL[i].bytes, rct::identity().bytes);
  }
  PERF_TIMER_END(PROVE_aLaR);

     
  // DEBUG: Test to ensure this recovers the value
#ifdef DEBUG_BP
  uint64_t test_aL = 0, test_aR = 0;
  for (size_t i = 0; i < N; ++i)
  {
    if (aL[i] == rct::identity())
      test_aL += ((uint64_t)1)<<i;
    if (aR[i] == rct::zero())
      test_aR += ((uint64_t)1)<<i;
  }
  CHECK_AND_ASSERT_THROW_MES(test_aL == v, "test_aL failed");
  CHECK_AND_ASSERT_THROW_MES(test_aR == v, "test_aR failed");
#endif
     
  PERF_TIMER_START(PROVE_step1);
  // PAPER LINES 38-39
  rct::key alpha = rct::skGen();
  rct::key ve = vector_exponent(aL, aR);
  rct::key A;
  rct::addKeys(A, ve, rct::scalarmultKey(rct::H, alpha));
     
  // PAPER LINES 40-42
  rct::keyV sL = rct::skvGen(N), sR = rct::skvGen(N);
  rct::key rho = rct::skGen();
  ve = vector_exponent(sL, sR);
  rct::key S;
  rct::addKeys(S, ve, rct::scalarmultKey(rct::H, rho));
     
  // PAPER LINES 43-45
  rct::keyV hashed;
  hashed.push_back(A);
  hashed.push_back(S);
  rct::key y = rct::hash_to_scalar(hashed);
  rct::key z = rct::hash_to_scalar(y);
     
  // Polynomial construction before PAPER LINE 46
  rct::key t0 = rct::zero();
  rct::key t1 = rct::zero();
  rct::key t2 = rct::zero();
           
  static const auto oneN = vector_powers(rct::identity(), N);
  static const auto twoN = vector_powers(two(), N);
  const auto yN = vector_powers(y, N);

  rct::key ip1y = inner_product(oneN, yN);
  rct::key tmp;
  sc_mul(tmp.bytes, z.bytes, ip1y.bytes);
  sc_add(t0.bytes, t0.bytes, tmp.bytes);

  rct::key zsq;
  sc_mul(zsq.bytes, z.bytes, z.bytes);
  sc_mul(tmp.bytes, zsq.bytes, sv.bytes);
  sc_add(t0.bytes, t0.bytes, tmp.bytes);

  rct::key k = rct::zero();
  sc_mul(tmp.bytes, zsq.bytes, ip1y.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);

  rct::key zcu;
  sc_mul(zcu.bytes, zsq.bytes, z.bytes);
  static const rct::key ip12 = inner_product(oneN, twoN);
  sc_mul(tmp.bytes, zcu.bytes, ip12.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
  sc_add(t0.bytes, t0.bytes, k.bytes);

  // DEBUG: Test the value of t0 has the correct form
#ifdef DEBUG_BP
  rct::key test_t0 = rct::zero();
  rct::key iph = inner_product(aL, hadamard(aR, yN));
  sc_add(test_t0.bytes, test_t0.bytes, iph.bytes);
  rct::key ips = inner_product(vector_subtract(aL, aR), yN);
  sc_mul(tmp.bytes, z.bytes, ips.bytes);
  sc_add(test_t0.bytes, test_t0.bytes, tmp.bytes);
  rct::key ipt = inner_product(twoN, aL);
  sc_mul(tmp.bytes, zsq.bytes, ipt.bytes);
  sc_add(test_t0.bytes, test_t0.bytes, tmp.bytes);
  sc_add(test_t0.bytes, test_t0.bytes, k.bytes);
  CHECK_AND_ASSERT_THROW_MES(t0 == test_t0, "t0 check failed");
#endif
  PERF_TIMER_END(PROVE_step1);
     
  PERF_TIMER_START(PROVE_step2);
  const auto HyNsR = hadamard(yN, sR);
  const auto vpIz = vector_scalar(oneN, z);

  rct::key ip1 = inner_product(vector_subtract(aL, vpIz), HyNsR);
  sc_add(t1.bytes, t1.bytes, ip1.bytes);

  rct::key ip2 = inner_product(sL, vector_add(hadamard(yN, vector_add(aR, vpIz)), vector_scalar(twoN, zsq)));
  sc_add(t1.bytes, t1.bytes, ip2.bytes);

  rct::key ip3 = inner_product(sL, HyNsR);
  sc_add(t2.bytes, t2.bytes, ip3.bytes);
     
  // PAPER LINES 47-48
  rct::key tau1 = rct::skGen(), tau2 = rct::skGen();

  rct::key T1 = rct::addKeys(rct::scalarmultBase(t1), rct::scalarmultKey(rct::H, tau1));
  rct::key T2 = rct::addKeys(rct::scalarmultBase(t2), rct::scalarmultKey(rct::H, tau2));
     
  // PAPER LINES 49-51
  hashed.clear();
  hashed.push_back(z);
  hashed.push_back(T1);
  hashed.push_back(T2);
  rct::key x = rct::hash_to_scalar(hashed);
           
  // PAPER LINES 52-53
  rct::key taux = rct::zero();
  sc_mul(taux.bytes, tau1.bytes, x.bytes);
  rct::key xsq;
  sc_mul(xsq.bytes, x.bytes, x.bytes);
  sc_mul(tmp.bytes, tau2.bytes, xsq.bytes);
  sc_add(taux.bytes, taux.bytes, tmp.bytes);
  sc_mul(tmp.bytes, gamma.bytes, zsq.bytes);
  sc_add(taux.bytes, taux.bytes, tmp.bytes);
  rct::key mu;
  sc_mul(tmp.bytes, x.bytes, rho.bytes);
  sc_add(mu.bytes, tmp.bytes, alpha.bytes);
     
  // PAPER LINES 54-57
  rct::keyV l = vector_add(vector_subtract(aL, vpIz), vector_scalar(sL, x));
  rct::keyV r = vector_add(hadamard(yN, vector_add(aR, vector_add(vpIz, vector_scalar(sR, x)))), vector_scalar(twoN, zsq));
  PERF_TIMER_END(PROVE_step2);
     
  PERF_TIMER_START(PROVE_step3);
  rct::key t = inner_product(l, r);

  // DEBUG: Test if the l and r vectors match the polynomial forms
#ifdef DEBUG_BP
  rct::key test_t = rct::zero();
  sc_mul(tmp.bytes, t1.bytes, x.bytes);
  sc_add(tmp.bytes, tmp.bytes, t0.bytes);
  sc_add(test_t.bytes, test_t.bytes, tmp.bytes);
  sc_mul(tmp.bytes, t2.bytes, xsq.bytes);
  sc_add(test_t.bytes, test_t.bytes, tmp.bytes);
  CHECK_AND_ASSERT_THROW_MES(test_t == t, "test_t check failed");
#endif

  // PAPER LINES 32-33
  hashed.clear();
  hashed.push_back(x);
  hashed.push_back(taux);
  hashed.push_back(mu);
  hashed.push_back(t);
  rct::key x_ip = rct::hash_to_scalar(hashed);

  // These are used in the inner product rounds
  size_t nprime = N;
  rct::keyV Gprime(N);
  rct::keyV Hprime(N);
  rct::keyV aprime(N);
  rct::keyV bprime(N);
  const rct::key yinv = invert(y);
  rct::key yinvpow = rct::identity();
  for (size_t i = 0; i < N; ++i)
  {
    Gprime[i] = Gi[i];
    Hprime[i] = scalarmultKey(Hi[i], yinvpow);
    sc_mul(yinvpow.bytes, yinvpow.bytes, yinv.bytes);
    aprime[i] = l[i];
    bprime[i] = r[i];
  }
  rct::keyV L(logN);
  rct::keyV R(logN);
  int round = 0;
  rct::keyV w(logN); // this is the challenge x in the inner product protocol
  PERF_TIMER_END(PROVE_step3);

  PERF_TIMER_START(PROVE_step4);
  // PAPER LINE 13
  while (nprime > 1)
  {
    // PAPER LINE 15
    nprime /= 2;

    // PAPER LINES 16-17
    rct::key cL = inner_product(slice(aprime, 0, nprime), slice(bprime, nprime, bprime.size()));
    rct::key cR = inner_product(slice(aprime, nprime, aprime.size()), slice(bprime, 0, nprime));

    // PAPER LINES 18-19
    L[round] = vector_exponent_custom(slice(Gprime, nprime, Gprime.size()), slice(Hprime, 0, nprime), slice(aprime, 0, nprime), slice(bprime, nprime, bprime.size()));
    sc_mul(tmp.bytes, cL.bytes, x_ip.bytes);
    rct::addKeys(L[round], L[round], rct::scalarmultBase(tmp));
    R[round] = vector_exponent_custom(slice(Gprime, 0, nprime), slice(Hprime, nprime, Hprime.size()), slice(aprime, nprime, aprime.size()), slice(bprime, 0, nprime));
    sc_mul(tmp.bytes, cR.bytes, x_ip.bytes);
    rct::addKeys(R[round], R[round], rct::scalarmultBase(tmp));

    // PAPER LINES 21-22
    hashed.clear();
    if (round == 0)
    {
      hashed.push_back(L[0]);
      hashed.push_back(R[0]);
      w[0] = rct::hash_to_scalar(hashed);
    }
    else
    {
      hashed.push_back(w[round - 1]);
      hashed.push_back(L[round]);
      hashed.push_back(R[round]);
      w[round] = rct::hash_to_scalar(hashed);
    }

    // PAPER LINES 24-25
    const rct::key winv = invert(w[round]);
    Gprime = hadamard2(vector_scalar2(slice(Gprime, 0, nprime), winv), vector_scalar2(slice(Gprime, nprime, Gprime.size()), w[round]));
    Hprime = hadamard2(vector_scalar2(slice(Hprime, 0, nprime), w[round]), vector_scalar2(slice(Hprime, nprime, Hprime.size()), winv));

    // PAPER LINES 28-29
    aprime = vector_add(vector_scalar(slice(aprime, 0, nprime), w[round]), vector_scalar(slice(aprime, nprime, aprime.size()), winv));
    bprime = vector_add(vector_scalar(slice(bprime, 0, nprime), winv), vector_scalar(slice(bprime, nprime, bprime.size()), w[round]));

    ++round;
  }
  PERF_TIMER_END(PROVE_step4);

  // PAPER LINE 58 (with inclusions from PAPER LINE 8 and PAPER LINE 20)
  return ProofTuple(V, A, S, T1, T2, taux, mu, L, R, aprime[0], bprime[0], t);
}
     
/* Given a range proof, determine if it is valid */
static bool VERIFY(const ProofTuple &proof)
{
  // Reconstruct the challenges
  PERF_TIMER_START(VERIFY);
  PERF_TIMER_START(VERIFY_start);
  rct::keyV hashed;
  hashed.push_back(proof.A);
  hashed.push_back(proof.S);
  rct::key y = rct::hash_to_scalar(hashed);
  rct::key z = rct::hash_to_scalar(y);
  hashed.clear();
  hashed.push_back(z);
  hashed.push_back(proof.T1);
  hashed.push_back(proof.T2);
  rct::key x = rct::hash_to_scalar(hashed);
  PERF_TIMER_END(VERIFY_start);
     
  PERF_TIMER_START(VERIFY_precalc);
  const auto yN = vector_powers(y, 64);
  static const rct::keyV oneN = vector_powers(rct::identity(), 64);
  static const rct::keyV twoN = vector_powers(two(), 64);
  static const rct::key ip12 = inner_product(oneN, twoN);
  PERF_TIMER_END(VERIFY_precalc);

  PERF_TIMER_START(VERIFY_line_60);
  // Reconstruct the challenges
  hashed.clear();
  hashed.push_back(x);
  hashed.push_back(proof.taux);
  hashed.push_back(proof.mu);
  hashed.push_back(proof.t);
  rct::key x_ip = hash_to_scalar(hashed);
  PERF_TIMER_END(VERIFY_line_60);
     
  PERF_TIMER_START(VERIFY_line_61);
  // PAPER LINE 61
  rct::key L61Left = rct::addKeys(rct::scalarmultKey(rct::H, proof.taux), rct::scalarmultBase(proof.t));
     
  rct::key k = rct::zero();
  rct::key ip1y = inner_product(oneN, yN);
  rct::key zsq;
  sc_mul(zsq.bytes, z.bytes, z.bytes);
  rct::key tmp, tmp2;
  sc_mul(tmp.bytes, zsq.bytes, ip1y.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
  rct::key zcu;
  sc_mul(zcu.bytes, zsq.bytes, z.bytes);
  sc_mul(tmp.bytes, zcu.bytes, ip12.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
  PERF_TIMER_END(VERIFY_line_61);
     
  PERF_TIMER_START(VERIFY_line_61rl);
  sc_mul(tmp.bytes, z.bytes, ip1y.bytes);
  sc_add(tmp.bytes, k.bytes, tmp.bytes);
  rct::key L61Right = rct::scalarmultBase(tmp);

  tmp = rct::scalarmultKey(proof.V, zsq);
  rct::addKeys(L61Right, L61Right, tmp);

  tmp = rct::scalarmultKey(proof.T1, x);
  rct::addKeys(L61Right, L61Right, tmp);

  rct::key xsq;
  sc_mul(xsq.bytes, x.bytes, x.bytes);
  tmp = rct::scalarmultKey(proof.T2, xsq);
  rct::addKeys(L61Right, L61Right, tmp);
  PERF_TIMER_END(VERIFY_line_61rl);

  if (!(L61Right == L61Left))
  {
    MGINFO("VERIFY failure at step 1");
    return false;
  }
     
  PERF_TIMER_START(VERIFY_line_62);
  // PAPER LINE 62
  rct::key P = rct::addKeys(proof.A, rct::scalarmultKey(proof.S, x));
  PERF_TIMER_END(VERIFY_line_62);

  // Compute the number of rounds for the inner product
  size_t rounds = proof.L.size();
  CHECK_AND_ASSERT_THROW_MES(rounds > 0, "Zero rounds");

  PERF_TIMER_START(VERIFY_line_21_22);
  // PAPER LINES 21-22
  // The inner product challenges are computed per round
  rct::keyV w(rounds);
  hashed.clear();
  hashed.push_back(proof.L[0]);
  hashed.push_back(proof.R[0]);
  w[0] = rct::hash_to_scalar(hashed);
  for (size_t i = 1; i < rounds; ++i)
  {
    hashed.clear();
    hashed.push_back(w[i-1]);
    hashed.push_back(proof.L[i]);
    hashed.push_back(proof.R[i]);
    w[i] = rct::hash_to_scalar(hashed);
  }
  PERF_TIMER_END(VERIFY_line_21_22);

  PERF_TIMER_START(VERIFY_line_24_25);
  // Basically PAPER LINES 24-25
  // Compute the curvepoints from G[i] and H[i]
  rct::key inner_prod_G = rct::identity();
  rct::key inner_prod_H = rct::identity();
  rct::key yinvpow = rct::identity();
  rct::key ypow = rct::identity();

  PERF_TIMER_START(VERIFY_line_24_25_invert);
  const rct::key yinv = invert(y);
  rct::keyV winv(rounds);
  for (size_t i = 0; i < rounds; ++i)
    winv[i] = invert(w[i]);
  PERF_TIMER_END(VERIFY_line_24_25_invert);

  rct::keyV g_scalar_v(N), h_scalar_v(N);
  for (size_t i = 0; i < N; ++i)
  {
    // Convert the index to binary IN REVERSE and construct the scalar exponent
    rct::key g_scalar = proof.a;
    rct::key h_scalar;
    sc_mul(h_scalar.bytes, proof.b.bytes, yinvpow.bytes);

    for (size_t j = rounds; j-- > 0; )
    {
      size_t J = w.size() - j - 1;

      if ((i & (((size_t)1)<<j)) == 0)
      {
        sc_mul(g_scalar.bytes, g_scalar.bytes, winv[J].bytes);
        sc_mul(h_scalar.bytes, h_scalar.bytes, w[J].bytes);
      }
      else
      {
        sc_mul(g_scalar.bytes, g_scalar.bytes, w[J].bytes);
        sc_mul(h_scalar.bytes, h_scalar.bytes, winv[J].bytes);
      }
    }

    // Adjust the scalars using the exponents from PAPER LINE 62
    sc_add(g_scalar.bytes, g_scalar.bytes, z.bytes);
    sc_mul(tmp2.bytes, zsq.bytes, twoN[i].bytes);
    sc_mul(tmp.bytes, z.bytes, ypow.bytes);
    sc_add(tmp.bytes, tmp.bytes, tmp2.bytes);
    sc_mul(tmp.bytes, tmp.bytes, yinvpow.bytes);
    sc_sub(h_scalar.bytes, h_scalar.bytes, tmp.bytes);

    g_scalar_v[i] = g_scalar;
    h_scalar_v[i] = h_scalar;

    sc_mul(yinvpow.bytes, yinvpow.bytes, yinv.bytes);
    sc_mul(ypow.bytes, ypow.bytes, y.bytes);
  }
  PERF_TIMER_END(VERIFY_line_24_25);

  PERF_TIMER_START(VERIFY_line_24_25_ak3);
  // Now compute the basepoint's scalar multiplication
  // Each of these could be written as a multiexp operation instead
  CHECK_AND_ASSERT_THROW_MES((N&1)==0, "N must be even");
  for (size_t i = 0; i < N; i += 2)
  {
    rct::addKeys3(tmp, g_scalar_v[i], Gprecomp[i], g_scalar_v[i+1], Gprecomp[i+1]);
    rct::addKeys(inner_prod_G, inner_prod_G, tmp);
    rct::addKeys3(tmp, h_scalar_v[i], Hprecomp[i], h_scalar_v[i+1], Hprecomp[i+1]);
    rct::addKeys(inner_prod_H, inner_prod_H, tmp);
  }
  PERF_TIMER_END(VERIFY_line_24_25_ak3);

  PERF_TIMER_START(VERIFY_line_26);
  // PAPER LINE 26
  rct::key pprime;
  sc_sub(tmp.bytes, rct::zero().bytes, proof.mu.bytes);
  rct::addKeys(pprime, P, rct::scalarmultKey(rct::H, tmp));

  for (size_t i = 0; i < rounds; ++i)
  {
    sc_mul(tmp.bytes, w[i].bytes, w[i].bytes);
    sc_mul(tmp2.bytes, winv[i].bytes, winv[i].bytes);
#if 1
    ge_dsmp cacheL, cacheR;
    rct::precomp(cacheL, proof.L[i]);
    rct::precomp(cacheR, proof.R[i]);
    rct::addKeys3(tmp, tmp, cacheL, tmp2, cacheR);
    rct::addKeys(pprime, pprime, tmp);
#else
    rct::addKeys(pprime, pprime, rct::scalarmultKey(proof.L[i], tmp));
    rct::addKeys(pprime, pprime, rct::scalarmultKey(proof.R[i], tmp2));
#endif
  }
  sc_mul(tmp.bytes, proof.t.bytes, x_ip.bytes);
  rct::addKeys(pprime, pprime, rct::scalarmultBase(tmp));
  PERF_TIMER_END(VERIFY_line_26);

  PERF_TIMER_START(VERIFY_step2_check);
  sc_mul(tmp.bytes, proof.a.bytes, proof.b.bytes);
  sc_mul(tmp.bytes, tmp.bytes, x_ip.bytes);
  tmp = rct::scalarmultBase(tmp);
  rct::addKeys(tmp, tmp, inner_prod_G);
  rct::addKeys(tmp, tmp, inner_prod_H);
  PERF_TIMER_END(VERIFY_step2_check);
  if (!(pprime == tmp))
  {
    MGINFO("VERIFY failure at step 2");
    return false;
  }

  PERF_TIMER_END(VERIFY);
  return true;
}

static rct::key get_exponent(const rct::key &base, size_t idx)
{
  std::string hashed = std::string((const char*)base.bytes, sizeof(base)) + tools::get_varint_data(idx);
  return rct::hashToPoint(rct::hash2rct(crypto::cn_fast_hash(hashed.data(), hashed.size())));
}


static void test_borromean()
{
  using namespace rct;

  //Tests for Borromean signatures
  //#boro true one, false one, C != sum Ci, and one out of the range..
  key64 xv;
  key64 P1v;
  key64 P2v;
  bits indi;

  for (size_t j = 0 ; j < N ; j++) {
    indi[j] = (int)randXmrAmount(2);

    xv[j] = skGen();
    if ( (int)indi[j] == 0 ) {
      scalarmultBase(P1v[j], xv[j]);
    } else {
      addKeys1(P1v[j], xv[j], H2[j]);
    }
    subKeys(P2v[j], P1v[j], H2[j]);
  }

  //#true one
  boroSig bb = genBorromean(xv, P1v, P2v, indi);
  PERF_TIMER_START(borromean);
  (verifyBorromean(bb, P1v, P2v));
  PERF_TIMER_END(borromean);
}

     
int main(int argc, char **argv)
{
  mlog_configure("bulletproofs", true);
  MGINFO("precomp");
     
  for (size_t i = 0; i < N; ++i)
  {
    Hi[i] = get_exponent(rct::H, i);
    rct::precomp(Hprecomp[i], Hi[i]);
    Gi[i] = get_exponent(rct::H, i+4096);
    rct::precomp(Gprecomp[i], Gi[i]);
  }

  int TRIALS = 250;
  if (argc > 1)
  {
    unsigned n = atoi(argv[1]);
    if (n > 0) TRIALS = n;
  }

  MGINFO("proving");
  int count = 0;
  ProofTuple proofs[TRIALS];
  while (count < TRIALS)
  {
    uint64_t amount = crypto::rand<uint64_t>();
    proofs[count++] = PROVE(amount, rct::skGen());
  }

  MGINFO("verifying");
  count = 0;
  while (count < TRIALS)
  {
    if (!VERIFY(proofs[count++]))
      printf("Test failed\n");
  }

  test_borromean();

  MGINFO("done");
  return 0;
}
