/*
g++ -O2 -g -fPIC -DPIC -Wall -W -Wno-unused -I src -I contrib/epee/include/ -I external/easylogging++/ bulletproofs.cc build/debug/src/ringct/libringct.so build/debug/src/crypto/libcncrypto.so build/debug/src/common/libcommon.so build/debug/contrib/epee/src/libepee*.a build/debug/external/easylogging++/libeasylogging.so -lboost_system -lboost_thread -lgmp  -lstdc++
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

//#define DEBUG_BP
#define trace() printf("trace: %u\n", __LINE__)
#define PERF_TIMER_START(name) tools::PerformanceTimer *pt_##name = new tools::PerformanceTimer(#name, 1000000, el::Level::Info)
#define PERF_TIMER_END(name) do { delete pt_##name; pt_##name = NULL; } while(0)

    //package how.monero.hodl.ringSignature;
     
    //import how.monero.hodl.crypto.Curve25519Point;
    //import how.monero.hodl.crypto.Scalar;
    //import how.monero.hodl.crypto.CryptoUtil;
    //import how.monero.hodl.util.ByteUtil;
    //import java.math.BigInteger;
    //import how.monero.hodl.util.VarInt;
    //import java.util.Random;
     
    //import static how.monero.hodl.crypto.Scalar.randomScalar;
    //import static how.monero.hodl.crypto.CryptoUtil.*;
    //import static how.monero.hodl.util.ByteUtil.*;
     
    //public class Bulletproof
    //{
    //    private static int N;
static constexpr size_t N = 64;
    //    private static Curve25519Point G;
    //    private static Curve25519Point H;
    //    private static Curve25519Point[] Gi;
    //    private static Curve25519Point[] Hi;
static rct::key Hi[N], Gi[N];
static ge_dsmp Gprecomp[64], Hprecomp[64];
static rct::keyV p2;
static rct::key two()
{
  static const rct::key TWO = { {0x02, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00 , 0x00, 0x00, 0x00,0x00  } };
  return TWO;
}
     
    //    public static class ProofTuple
    //    {
    //        private Curve25519Point V;
    //        private Curve25519Point A;
    //        private Curve25519Point S;
    //        private Curve25519Point T1;
    //        private Curve25519Point T2;
    //        private Scalar taux;
    //        private Scalar mu;
    //        private Scalar[] l;
    //        private Scalar[] r;
     
    //        public ProofTuple(Curve25519Point V, Curve25519Point A, Curve25519Point S, Curve25519Point T1, Curve25519Point T2, Scalar taux, Scalar mu, Scalar[] l, Scalar[] r)
    //        {
    //            this.V = V;
    //            this.A = A;
    //            this.S = S;
    //            this.T1 = T1;
    //            this.T2 = T2;
    //            this.taux = taux;
    //            this.mu = mu;
    //            this.l = l;
    //            this.r = r;
    //        }
    //    }

struct ProofTuple
{
  rct::key V, A, S, T1, T2;
  rct::key taux, mu;
  rct::keyV l, r;

  ProofTuple() {}
  ProofTuple(const rct::key &V, const rct::key &A, const rct::key &S, const rct::key &T1, const rct::key &T2, const rct::key &taux, const rct::key &mu, const rct::keyV &l, const rct::keyV &r):
    V(V), A(A), S(S), T1(T1), T2(T2), taux(taux), mu(mu), l(l), r(r) {}
};
     
    //    /* Given two scalar arrays, construct a vector commitment */
    //    public static Curve25519Point VectorExponent(Scalar[] a, Scalar[] b)
    //    {
    //        Curve25519Point Result = Curve25519Point.ZERO;
    //        for (int i = 0; i < N; i++)
    //        {
    //            Result = Result.add(Gi[i].scalarMultiply(a[i]));
    //            Result = Result.add(Hi[i].scalarMultiply(b[i]));
    //        }
    //        return Result;
    //    }
     
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

    //    /* Given a scalar, construct a vector of powers */
    //    public static Scalar[] VectorPowers(Scalar x)
    //    {
    //        Scalar[] result = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            result[i] = x.pow(i);
    //        }
    //        return result;
    //    }
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
     
    //    /* Given two scalar arrays, construct the inner product */
    //    public static Scalar InnerProduct(Scalar[] a, Scalar[] b)
    //    {
    //        Scalar result = Scalar.ZERO;
    //        for (int i = 0; i < N; i++)
    //        {
    //            result = result.add(a[i].mul(b[i]));
    //        }
    //        return result;
    //    }
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
     
    //    /* Given two scalar arrays, construct the Hadamard product */
    //    public static Scalar[] Hadamard(Scalar[] a, Scalar[] b)
    //    {
    //        Scalar[] result = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            result[i] = a[i].mul(b[i]);
    //        }
    //        return result;
    //    }
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
     
    //    /* Add two vectors */
    //    public static Scalar[] VectorAdd(Scalar[] a, Scalar[] b)
    //    {
    //        Scalar[] result = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            result[i] = a[i].add(b[i]);
    //        }
    //        return result;
    //    }
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
     
    //    /* Subtract two vectors */
    //    public static Scalar[] VectorSubtract(Scalar[] a, Scalar[] b)
    //    {
    //        Scalar[] result = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            result[i] = a[i].sub(b[i]);
    //        }
    //        return result;
    //    }
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
     
    //    /* Multiply a scalar and a vector */
    //    public static Scalar[] VectorScalar(Scalar[] a, Scalar x)
    //    {
    //        Scalar[] result = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            result[i] = a[i].mul(x);
    //        }
    //        return result;
    //    }
static rct::keyV vector_scalar(const rct::keyV &a, const rct::key &x)
{
  rct::keyV res(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    sc_mul(res[i].bytes, a[i].bytes, x.bytes);
  }
  return res;
}
     
    //    /* Compute the inverse of a scalar, the stupid way */
    //    public static Scalar Invert(Scalar x)
    //    {
    //        Scalar inverse = new Scalar(x.toBigInteger().modInverse(CryptoUtil.l));
    //        assert x.mul(inverse).equals(Scalar.ONE);
     
    //        return inverse;
    //    }
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
     
    //    /* Given a value v (0..2^N-1) and a mask gamma, construct a range proof */
    //    public static ProofTuple PROVE(Scalar v, Scalar gamma)
    //    {
    //        Curve25519Point V = G.scalarMultiply(v).add(H.scalarMultiply(gamma));
    //       
    //        // PAPER LINES 36-37
    //        Scalar[] aL = new Scalar[N];
    //        Scalar[] aR = new Scalar[N];
static ProofTuple PROVE(uint64_t v, const rct::key &gamma)
{
  rct::key V;
  rct::keyV aL(N), aR(N);

  // vG + gammaH
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
     
    //        BigInteger tempV = v.toBigInteger();
    //        for (int i = N-1; i >= 0; i--)
    //        {
    //            BigInteger basePow = BigInteger.valueOf(2).pow(i);
    //            if (tempV.divide(basePow).equals(BigInteger.ZERO))
    //            {
    //                aL[i] = Scalar.ZERO;
    //            }
    //            else
    //            {
    //                aL[i] = Scalar.ONE;
    //                tempV = tempV.subtract(basePow);
    //            }
    //            aR[i] = aL[i].sub(Scalar.ONE);
    //        }
#if 0
  uint64_t tempV = v;
  for (size_t i = N; i-- > 0; )
  {
    //MGINFO("try1("<<v<<"): "<< i << ": " << (tempV / (((uint64_t)1)<<i) ));
    if (tempV / (((uint64_t)1)<<i) == 0) // TODO: replace with &
    {
      aL[i] = rct::zero();
    }
    else
    {
      aL[i] = rct::identity();
      tempV -= ((uint64_t)1) << i;
    }
    sc_sub(aR[i].bytes, aL[i].bytes, rct::identity().bytes);
  }
#endif

  for (size_t i = N; i-- > 0; )
  {
    //MGINFO("try2("<<v<<"): "<< i << ": " << (v&(((uint64_t)1)<<i)));
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

     
    //        // DEBUG: Test to ensure this recovers the value
    //        BigInteger test_aL = BigInteger.ZERO;
    //        BigInteger test_aR = BigInteger.ZERO;
    //        for (int i = 0; i < N; i++)
    //        {
    //            if (aL[i].equals(Scalar.ONE))
    //                test_aL = test_aL.add(BigInteger.valueOf(2).pow(i));
    //            if (aR[i].equals(Scalar.ZERO))
    //                test_aR = test_aR.add(BigInteger.valueOf(2).pow(i));
    //        }
    //        assert test_aL.equals(v.toBigInteger());
    //        assert test_aR.equals(v.toBigInteger());
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
     
    //        // PAPER LINES 38-39
    //        Scalar alpha = randomScalar();
    //        Curve25519Point A = VectorExponent(aL,aR).add(H.scalarMultiply(alpha));
  rct::key alpha = rct::skGen();
  rct::key ve = vector_exponent(aL, aR);
  rct::key A;
  rct::addKeys(A, ve, rct::scalarmultKey(rct::H, alpha));
     
    //        // PAPER LINES 40-42
    //        Scalar[] sL = new Scalar[N];
    //        Scalar[] sR = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            sL[i] = randomScalar();
    //            sR[i] = randomScalar();
    //        }
    //        Scalar rho = randomScalar();
    //        Curve25519Point S = VectorExponent(sL,sR).add(H.scalarMultiply(rho));
  rct::keyV sL = rct::skvGen(N), sR = rct::skvGen(N);
  rct::key rho = rct::skGen();
  ve = vector_exponent(sL, sR);
  rct::key S;
  rct::addKeys(S, ve, rct::scalarmultKey(rct::H, rho));
     
    //        // PAPER LINES 43-45
    //        Scalar y = hashToScalar(concat(A.toBytes(),S.toBytes()));
    //        Scalar z = hashToScalar(y.bytes);
  rct::keyV hashed;
  hashed.push_back(A);
  hashed.push_back(S);
  rct::key y = rct::hash_to_scalar(hashed);
  rct::key z = rct::hash_to_scalar(y);
     
    //        Scalar t0 = Scalar.ZERO;
    //        Scalar t1 = Scalar.ZERO;
    //        Scalar t2 = Scalar.ZERO;
  rct::key t0 = rct::zero();
  rct::key t1 = rct::zero();
  rct::key t2 = rct::zero();
           
  static const auto oneN = vector_powers(rct::identity(), N);
  static const auto twoN = vector_powers(two(), N);
  const auto yN = vector_powers(y, N);

    //        t0 = t0.add(z.mul(InnerProduct(VectorPowers(Scalar.ONE),VectorPowers(y))));
  rct::key ip1y = inner_product(oneN, yN);
  rct::key tmp;
  sc_mul(tmp.bytes, z.bytes, ip1y.bytes);
  sc_add(t0.bytes, t0.bytes, tmp.bytes);
    //        t0 = t0.add(z.sq().mul(v));
  rct::key zsq;
  sc_mul(zsq.bytes, z.bytes, z.bytes);
  sc_mul(tmp.bytes, zsq.bytes, sv.bytes);
  sc_add(t0.bytes, t0.bytes, tmp.bytes);
    //        Scalar k = Scalar.ZERO;
    //        k = k.sub(z.sq().mul(InnerProduct(VectorPowers(Scalar.ONE),VectorPowers(y))));
  rct::key k = rct::zero();
  sc_mul(tmp.bytes, zsq.bytes, ip1y.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
    //        k = k.sub(z.pow(3).mul(InnerProduct(VectorPowers(Scalar.ONE),VectorPowers(Scalar.TWO))));
  rct::key zcu;
  sc_mul(zcu.bytes, zsq.bytes, z.bytes);
  static const rct::key ip12 = inner_product(oneN, twoN);
  sc_mul(tmp.bytes, zcu.bytes, ip12.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
    //        t0 = t0.add(k);
  sc_add(t0.bytes, t0.bytes, k.bytes);

    //        // DEBUG: Test the value of t0 has the correct form
    //        Scalar test_t0 = Scalar.ZERO;
    //        test_t0 = test_t0.add(InnerProduct(aL,Hadamard(aR,VectorPowers(y))));
    //        test_t0 = test_t0.add(z.mul(InnerProduct(VectorSubtract(aL,aR),VectorPowers(y))));
    //        test_t0 = test_t0.add(z.sq().mul(InnerProduct(VectorPowers(Scalar.TWO),aL)));
    //        test_t0 = test_t0.add(k);
    //        assert test_t0.equals(t0);
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
     
  const auto HyNsR = hadamard(yN, sR);
  const auto vpIz = vector_scalar(oneN, z);

    //        t1 = t1.add(InnerProduct(VectorSubtract(aL,VectorScalar(VectorPowers(Scalar.ONE),z)),Hadamard(VectorPowers(y),sR)));
  rct::key ip1 = inner_product(vector_subtract(aL, vpIz), HyNsR);
  sc_add(t1.bytes, t1.bytes, ip1.bytes);

    //        t1 = t1.add(InnerProduct(sL,VectorAdd(Hadamard(VectorPowers(y),VectorAdd(aR,VectorScalar(VectorPowers(Scalar.ONE),z))),VectorScalar(VectorPowers(Scalar.TWO),z.sq()))));
  rct::key ip2 = inner_product(sL, vector_add(hadamard(yN, vector_add(aR, vpIz)), vector_scalar(twoN, zsq)));
  sc_add(t1.bytes, t1.bytes, ip2.bytes);

    //        t2 = t2.add(InnerProduct(sL,Hadamard(VectorPowers(y),sR)));
  rct::key ip3 = inner_product(sL, HyNsR);
  sc_add(t2.bytes, t2.bytes, ip3.bytes);
     
    //        // PAPER LINES 47-48
    //        Scalar tau1 = randomScalar();
    //        Scalar tau2 = randomScalar();
  rct::key tau1 = rct::skGen(), tau2 = rct::skGen();

    //        Curve25519Point T1 = G.scalarMultiply(t1).add(H.scalarMultiply(tau1));
    //        Curve25519Point T2 = G.scalarMultiply(t2).add(H.scalarMultiply(tau2));
  rct::key T1 = rct::addKeys(rct::scalarmultBase(t1), rct::scalarmultKey(rct::H, tau1));
  rct::key T2 = rct::addKeys(rct::scalarmultBase(t2), rct::scalarmultKey(rct::H, tau2));
     
    //        // PAPER LINES 49-51
    //        Scalar x = hashToScalar(concat(z.bytes,T1.toBytes(),T2.toBytes()));
  hashed.clear();
  hashed.push_back(z);
  hashed.push_back(T1);
  hashed.push_back(T2);
  rct::key x = rct::hash_to_scalar(hashed);
           
    //        // PAPER LINES 52-53
    //        Scalar taux = Scalar.ZERO;
    //        taux = tau1.mul(x);
    //        taux = taux.add(tau2.mul(x.sq()));
    //        taux = taux.add(gamma.mul(z.sq()));
    //        Scalar mu = x.mul(rho).add(alpha);
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
     
    //        // PAPER LINES 54-57
    //        Scalar[] l = new Scalar[N];
    //        Scalar[] r = new Scalar[N];
     
    //        l = VectorAdd(VectorSubtract(aL,VectorScalar(VectorPowers(Scalar.ONE),z)),VectorScalar(sL,x));
    //        r = VectorAdd(Hadamard(VectorPowers(y),VectorAdd(aR,VectorAdd(VectorScalar(VectorPowers(Scalar.ONE),z),VectorScalar(sR,x)))),VectorScalar(VectorPowers(Scalar.TWO),z.sq()));
  rct::keyV l = vector_add(vector_subtract(aL, vpIz), vector_scalar(sL, x));
  rct::keyV r = vector_add(hadamard(yN, vector_add(aR, vector_add(vpIz, vector_scalar(sR, x)))), vector_scalar(twoN, zsq));
     
    //        // DEBUG: Test if the l and r vectors match the polynomial forms
    //        Scalar test_t = Scalar.ZERO;
    //        test_t = test_t.add(t0).add(t1.mul(x));
    //        test_t = test_t.add(t2.mul(x.sq()));
    //        assert test_t.equals(InnerProduct(l,r));
#ifdef DEBUG_BP
  rct::key test_t = rct::zero();
  sc_mul(tmp.bytes, t1.bytes, x.bytes);
  sc_add(tmp.bytes, tmp.bytes, t0.bytes);
  sc_add(test_t.bytes, test_t.bytes, tmp.bytes);
  sc_mul(tmp.bytes, t2.bytes, xsq.bytes);
  sc_add(test_t.bytes, test_t.bytes, tmp.bytes);
  tmp = inner_product(l, r);
  sc_sub(tmp.bytes, test_t.bytes, tmp.bytes);
  CHECK_AND_ASSERT_THROW_MES(tmp == rct::zero(), "test_t check failed");
#endif

            // PAPER LINE 58
    //        return new ProofTuple(V,A,S,T1,T2,taux,mu,l,r);
    //    }
  return ProofTuple(V, A, S, T1, T2, taux, mu, l, r);
}
     
    //    /* Given a range proof, determine if it is valid */
    //    public static boolean VERIFY(ProofTuple proof)
    //    {
static bool VERIFY(const ProofTuple &proof)
{
    //        // Reconstruct the challenges
    //        Scalar y = hashToScalar(concat(proof.A.toBytes(),proof.S.toBytes()));
    //        Scalar z = hashToScalar(y.bytes);
    //        Scalar x = hashToScalar(concat(z.bytes,proof.T1.toBytes(),proof.T2.toBytes()));
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
    //        // PAPER LINE 60
    //        Scalar t = InnerProduct(proof.l,proof.r);
  rct::key t = inner_product(proof.l, proof.r);
  PERF_TIMER_END(VERIFY_line_60);
     
  PERF_TIMER_START(VERIFY_line_61);
    //        // PAPER LINE 61
    //        Curve25519Point L61Left = H.scalarMultiply(proof.taux).add(G.scalarMultiply(t));
  rct::key L61Left = rct::addKeys(rct::scalarmultKey(rct::H, proof.taux), rct::scalarmultBase(t));
     
    //        Scalar k = Scalar.ZERO;
    //        k = k.sub(z.sq().mul(InnerProduct(VectorPowers(Scalar.ONE),VectorPowers(y))));
    //        k = k.sub(z.pow(3).mul(InnerProduct(VectorPowers(Scalar.ONE),VectorPowers(Scalar.TWO))));
  rct::key k = rct::zero();
  rct::key ip1y = inner_product(oneN, yN);
  rct::key zsq;
  sc_mul(zsq.bytes, z.bytes, z.bytes);
  rct::key tmp;
  sc_mul(tmp.bytes, zsq.bytes, ip1y.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
  rct::key zcu;
  sc_mul(zcu.bytes, zsq.bytes, z.bytes);
  sc_mul(tmp.bytes, zcu.bytes, ip12.bytes);
  sc_sub(k.bytes, k.bytes, tmp.bytes);
  PERF_TIMER_END(VERIFY_line_61);
     
  PERF_TIMER_START(VERIFY_line_61rl);
    //        Curve25519Point L61Right = G.scalarMultiply(k.add(z.mul(InnerProduct(VectorPowers(Scalar.ONE),VectorPowers(y)))));
  sc_mul(tmp.bytes, z.bytes, ip1y.bytes);
  sc_add(tmp.bytes, k.bytes, tmp.bytes);
  rct::key L61Right = rct::scalarmultBase(tmp);

    //        L61Right = L61Right.add(proof.V.scalarMultiply(z.sq()));
  tmp = rct::scalarmultKey(proof.V, zsq);
  rct::addKeys(L61Right, L61Right, tmp);

    //        L61Right = L61Right.add(proof.T1.scalarMultiply(x));
  tmp = rct::scalarmultKey(proof.T1, x);
  rct::addKeys(L61Right, L61Right, tmp);

    //        L61Right = L61Right.add(proof.T2.scalarMultiply(x.sq()));
  rct::key xsq;
  sc_mul(xsq.bytes, x.bytes, x.bytes);
  tmp = rct::scalarmultKey(proof.T2, xsq);
  rct::addKeys(L61Right, L61Right, tmp);
  PERF_TIMER_END(VERIFY_line_61rl);

    //        if (!L61Right.equals(L61Left))
    //        {
    //            return false;
    //        }
  if (!(L61Right == L61Left))
  {
    MGINFO("VERIFY failure at step 1");
    return false;
  }
     
  PERF_TIMER_START(VERIFY_line_62);
    //        // PAPER LINE 62
    //        Curve25519Point P = Curve25519Point.ZERO;
    //        P = P.add(proof.A);
    //        P = P.add(proof.S.scalarMultiply(x));
  rct::key P = rct::identity();
  rct::addKeys(P, P, proof.A);
  rct::addKeys(P, P, rct::scalarmultKey(proof.S, x));
           
    //        Scalar[] Gexp = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //            Gexp[i] = Scalar.ZERO.sub(z);
  rct::keyV Gexp(N);
  for (size_t i = 0; i < N; ++i)
    sc_sub(Gexp[i].bytes, rct::zero().bytes, z.bytes);
  PERF_TIMER_END(VERIFY_line_62);
     
  PERF_TIMER_START(VERIFY_line_62_2);
    //        Scalar[] Hexp = new Scalar[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            Hexp[i] = Scalar.ZERO;
    //            Hexp[i] = Hexp[i].add(z.mul(y.pow(i)));
    //            Hexp[i] = Hexp[i].add(z.sq().mul(Scalar.TWO.pow(i)));
    //            Hexp[i] = Hexp[i].mul(Invert(y).pow(i));
    //        }
    //        P = P.add(VectorExponent(Gexp,Hexp));
  rct::keyV Hexp(N);
  rct::key pytmp = rct::identity();
  rct::keyV pyinv(N);
  pyinv[0] = rct::identity();
  const rct::key yinv = invert(y);
  PERF_TIMER_END(VERIFY_line_62_2);

  PERF_TIMER_START(VERIFY_line_62_3);
  for (size_t i = 0; i < N; ++i)
  {
    sc_mul(Hexp[i].bytes, z.bytes, pytmp.bytes);
    sc_mul(tmp.bytes, zsq.bytes, p2[i].bytes);
    sc_add(Hexp[i].bytes, Hexp[i].bytes, tmp.bytes);
    sc_mul(Hexp[i].bytes, Hexp[i].bytes, pyinv[i].bytes);

    if (i != N-1)
    {
      // TODO: can save a mul on the first step on those
      sc_mul(pytmp.bytes, pytmp.bytes, y.bytes);
      sc_mul(pyinv[i+1].bytes, pyinv[i].bytes, yinv.bytes);
    }
  }
  PERF_TIMER_END(VERIFY_line_62_3);

  PERF_TIMER_START(VERIFY_line_62_4);
  rct::addKeys(P, P, vector_exponent(Gexp, Hexp));
  PERF_TIMER_END(VERIFY_line_62_4);
     
  PERF_TIMER_START(VERIFY_line_63);
    //        // PAPER LINE 63
    //        for (int i = 0; i < N; i++)
    //        {
    //            Hexp[i] = Scalar.ZERO;
    //            Hexp[i] = Hexp[i].add(proof.r[i]);
    //            Hexp[i] = Hexp[i].mul(Invert(y).pow(i));
    //        }
  for (size_t i = 0; i < N; ++i)
  {
    sc_mul(Hexp[i].bytes, proof.r[i].bytes, pyinv[i].bytes);
  }
  PERF_TIMER_END(VERIFY_line_63);

  PERF_TIMER_START(VERIFY_line_63_2);
    //        Curve25519Point L63Right = VectorExponent(proof.l,Hexp).add(H.scalarMultiply(proof.mu));
  rct::key L63Right = rct::addKeys(vector_exponent(proof.l, Hexp), rct::scalarmultKey(rct::H, proof.mu));
  PERF_TIMER_END(VERIFY_line_63_2);
     
    //        if (!L63Right.equals(P))
    //        {
    //            return false;
    //        }
  if (!(L63Right == P))
  {
    MGINFO("VERIFY failure at step 2");
    return false;
  }
     
    //        return true;
    //    }
  PERF_TIMER_END(VERIFY);
  return true;
}
     
    //    public static void main(String[] args)
    //    {
    //        // Number of bits in the range
    //        N = 64;
int main(int argc, char **argv)
{
  MGINFO("precomp");
     
    //        // Set the curve base points
    //        G = Curve25519Point.G;
    //        H = Curve25519Point.hashToPoint(G);
    //        Gi = new Curve25519Point[N];
    //        Hi = new Curve25519Point[N];
    //        for (int i = 0; i < N; i++)
    //        {
    //            Gi[i] = getHpnGLookup(i);
    //            Hi[i] = getHpnGLookup(N+i);
    //        }
  rct::key generator = rct::H;
  for (size_t i = 0; i < N; ++i)
  {
    generator = rct::addKeys(generator, generator);
    Hi[i] = rct::scalarmultBase(generator);
    rct::precomp(Hprecomp[i], Hi[i]);
    generator = rct::addKeys(generator, generator);
    Gi[i] = rct::scalarmultBase(generator);
    rct::precomp(Gprecomp[i], Gi[i]);
  }

  rct::key p2tmp = rct::identity();
  for (size_t i = 0; i < N; ++i)
  {
    p2.push_back(p2tmp);
    if (i != N-1)
      sc_mul(p2tmp.bytes, p2tmp.bytes, two().bytes);
  }

  //rct::key M, I = rct::identity();
  //sc_mul(M.bytes, I.bytes, I.bytes);
  //MGINFO("I*I: " << M);

  //rct::key hinv = invert(rct::skGen());
  //rct::key hinv = invert(rct::identity());

    //        // Run a bunch of randomized trials
    //        Random rando = new Random();
    //        int TRIALS = 250;
    //        int count = 0;
     
    //        while (count < TRIALS)
    //        {
    //            long amount = rando.nextLong();
    //            if (amount > Math.pow(2,N)-1 || amount < 0)
    //                continue;
     
    //            ProofTuple proof = PROVE(new Scalar(BigInteger.valueOf(amount)),randomScalar());
    //            if (!VERIFY(proof))
    //                System.out.println("Test failed");
     
    //            count += 1;
    //        }
    //    }
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

    //} 
  MGINFO("done");
  return 0;
}
