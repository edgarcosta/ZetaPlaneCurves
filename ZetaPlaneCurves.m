freeze;
/***********************************************************************************************

  ZetaPlaneCurves.

  Computes LPolynomial of plane curves, potentially singular.

  Code is base on Theorem 3.1 of D. Harvey's "Computing zeta functions of arithmetic schemes"

  intrinsics in this file:

  - LPolynomial(f :: RngMPolElt[FldFin]: KnownFactor := false, corrections := false) -> RngUPolElt
  - LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt


  Edgar Costa, 2020
***********************************************************************************************/


// Given f defining a plane curve in P2
// computes the diferences of point counts between the plane model and its normalisation for GF(p^r)
// at the moment we only support nodal singularities
function NormalizationCorrections(f, up_to_r)
  R3 := Parent(f);
  assert Rank(R3) eq 3;
  assert IsPrimeField(BaseRing(R3));
  p := Characteristic(BaseRing(R3));
  //C := Curve(ProjectiveSpace(R3), f);
  //assert IsNodalCurve(C);
  R2<x,y> := PolynomialRing(BaseRing(R3), 2);

  proj2 := hom< R3 -> R2 | [x, y, 1]>;
  F := proj2(f);
  Fx := Derivative(F, x);
  Fy := Derivative(F, y);
  Re1 := Resultant(Fx, Fy, y);
  Re2 := Resultant(Fx, F, y);
  Re3 := Resultant(Fy, F, y);
  R1<x1> := PolynomialRing(BaseRing(R3));
  proj1 := hom< R2 -> R1 | [x1, 1]>;
  H := GCD([proj1(elt) : elt in [Re1, Re2, Re3]]);

  singular_points := [* *];
  if H ne 1 then
    // the curve is not smooth in the affine coordinate z = 1
    for fac in Factorisation(H) do
      pol, mul := Explode(fac);
      // I only wrote the code to handle nodes at the moment
      assert mul eq 1;
      deg := Degree(pol);
      K<a> := GF(p^deg);
      Ku<u> := PolynomialRing(K);
      K2<lx,ly> := PolynomialRing(K, 2);
      for xroot in Roots(Evaluate(pol, u)) do
        Px, Pxmul := Explode(xroot);
        assert Pxmul eq 1;
        // determine the Y coordinate of the singularity
        possiblePy := Roots(Evaluate(F, [Px, u]));
        for yroot in possiblePy do
          Py, Pymul := Explode(yroot);
          if Pymul gt 1 then
            FatP := Evaluate(F, [lx+ Px, ly + Py]);
            if MonomialCoefficient(FatP, lx) eq 0 and MonomialCoefficient(FatP, ly) eq 0 then
              FatPdegree2 := MonomialCoefficient(FatP, lx^2) * lx^2;
              FatPdegree2 +:= MonomialCoefficient(FatP, ly^2) * ly^2;
              FatPdegree2 +:= MonomialCoefficient(FatP, lx * ly) * lx * ly;
              normal_deg := (3 - #Factorisation(FatPdegree2)) * deg;
              Append(~singular_points, <[Px, Py, 1], deg, normal_deg>);
            end if;
          end if;
        end for;
      end for;
    end for;
  end if;
  // now deal with singularities at infinity
  // x := 1
  // z := 0
  deg := 1;
  K<a> := GF(p^deg);
  Ku<u> := PolynomialRing(K);
  K2<lx, ly> := PolynomialRing(K, 2); // lx will play the role of lz
  possiblePy := Roots(Evaluate(f, [1, u, 0]));
  // a bit repetitive code, but too lazy do anything different
  for yroot in possiblePy do
    Py, Pymul := Explode(yroot);
    if Pymul gt 1 then
      FatP := Evaluate(f, [1, ly + Py, lx + 0]);
      if MonomialCoefficient(FatP, lx) eq 0 and MonomialCoefficient(FatP, ly) eq 0 then
        FatPdegree2 := MonomialCoefficient(FatP, lx^2) * lx^2;
        FatPdegree2 +:= MonomialCoefficient(FatP, ly^2) * ly^2;
        FatPdegree2 +:= MonomialCoefficient(FatP, lx * ly) * lx * ly;
        normal_deg := (3 - #Factorisation(FatPdegree2)) * deg;
        Append(~singular_points, <[1, Py, 0], deg, normal_deg>);
      end if;
    end if;
  end for;
  // now compute the corrections
  corrections := [0 : _ in [1..up_to_r]];
  for pt in singular_points do
    _, deg, ndeg := Explode(pt);
    extension := deg;
    while extension le up_to_r do
      corrections[extension] -:= 1;
      extension +:= deg;
    end while;
    extension := ndeg;
    while extension le up_to_r do
      corrections[extension ] +:= 2;
      extension +:= ndeg;
    end while;
  end for;
  return corrections;
end function;

// returns the LPolynomial of f(x,y,z)=0 wth f homogeneous in GF(p)[x,y,z]
function ZetaPlaneCurve(f)
  return LPolynomial(Curve(ProjectiveSpace(Parent(f)), f));
end function;

// return the LPolynomial of f(x,y,z)=0 mod p with f homogeneous in Z[x,y,z]
function ZetaPlaneCurveMod(f, p)
  return ZetaPlaneCurve(ChangeRing(f, GF(p)));
end function;


// Counts points on f(x,y,z)=0 with one or more of the coordinates zero
// by Andrew V. Sutherland
function MissingPoints(f)
  assert IsHomogeneous(f) and Rank(Parent(f)) eq 3;
  R<t> := PolynomialRing(BaseRing(Parent(f)));
  rot := func<v, i | [v[1 + ((i+j) mod 3)] : j in [0..2]]>;
  N := &+[#[a : a in Roots(Evaluate(f,rot([0, t, 1], i))) | a[1] ne 0] : i in [0..2]];
  N +:= &+[Evaluate(f,rot([0,0,1],i)) eq 0 select 1 else 0 :i in [0..2]];
  return N;
end function;

// Return the L-polynomial of the nice curve X/Fp given [p^r+1-#C(Fp^r):r in [1..g]]
function TracesToLPolynomial(ts, p, g: KnownFactor:=1)
  R<T> := PolynomialRing(Integers());
  possible := FrobeniusTracesToWeilPolynomials(ts, p, 1, 2*g :KnownFactor:=ReciprocalPolynomial(R!KnownFactor));
  assert #possible eq 1;
  return ReciprocalPolynomial(R!possible[1]);
end function;



// Compute the matrices A_{F^s} = M_s mod p^e (a=1) for s in [0..e], via Lemma 3.2 using brute force
// The dimension of the returned matrices is Binomial(d*s+2,2), where d=Degree(f)
function mats(f, p, e)
    R1 := PolynomialRing(Integers(p^e));
    R2 := PolynomialRing(R1);
    F := Evaluate(ChangeRing(ChangeRing(f,Integers()),Integers(p^e)), [R1.1, R2.1, 1]);
    // compute all the powers of f necessary
    F := F^(p-1);
    fp1s := [1, F, F^2];
    for s in [3..e] do
        Append(~fp1s, F*fp1s[s]);
    end for;
    assert #fp1s eq e + 1;

    // list of exponent vectors of monomial basis for R[x,y,z]_(d*s)
    R := PolynomialRing(Integers(),3);
    d := Degree(f);
    Bds := [[Vector([Degree(m,R.i):i in [1..3]]):m in MonomialsOfDegree(R,d*s)] : s in [0..e]];


    // for g := f^((p-1)s) return M_s
    function mat(g, s, B, p)
        assert s le e;
        cc := [Coefficients(a) : a in Coefficients(g)];
        M :=  Matrix([[
                    (Min(ee) lt 0) or (ee[1]+1 gt #cc) or (ee[2] + 1 gt #cc[ee[1]+1])
                          select 0 else cc[ee[1]+1][ee[2]+1]
                          where ee is [p * u[k] - v[k] : k in [1..3]]
                      : u in B] : v in B]);
        return M;
    end function;
    return [*mat(fp1s[s+1], s, Bds[s+1], p) : s in [0..e] *];
end function;




// Compute #X_f(F_p^r) mod p^e, where X_f is the plane curve f(x,y,z) = 0 using the trace formula (Theorem 3.1)
// see ApproximateNumberOfPoints for possible improvments
function points_trace_formula(f, p, r, e : Ms:=[])
  assert p gt 1 + e/r;
  if Ms eq [* *] then
      Ms := mats(f, p, e);
  end if;
  R := Integers(p^e);
  N := R!MissingPoints(ChangeRing(f,GF(p^r))) +
         (R!p^r-R!1)^2*&+[(-1)^s*R!Binomial(e,s)*Trace(Ms[s+1]^r):s in [0..e]];
  return N;
end function;

// Computes the zeta function of f(x,y,z) = 0 over Fp
intrinsic LPolynomial(f :: RngMPolElt: KnownFactor := false, corrections := false) -> RngUPolElt
{
  The L-polynomial of the of the projective normalisation of the curve C defined by the zero locus of f in P^3_Q, where C has at most nodal singularities.
  The L-polynomial is obtained by point counts of normalisation, deduced from the point counts on the plane model obtained via Harvey's trace formula, Theorem 3.1 in "Computing zeta functions of arithmetic schemes". where the matrices are computed in a naive fashion.
  }
  R3 := Parent(f);
  require Rank(R3) eq 3 : "f is expected to be a multivariate polynomial in 3 variables";
  require IsPrimeField(BaseRing(R3)) : "curve must be over prime field";
  p := Characteristic(BaseRing(R3));
  require not IsDivisibleBy(Degree(f), p) : "p divides the degree of f"; // possible to work around
  C := Curve(ProjectiveSpace(R3), f);
  g := Genus(Curve(ProjectiveSpace(Parent(f)), f));
 if KnownFactor cmpeq false then KnownFactor := PolynomialRing(Integers())!1; end if;
  up_to_r := Ceiling((2*g - Degree(KnownFactor))/2);
  if corrections cmpeq false then
    require IsNodalCurve(C) : "the plane curve C does not have  only nodes as singularities, please provide point counts corrections for the normalisation";
    corrections := NormalizationCorrections(ChangeRing(f, GF(p)), g);
  else
    require #corrections ge up_to_r : "need corrections up to GF(p^" cat IntegerToString(up_to_r) cat ")";
  end if;
  // one could decrease 2*g to 2*up_to_r by correcting the point counts with known factor and just computing the unkown factor
  e := Ceiling(Log(p,2*g*p^(up_to_r/2)));
  require p gt 1 : "p is too small"; // only implemented the simpler version of the trace formula
  // Compute M_s for s in [0..e]
  time Ms := mats(f, p, e);
  tmodpe := [p^r+1-Integers()!points_trace_formula(f, p, r, e : Ms:=Ms) - corrections[r]
    : r in [1..up_to_r]];
  lift := func<m, n| 2*a gt n select a-n else a where a:= (m mod n)>;
  t := [lift(elt, p^e) : elt in tmodpe];
  return TracesToLPolynomial(t, p, g: KnownFactor:=KnownFactor);
end intrinsic;


intrinsic LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt
{The L-polynomial of the of the projective normalisation of the curve C defined by the zero locus of f in P^3_GF(p), where C has at most nodal singularities and f in Q[x,y,z].}
 return  LPolynomial(PolynomialRing(GF(p),3)!f: KnownFactor := KnownFactor, corrections := false);
end intrinsic;

