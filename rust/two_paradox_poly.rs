// two_paradox_poly.rs -- LABELED tournaments with domination number gamma >= 3,
//                        computed in POLYNOMIAL time (no caterpillar enumeration).
//
// Computes  L(n) = # labeled tournaments on n vertices with domination number
// gamma >= 3  ( = no dominating pair ), together with the domination-number
// distribution
//     gamma = 1 : n * 2^C(n-1,2)                       (a single dominating vertex)
//     gamma = 3 : L(n)          (== gamma >= 3 exactly for n <= 18)
//     gamma = 2 : 2^C(n,2) - (gamma=1) - (gamma=3)      (by complement)
// The smallest tournament with gamma = 4 has 19 vertices, so the "gamma = 3"
// label equals "gamma >= 3" only for n <= 18; the polynomial L(n) computation
// itself is valid for larger n (used below for the n=16/20/24 timing).
//
// METHOD -- the POLYNOMIAL construction (contrast with two_paradox_count.rs,
// which is single-exponential: it enumerates the Theta(2^k) caterpillars/suns).
// This program computes the SAME per-size weights
//     w_k(t) = sum_{caterpillars C,|C|=k} (k!/|Aut C|)(-1)^{k-1} Within(C) I(C)^t
//     oc_v(t)= sum_{odd-cycle suns S,|S|=v}(v!/|Aut S|)(-1)^{e(S)} Within(S) I(S)^t
// with NO enumeration of caterpillars / forests / suns.  It is a faithful port of
// the validated scratchpad/poly_aggregate.py.  Three ideas:
//   1. Symmetric-power I(C)^t via COUNT-COMPRESSION.  I(C)^t = #independent sets
//      of t disjoint copies of C.  Run the spine left->right tracking only the
//      COUNT c in {0..t} of copies whose current spine vertex is "in" (collapsing
//      2^t -> t+1).  Spine-edge count-transition multiplicity C(t-c', c); leaf
//      factor 2^{l(t-c)}.  This (t+1)-state DP runs in LOCKSTEP with the 4-state
//      Within transfer -- no I(C)^t is ever expanded.
//   2. Integer-EGF labelling.  A reduced spine encoding makes k!/|Aut C| fall out
//      with NO automorphism computation: the m>=2 sum is halved (spine reversal),
//      and the multinomial k!/prod l_i! is produced telescopically by the integer
//      per-vertex weight (k-u)*C(k-u-1, l).  K_2 is the one special base case.
//   3. Cyclic-seam DP for odd-cycle suns.  Fix the seed in-set (size c0); track the
//      pair (a,b) = in-copies inside/outside it around the cycle; close the odd
//      cycle by forcing a_{g-1}=0; divide the labelling by 2g (dihedral D_g).
// The outer forest recurrence D(v,j) (O(v^2), with the 2^j inter-component
// coupling) and L(n) = sum_v C(n,v) 2^C(n-v,2) A_v(n-v) are the already-polynomial
// parts.  The Within 4-state transfer (matrices G, adj(G)=144 G^{-1}, the DS/triple
// leaf channels) and the signed bignum `Big` are reused VERBATIM from
// rust/caterpillar_within.rs (via rust/two_paradox_count.rs).
//
// DP STATE DIMENSIONS (every loop bound is polynomial in n, t; NO 2^k factor):
//   * caterpillar S_k(t): state (j in 0..4, c in 0..=t, u in 0..=k)
//       -> 4*(t+1)*(k+1) = O(k t) states; O(k) spine positions;
//          each transition loops l<=k, c<=t, j'<4  =>  O(k^3 t^2) big-int ops.
//   * odd-cycle oc_v(t): outer c0<=t and g odd <=v; state (a in 0..=t, b in 0..=t,
//       u in 0..=v) -> O(v t^2) states; transition loops l<=v, a<=t, b<=t
//          =>  O(v^4 t^5) big-int ops.
//   * outer forest D(v,j): O(v^2) memoised recurrence.
// Whole of L(n) is poly(n) big-int operations on integers of O(n^2) bits (values
// reach 2^C(18,2)=2^153, so i128 is insufficient -- `Big` is used throughout).
//
// Reproduces L(7)=240, L(8)=109440, L(9)=77274880, L(10)=107731095040,
// L(11)=277644146489600, L(12)=1259251974345052160,
// L(13)=10059739307720354796544, L(14)=144025884971245392244400128,
// L(15)=3764299177453034811662226208768,
// L(16)=182456135219244423561589552542810112.
//
// Usage: two_paradox_poly N        (1 <= N <= 30; gamma=3 label exact for N<=18)
// Build:  rustc -O -C target-cpu=native -C overflow-checks=on \
//             rust/two_paradox_poly.rs -o /tmp/tpp

use std::collections::HashMap;

// ======================= minimal signed big integer =======================
// (verbatim from rust/caterpillar_within.rs, via rust/two_paradox_count.rs.)
#[derive(Clone)]
struct Big { neg: bool, mag: Vec<u64> } // little-endian base 2^64, no trailing zeros
impl Big {
    fn zero() -> Big { Big{neg:false, mag:vec![]} }
    fn norm(mut self) -> Big { while let Some(&0)=self.mag.last(){ self.mag.pop(); } if self.mag.is_empty(){self.neg=false;} self }
    fn from_i128(x: i128) -> Big {
        let neg = x<0; let mut u = (if neg { x.unsigned_abs() } else { x as u128 }) as u128;
        let mut m=vec![]; while u>0 { m.push((u & 0xFFFFFFFFFFFFFFFF) as u64); u>>=64; }
        Big{neg, mag:m}.norm()
    }
    fn is_zero(&self)->bool{ self.mag.is_empty() }
    fn cmp_mag(a:&[u64],b:&[u64])->std::cmp::Ordering{
        use std::cmp::Ordering::*;
        if a.len()!=b.len(){ return a.len().cmp(&b.len()); }
        for i in (0..a.len()).rev(){ if a[i]!=b[i]{ return a[i].cmp(&b[i]); } }
        Equal
    }
    fn add_mag(a:&[u64],b:&[u64])->Vec<u64>{
        let mut r=vec![]; let mut carry=0u128;
        for i in 0..a.len().max(b.len()){
            let av=*a.get(i).unwrap_or(&0) as u128; let bv=*b.get(i).unwrap_or(&0) as u128;
            let s=av+bv+carry; r.push((s & 0xFFFFFFFFFFFFFFFF) as u64); carry=s>>64;
        }
        if carry>0 { r.push(carry as u64); }
        r
    }
    fn sub_mag(a:&[u64],b:&[u64])->Vec<u64>{ // a>=b
        let mut r=vec![]; let mut borrow=0i128;
        for i in 0..a.len(){
            let av=a[i] as i128; let bv=*b.get(i).unwrap_or(&0) as i128;
            let mut s=av-bv-borrow;
            if s<0 { s+=1i128<<64; borrow=1; } else { borrow=0; }
            r.push(s as u64);
        }
        r
    }
    fn add(&self,o:&Big)->Big{
        if self.neg==o.neg {
            Big{neg:self.neg, mag:Big::add_mag(&self.mag,&o.mag)}.norm()
        } else {
            match Big::cmp_mag(&self.mag,&o.mag){
                std::cmp::Ordering::Equal => Big::zero(),
                std::cmp::Ordering::Greater => Big{neg:self.neg, mag:Big::sub_mag(&self.mag,&o.mag)}.norm(),
                std::cmp::Ordering::Less => Big{neg:o.neg, mag:Big::sub_mag(&o.mag,&self.mag)}.norm(),
            }
        }
    }
    fn mul_i128(&self,x:i128)->Big{
        if x==0 || self.is_zero() { return Big::zero(); }
        let neg = self.neg ^ (x<0);
        let mut u = x.unsigned_abs(); // u128
        let lo = (u & 0xFFFFFFFFFFFFFFFF) as u64; u>>=64; let hi = u as u64;
        let mul64 = |m:&[u64], f:u64| -> Vec<u64> {
            let mut r=vec![]; let mut carry=0u128;
            for &mi in m { let p=(mi as u128)*(f as u128)+carry; r.push((p & 0xFFFFFFFFFFFFFFFF) as u64); carry=p>>64; }
            if carry>0 { r.push(carry as u64); }
            r
        };
        let p_lo = mul64(&self.mag, lo);
        let mut res = Big{neg:false, mag:p_lo};
        if hi!=0 {
            let mut p_hi = mul64(&self.mag, hi);
            p_hi.insert(0,0);
            res = res.add(&Big{neg:false, mag:p_hi});
        }
        res.neg=neg; res.norm()
    }
    fn div_small_exact(&self, d:u64)->Big{ // exact division, d>0
        let mut rem=0u128; let mut out=vec![0u64;self.mag.len()];
        for i in (0..self.mag.len()).rev(){
            let cur=(rem<<64) | (self.mag[i] as u128);
            out[i]=(cur/(d as u128)) as u64; rem=cur%(d as u128);
        }
        Big{neg:self.neg, mag:out}.norm()
    }
    fn mul_u64(&self,f:u64)->Big{
        if f==0||self.is_zero(){return Big::zero();}
        let mut r=vec![]; let mut carry=0u128;
        for &mi in &self.mag { let p=(mi as u128)*(f as u128)+carry; r.push((p&0xFFFFFFFFFFFFFFFF) as u64); carry=p>>64; }
        if carry>0 { r.push(carry as u64); }
        Big{neg:self.neg,mag:r}.norm()
    }
    fn to_string(&self)->String{
        if self.is_zero(){ return "0".into(); }
        let mut m=self.mag.clone(); let mut digits=vec![];
        while !m.is_empty(){
            let mut rem=0u128;
            for i in (0..m.len()).rev(){ let cur=(rem<<64)|(m[i] as u128); m[i]=(cur/1000000000000000000u128) as u64; rem=cur%1000000000000000000u128; }
            while let Some(&0)=m.last(){ m.pop(); }
            digits.push(rem as u64);
        }
        let mut s=String::new(); if self.neg { s.push('-'); }
        s.push_str(&digits.last().unwrap().to_string());
        for i in (0..digits.len()-1).rev(){ s.push_str(&format!("{:018}", digits[i])); }
        s
    }
}
// big*big multiply
fn mul_big(a:&Big,b:&Big)->Big{
    if a.is_zero()||b.is_zero(){ return Big::zero(); }
    let mut res=Big::zero();
    for (i,&limb) in b.mag.iter().enumerate(){
        if limb==0 { continue; }
        let mut part=a.mul_u64(limb);
        for _ in 0..i { part.mag.insert(0,0); }
        res=res.add(&part);
    }
    res.neg = a.neg ^ b.neg; res.norm()
}

fn c2(n:i128)->u32{ (n*(n-1)/2) as u32 } // C(n,2) for n>=0
fn pow2(e:u32)->Big{ // 2^e as Big
    let mut m=vec![0u64; (e/64+1) as usize]; m[(e/64) as usize]=1u64<<(e%64); Big{neg:false,mag:m}.norm()
}
fn f_p(l:i128)->Big{ pow2(c2(l)) }               // 2^C(l,2)
fn f_q(l:i128)->Big{ if l<=0 { Big::zero() } else { pow2(c2(l-1)).mul_i128(l) } } // l*2^C(l-1,2)

// ================= Within(caterpillar) transfer matrix =================
// (verbatim from rust/caterpillar_within.rs; leaf-count-word spec.  The DS/triple
//  leaf channels and adj(G) are the black-box Within primitive.)
fn chan(l:i128, c:[i128;8])->Big{
    let a_s=[1i128,2,4,8];
    let mut acc=Big::zero();
    for k in 0..4 {
        if c[k]!=0 {
            let mut t=f_p(l);
            let mut al=1i128; for _ in 0..l { al*=a_s[k]; }
            t=t.mul_i128(al).mul_i128(c[k]);
            acc=acc.add(&t);
        }
        if c[4+k]!=0 {
            let mut t=f_q(l);
            let mut al=1i128; for _ in 0..(l-1).max(0) { al*=a_s[k]; }
            if l>=1 { t=t.mul_i128(al).mul_i128(c[4+k]); acc=acc.add(&t); }
        }
    }
    acc
}
fn star(m:i128)->Big{ f_p(m).add(&f_q(m)) } // 2^C(m,2)+m 2^C(m-1,2)  = Within(K_{1,m})
const D1:[i128;8]=[3,0,0,0, 1,0,0,0];
const D2:[i128;8]=[6,0,0,0, 2,0,0,0];
const D3:[i128;8]=[22,0,0,0, 8,0,0,0];
fn ds(p:i128,q:i128)->Big{   // DS(p,q) = Within([p,q])
    if p==0 { return star(q+1); }
    if q==0 { return star(p+1); }
    let (a,b)= if q<4 {(p,q)} else {(q,p)};
    match b { 1=>chan(a,D1),2=>chan(a,D2),3=>chan(a,D3), _=>unreachable!() }
}
const G00:[i128;8]=[0,2,2,0, 0,0,2,0];
const G01:[i128;8]=[1,3,0,0, 0,1,0,0];
const G02:[i128;8]=[2,6,0,0, 0,2,0,0];
const G03:[i128;8]=[8,22,0,0, 0,8,0,0];
const G11:[i128;8]=[4,0,0,0, 0,0,0,0];
const G12:[i128;8]=[8,0,0,0, 0,0,0,0];
const G13:[i128;8]=[30,0,0,0, 0,0,0,0];
const G22:[i128;8]=[16,0,0,0, 0,0,0,0];
const G23:[i128;8]=[60,0,0,0, 0,0,0,0];
const G33:[i128;8]=[224,0,0,0, 0,0,0,0];
fn triple(a:i128,l:i128,c:i128)->Big{  // Within([a,l,c]), a,c in 0..3
    let (x,y)=(a.min(c),a.max(c));
    let co = match (x,y) {
        (0,0)=>G00,(0,1)=>G01,(0,2)=>G02,(0,3)=>G03,(1,1)=>G11,(1,2)=>G12,(1,3)=>G13,(2,2)=>G22,(2,3)=>G23,(3,3)=>G33,_=>unreachable!()
    };
    chan(l,co)
}
const ADJ:[[i128;4];4]=[ // adj(G) = 144 * G^{-1}
    [0,-48,24,0],
    [-48,-30112,22928,-2088],
    [24,22928,-17416,1584],
    [0,-2088,1584,-144],
];
const DETG:u64=144;

// --- per-spine-vertex operators as functions of the leaf count l ---
fn r_init(l0:i128)->[Big;4]{                 // r(l0)[j] = DS(l0,j)
    [ds(l0,0),ds(l0,1),ds(l0,2),ds(l0,3)]
}
fn w_til(l:i128)->[[Big;4];4]{               // Wtil(l)[j][d] = sum_a triple(a,l,j)*adj(G)[d][a]
    std::array::from_fn(|j| std::array::from_fn(|d| {
        let mut acc=Big::zero();
        for a in 0..4 { acc=acc.add(&triple(a as i128,l,j as i128).mul_i128(ADJ[d][a])); }
        acc
    }))
}
fn rho_til(x:i128)->[Big;4]{                  // rho(x)[d] = sum_a DS(a,x)*adj(G)[d][a]
    std::array::from_fn(|d| {
        let mut acc=Big::zero();
        for a in 0..4 { acc=acc.add(&ds(a as i128,x).mul_i128(ADJ[d][a])); }
        acc
    })
}

// ============================ small integer helpers ============================
fn comb(n:i128,k:i128)->i128{
    if k<0 || k>n { return 0; }
    let k = k.min(n-k);
    let mut r:i128=1;
    for i in 0..k { r = r*(n-i)/(i+1); }
    r
}
fn ipow(base:i128, t:i128)->Big{ // base^t as Big (t>=0)
    let mut r=Big::from_i128(1);
    for _ in 0..t { r=r.mul_i128(base); }
    r
}
// V[key] += val   (sparse accumulate)
fn dp_add(map:&mut HashMap<(i128,i128,i128),Big>, key:(i128,i128,i128), val:Big){
    if val.is_zero(){ return; }
    let e = map.entry(key).or_insert_with(Big::zero);
    let nv = e.add(&val);
    *e = nv;
}

// =====================================================================
//  Part B.  Fused caterpillar DP:  S_k(t) = sum_{labeled caterpillars,|C|=k}
//                                            Within(C) I(C)^t ;  w_k=(-1)^{k-1}S_k.
//  State (j, c, u):  j = Within R-index 0..3, c = count 0..t, u = labels used 0..k.
// =====================================================================
fn sk(k:i128, t:i128)->Big{
    if k<2 { return Big::zero(); }
    if k==2 {                                   // single edge K_2: Within=2, I=3
        return ipow(3, t).mul_i128(2);
    }
    let mut total=Big::zero();
    // ---- m=1 stars K_{1,k-1}: k labelings, I(K_{1,l})=2^l+1 ----
    let l0=k-1;                                  // >= 2 since k>=3
    let base=(1i128<<l0)+1;
    total=total.add(&mul_big(&star(l0).mul_i128(k), &ipow(base, t)));
    // ---- m>=2 caterpillars (needs k>=4) ----
    if k>=4 { total=total.add(&cat_mge2(k, t)); }
    total
}

fn wk(k:i128, t:i128)->Big{
    let s=sk(k,t);
    if (k-1)%2==0 { s } else { s.mul_i128(-1) }
}

fn cat_mge2(k:i128, t:i128)->Big{
    // precompute leaf-dependent objects for l = 0..k-1
    let ku=k as usize;
    let wt:Vec<[[Big;4];4]>=(0..ku).map(|l|w_til(l as i128)).collect();
    let rin:Vec<[Big;4]>=(0..ku).map(|l|r_init(l as i128)).collect();
    let rho:Vec<[Big;4]>=(0..ku).map(|l|rho_til(l as i128)).collect();
    let binom_t:Vec<i128>=(0..=t).map(|c|comb(t,c)).collect();

    // vertex 0 : l0 >= 1, leave room for >=1 more vertex => u0 <= k-1
    let mut v:HashMap<(i128,i128,i128),Big>=HashMap::new();
    let mut l0=1i128;
    while l0<k {
        let u0=1+l0;
        if u0>k-1 { break; }
        let lab0=k*comb(k-1,l0);
        for c0 in 0..=t {
            let sc=pow2((l0*(t-c0)) as u32).mul_i128(binom_t[c0 as usize]).mul_i128(lab0);
            if sc.is_zero() { continue; }
            let ri=&rin[l0 as usize];
            for j in 0..4 {
                if !ri[j].is_zero() { dp_add(&mut v, (j as i128, c0, u0), mul_big(&ri[j], &sc)); }
            }
        }
        l0+=1;
    }

    let mut total=Big::zero();
    total=total.add(&cat_readout(&v, k, t, &rho, 1));   // m=2
    let mut steps=0i128;
    while !v.is_empty() {
        v=cat_interior(&v, k, t, &wt);
        if v.is_empty() { break; }
        steps+=1;
        total=total.add(&cat_readout(&v, k, t, &rho, steps+1)); // m = steps+2
        if v.keys().all(|&(_,_,u)| u>k-2) { break; }
    }
    // reversal symmetry factor 1/2
    total.div_small_exact(2)
}

fn cat_interior(v:&HashMap<(i128,i128,i128),Big>, k:i128, t:i128, wt:&[[[Big;4];4]])->HashMap<(i128,i128,i128),Big>{
    let mut nv:HashMap<(i128,i128,i128),Big>=HashMap::new();
    for (&(jp,cp,u), val) in v.iter() {
        let lmax=k-2-u;                 // leave room for the read-out vertex
        if lmax<0 { continue; }
        for l in 0..=lmax {
            let up=u+1+l;
            let lab=(k-u)*comb(k-u-1,l);
            let wl=&wt[l as usize];
            let base_l=val.mul_i128(lab);
            for c in 0..=t {
                let cmb=comb(t-cp,c);
                if cmb==0 { continue; }
                let sc=pow2((l*(t-c)) as u32).mul_i128(cmb);
                let b=mul_big(&base_l,&sc);
                for j in 0..4 {
                    let w=&wl[j][jp as usize];
                    if !w.is_zero() { dp_add(&mut nv, (j as i128, c, up), mul_big(w,&b)); }
                }
            }
        }
    }
    nv
}

fn cat_readout(v:&HashMap<(i128,i128,i128),Big>, k:i128, t:i128, rho:&[[Big;4]], power:i128)->Big{
    let mut s=Big::zero();
    for (&(jp,cp,u), val) in v.iter() {
        let l=k-1-u;
        if l<1 { continue; }
        let lab=(k-u)*comb(k-u-1,l);      // C(k-u-1,l)=1, lab=k-u
        let rl=&rho[l as usize][jp as usize];
        if rl.is_zero() { continue; }
        let base=mul_big(&val.mul_i128(lab), rl);
        for c in 0..=t {
            let cmb=comb(t-cp,c);
            if cmb==0 { continue; }
            let sc=pow2((l*(t-c)) as u32).mul_i128(cmb);
            s=s.add(&mul_big(&base,&sc));
        }
    }
    // clear 144^power exactly (one factor per interior step + the read-out rho)
    for _ in 0..power { s=s.div_small_exact(DETG); }
    s
}

// =====================================================================
//  Part C.  Fused odd-cycle-sun DP:  oc_v(t) = (-1)^v * sum over labeled suns.
//  Within(sun) = 2^{1 + sum C(l_i,2)} ; labelling / cyclic seam per the NOTE.
//  State (a, b, u): a = in-copies inside seed S0, b = outside, u = labels used.
// =====================================================================
fn oc_v(v:i128, t:i128)->Big{
    let mut total=Big::zero();
    let mut g=3i128;
    while g<=v {
        let gsum=sun_cycle(v, t, g);
        // (2*gsum)/(2g) ; the global Within factor 2 cancels the dihedral 2
        total=total.add(&gsum.div_small_exact(g as u64));
        g+=2;
    }
    if v%2==0 { total } else { total.mul_i128(-1) }
}

fn sun_cycle(v:i128, t:i128, g:i128)->Big{
    let mut gsum=Big::zero();
    for c0 in 0..=t {
        let top=t-c0;
        // ---- vertex 0 : seeds (a0,b0)=(c0,0) ----
        let mut vmap:HashMap<(i128,i128,i128),Big>=HashMap::new();
        let mult=comb(t,c0);
        let mut l0=0i128;
        while l0 < v-g+2 {                  // need room for g-1 further vertices
            let u0=1+l0;
            if u0 > v-(g-1) { break; }
            let lf=pow2((l0*(t-c0)) as u32);
            let wf=pow2(c2(l0));
            let lab=v*comb(v-1,l0);
            let term=mul_big(&lf,&wf).mul_i128(lab).mul_i128(mult);
            dp_add(&mut vmap, (c0, 0, u0), term);
            l0+=1;
        }
        // ---- vertices 1..g-1 ----
        for step in 1..g {
            let last = step==g-1;
            let rem_after=g-1-step;
            let mut nv:HashMap<(i128,i128,i128),Big>=HashMap::new();
            for (&(a2,b2,u), val) in vmap.iter() {
                let lmax=v-u-1-rem_after;       // leave >=1 label per remaining vertex
                if lmax<0 { continue; }
                let avail_a=c0-a2;
                let avail_b=top-b2;
                for l in 0..=lmax {
                    let up=u+1+l;
                    let lab=(v-u)*comb(v-u-1,l);
                    let wf=pow2(c2(l));
                    let base_l=mul_big(&val.mul_i128(lab), &wf);
                    let a_hi = if last { 0 } else { avail_a };
                    let mut a=0i128;
                    while a<=a_hi {
                        let ca=comb(avail_a,a);
                        if ca==0 { a+=1; continue; }
                        let bl=base_l.mul_i128(ca);
                        for b in 0..=avail_b {
                            let cb=comb(avail_b,b);
                            if cb==0 { continue; }
                            let lf=pow2((l*(t-a-b)) as u32);
                            let term=mul_big(&bl.mul_i128(cb), &lf);
                            dp_add(&mut nv, (a, b, up), term);
                        }
                        a+=1;
                    }
                }
            }
            vmap=nv;
        }
        for (&(_,_,u), val) in vmap.iter() {
            if u==v { gsum=gsum.add(val); }
        }
    }
    gsum
}

// =====================================================================
//  Part D.  Outer engine  A_v(t) and L(n)  (the already-polynomial O(v^2) part).
// =====================================================================
fn a_forest(v:i128, t:i128)->Big{
    let wc:Vec<Big>=(0..=v).map(|k| if k>=1 { wk(k,t) } else { Big::zero() }).collect();
    let mut memo:HashMap<(i128,i128),Big>=HashMap::new();
    d_rec(v, 0, &wc, &mut memo)
}
fn d_rec(vv:i128, j:i128, wc:&[Big], memo:&mut HashMap<(i128,i128),Big>)->Big{
    if vv==0 { return Big::from_i128(1); }
    if let Some(x)=memo.get(&(vv,j)) { return x.clone(); }
    let mut s=Big::zero();
    for k in 1..=vv {
        let w=&wc[k as usize];
        if w.is_zero() { continue; }
        let sub=d_rec(vv-k, j+1, wc, memo);
        let mut term=mul_big(w, &sub);
        term=term.mul_i128(comb(vv-1, k-1));
        term=mul_big(&term, &pow2(j as u32));
        s=s.add(&term);
    }
    memo.insert((vv,j), s.clone());
    s
}
fn a_v(v:i128, t:i128)->Big{ a_forest(v,t).add(&oc_v(v,t)) }

fn l_count(n:i128)->Big{
    let mut tot=pow2(c2(n));
    for v in 2..=n {
        let t=n-v;
        let mut term=a_v(v,t).mul_i128(comb(n,v));
        term=mul_big(&term, &pow2(c2(n-v)));
        tot=tot.add(&term);
    }
    tot
}

fn main(){
    let args:Vec<String>=std::env::args().collect();
    if args.len()<2 {
        eprintln!("usage: {} N   (1 <= N <= 30; gamma=3 label exact for N<=18)", args[0]);
        std::process::exit(1);
    }
    let n:i128 = match args[1].parse() { Ok(x)=>x, Err(_)=>{ eprintln!("bad N"); std::process::exit(1); } };
    if n<1 || n>30 {
        eprintln!("N must satisfy 1 <= N <= 30 (the polynomial L(n) is computed for all such N; \
                   the gamma=3 label equals gamma>=3 only for N<=18).");
        std::process::exit(1);
    }
    let l=l_count(n);

    // domination-number distribution
    let g1=pow2(c2(n-1)).mul_i128(n);
    let g3=l.clone();
    let ttl=pow2(c2(n));
    let g2=ttl.add(&g1.mul_i128(-1)).add(&g3.mul_i128(-1));

    println!("L({}) = {}", n, l.to_string());
    println!("domination-number distribution (n={}):", n);
    println!("  gamma=1: {}", g1.to_string());
    println!("  gamma=2: {}", g2.to_string());
    println!("  gamma=3: {}{}", g3.to_string(), if n>18 {"   (= gamma>=3; may exceed gamma=3 for n>18)"} else {""});
    println!("  total  : {}", ttl.to_string());
}
