// two_paradox_count.rs -- LABELED tournaments with domination number gamma >= 3.
//
// Computes  L(n) = # labeled tournaments on n vertices with no dominating pair
//                = # in which every pair of vertices has a common dominator
//                = # with domination number gamma >= 3.
// For n <= 18 every tournament has gamma <= 3 (the smallest tournament with
// gamma = 4 has 19 vertices), so gamma >= 3 == gamma = 3 there, and this tool
// also prints the full domination-number distribution:
//     gamma = 1 : n * 2^C(n-1,2)                       (a single dominating vertex)
//     gamma = 3 : L(n)
//     gamma = 2 : 2^C(n,2) - (gamma=1) - (gamma=3)      (by complement)
// The distribution labels are valid only for n <= 18; the tool is capped there.
//
// METHOD (a direct port of the validated Python pipeline pipe.py):
//   Master inclusion-exclusion over the "common-dominator" hypergraph:
//     L(n) = sum_{v=0..n} C(n,v) * 2^C(n-v,2) * A_v(n-v),
//     A_v(t) = sum_{H no-isolated on [v]} (-1)^{e(H)} * Within(H) * I(H)^t.
//   The domination-graph structure theorem restricts the H with Within(H) != 0
//   to (a) caterpillar-forests and (b) single connected odd-cycle-suns, so:
//     * caterpillar-forest sparsity: an O(v^3) forest recurrence
//         D(v,j) = sum_{k=1..v} C(v-1,k-1) * w_k(t) * 2^j * D(v-k,j+1),  D(0,.)=1
//       where the 2^j factor is the 2^C(c,2) inter-component coupling and
//       w_k(t) = sum over connected caterpillar iso-classes on k vertices of
//                (-1)^{k-1} (k!/|Aut|) Within(cat) I(cat)^t;
//     * Within(caterpillar) via the exact 4-state spine transfer matrix
//       (matrices G, adj(G)=144 G^{-1}, W(l); reused verbatim from
//       rust/caterpillar_within.rs, validated to zero mismatches vs brute);
//     * odd-cycle term: single odd-cycle-suns enumerated up to the dihedral
//       group D_g, with I = independent-set count, |Aut| = |Stab_{D_g}| * prod l_i!,
//       and Within = 2^(1 + sum C(l_i,2)).
//
// All arithmetic is exact via a minimal signed big integer (from
// caterpillar_within.rs); transfer-matrix intermediates exceed i128 even when
// the answer fits, so bignum is required. Range: 1 <= n <= 18.
//
// Reproduces L(7)=240, L(8)=109440, L(9)=77274880, L(10)=107731095040,
// L(11)=277644146489600, L(12)=1259251974345052160,
// L(13)=10059739307720354796544, L(14)=144025884971245392244400128,
// L(15)=3764299177453034811662226208768,
// L(16)=182456135219244423561589552542810112,
// L(17)=16611757357122373923247804720005425659904,
// L(18)=2870353532067896148226715034553628871332397056.
//
// Usage: two_paradox_count N        (1 <= N <= 18)
// Build:  rustc -O -C target-cpu=native -C overflow-checks=on \
//             rust/two_paradox_count.rs -o /tmp/tpc

use std::collections::{HashMap, HashSet};

// ======================= minimal signed big integer =======================
// (verbatim from rust/caterpillar_within.rs -- the Within transfer matrix needs
//  it, and reusing it keeps the whole pipeline exact.)
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
// (verbatim from rust/caterpillar_within.rs; leaf-count-word spec.)
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
fn star(m:i128)->Big{ f_p(m).add(&f_q(m)) } // 2^C(m,2)+m 2^C(m-1,2)
const D1:[i128;8]=[3,0,0,0, 1,0,0,0];
const D2:[i128;8]=[6,0,0,0, 2,0,0,0];
const D3:[i128;8]=[22,0,0,0, 8,0,0,0];
fn ds(p:i128,q:i128)->Big{
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
fn triple(a:i128,l:i128,c:i128)->Big{
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
fn within(seq:&[i128])->Big{
    let m=seq.len();
    if m==0 { return Big::from_i128(1); }
    if m==1 { return star(seq[0]); }
    let mut r:[Big;4]=[ds(seq[0],0),ds(seq[0],1),ds(seq[0],2),ds(seq[0],3)];
    for &l in &seq[1..m-1] {
        let mut s:[Big;4]=[Big::zero(),Big::zero(),Big::zero(),Big::zero()];
        for a in 0..4 { let mut acc=Big::zero(); for d in 0..4 { acc=acc.add(&r[d].mul_i128(ADJ[d][a])); } s[a]=acc; }
        let mut nr:[Big;4]=[Big::zero(),Big::zero(),Big::zero(),Big::zero()];
        for j in 0..4 {
            let mut acc=Big::zero();
            for a in 0..4 { let w=triple(a as i128,l,j as i128); acc=acc.add(&mul_big(&w,&s[a])); }
            nr[j]=acc.div_small_exact(DETG);
        }
        r=nr;
    }
    let x=seq[m-1];
    let mut acc=Big::zero();
    for d in 0..4 {
        let mut coeff=Big::zero();
        for a in 0..4 { coeff=coeff.add(&ds(a as i128,x).mul_i128(ADJ[d][a])); }
        acc=acc.add(&mul_big(&coeff,&r[d]));
    }
    acc.div_small_exact(DETG)
}

// ============================ small integer helpers ============================
fn comb(n:i128,k:i128)->i128{
    if k<0 || k>n { return 0; }
    let k = k.min(n-k);
    let mut r:i128=1;
    for i in 0..k { r = r*(n-i)/(i+1); }
    r
}
fn factorial(k:i128)->i128{ let mut r=1i128; for i in 2..=k { r*=i; } r }
fn ipow(base:i128, t:usize)->Big{
    let mut r=Big::from_i128(1);
    for _ in 0..t { r=r.mul_i128(base); }
    r
}

// ============================ independent-set count ============================
fn indep_rec(avail:u64, adjm:&[u64], memo:&mut HashMap<u64,i128>)->i128{
    if avail==0 { return 1; }
    if let Some(&x)=memo.get(&avail){ return x; }
    let x = avail.trailing_zeros() as usize;
    let without = avail & !(1u64<<x);
    let r = indep_rec(without, adjm, memo) + indep_rec(without & !adjm[x], adjm, memo);
    memo.insert(avail, r);
    r
}
fn indep_count(nv:usize, adjm:&[u64])->i128{
    let mut memo=HashMap::new();
    let full = if nv==64 { u64::MAX } else { (1u64<<nv)-1 };
    indep_rec(full, adjm, &mut memo)
}

// ============================ caterpillar enumeration ============================
fn build_caterpillar(seq:&[usize])->(usize, Vec<Vec<usize>>){
    let m=seq.len();
    let mut edges:Vec<(usize,usize)>=vec![];
    for i in 0..m-1 { edges.push((i,i+1)); }
    let mut nxt=m;
    for i in 0..m { for _ in 0..seq[i] { edges.push((i,nxt)); nxt+=1; } }
    let k=nxt;
    let mut adj=vec![vec![];k];
    for (a,b) in edges { adj[a].push(b); adj[b].push(a); }
    (k,adj)
}
fn adj_bitmask(k:usize, adj:&[Vec<usize>])->Vec<u64>{
    let mut m=vec![0u64;k];
    for v in 0..k { for &u in &adj[v] { m[v]|=1u64<<u; } }
    m
}
fn tree_centers(k:usize, adj:&[Vec<usize>])->Vec<usize>{
    if k==1 { return vec![0]; }
    let mut d:Vec<usize>=(0..k).map(|v|adj[v].len()).collect();
    let mut removed=vec![false;k];
    let mut cnt=k;
    let mut layer:Vec<usize>=(0..k).filter(|&v|d[v]==1).collect();
    while cnt>2 {
        let mut nxt=vec![];
        for &v in &layer {
            removed[v]=true; cnt-=1;
            for &u in &adj[v] {
                if !removed[u] {
                    d[u]-=1;
                    if d[u]==1 { nxt.push(u); }
                }
            }
        }
        layer=nxt;
    }
    (0..k).filter(|&v|!removed[v]).collect()
}
// canonical string + automorphism count of the subtree rooted at v (parent excluded)
fn rooted(v:usize, parent:isize, adj:&[Vec<usize>])->(String,i128){
    let mut kids:Vec<(String,i128)>=vec![];
    for &c in &adj[v] { if c as isize != parent { kids.push(rooted(c, v as isize, adj)); } }
    kids.sort_by(|a,b| a.0.cmp(&b.0));
    let mut aut:i128=1;
    for r in &kids { aut*=r.1; }
    // multiply by (multiplicity!) over groups of equal canonical strings
    let mut i=0;
    while i<kids.len(){
        let mut j=i+1;
        while j<kids.len() && kids[j].0==kids[i].0 { j+=1; }
        aut *= factorial((j-i) as i128);
        i=j;
    }
    let mut canon=String::from("(");
    for r in &kids { canon.push_str(&r.0); }
    canon.push(')');
    (canon, aut)
}
fn tree_canon_aut(k:usize, adj:&[Vec<usize>])->(String,i128){
    let ctr=tree_centers(k,adj);
    if ctr.len()==1 {
        rooted(ctr[0], -1, adj)
    } else {
        let (c1,c2v)=(ctr[0],ctr[1]);
        let r1=rooted(c1, c2v as isize, adj);
        let r2=rooted(c2v, c1 as isize, adj);
        let (lo,hi)= if r1.0<=r2.0 {(r1.0.clone(),r2.0.clone())} else {(r2.0.clone(),r1.0.clone())};
        let aut=r1.1*r2.1*(if r1.0==r2.0 {2} else {1});
        (format!("BI({}|{})",lo,hi), aut)
    }
}

// ============================ compositions ============================
fn compose(parts:usize, total:usize, cur:&mut Vec<usize>, out:&mut Vec<Vec<usize>>){
    if parts==0 { if total==0 { out.push(cur.clone()); } return; }
    for first in 0..=total { cur.push(first); compose(parts-1, total-first, cur, out); cur.pop(); }
}
fn compositions(total:usize, parts:usize)->Vec<Vec<usize>>{
    let mut out=vec![]; let mut cur=vec![]; compose(parts,total,&mut cur,&mut out); out
}
fn gen_seqs(k:usize)->Vec<Vec<usize>>{
    let mut out=vec![];
    for m in 1..=k { for c in compositions(k-m, m) { out.push(c); } }
    out
}

// ============================ odd-cycle-sun helpers ============================
fn dihedral_perm(g:usize, rot:bool, r:usize, i:usize)->usize{
    if rot { (i+r)%g } else { ((r as i64 - i as i64).rem_euclid(g as i64)) as usize }
}
fn dihedral_canon(seq:&[usize])->Vec<usize>{
    let g=seq.len();
    let mut best:Option<Vec<usize>>=None;
    for &rot in &[true,false] {
        for r in 0..g {
            let cand:Vec<usize>=(0..g).map(|i|seq[dihedral_perm(g,rot,r,i)]).collect();
            if best.as_ref().map_or(true,|b|cand<*b) { best=Some(cand); }
        }
    }
    best.unwrap()
}
fn aut_oc(g:usize, leaves:&[usize])->i128{
    let mut stab=0i128;
    for &rot in &[true,false] {
        for r in 0..g {
            if (0..g).all(|i| leaves[dihedral_perm(g,rot,r,i)]==leaves[i]) { stab+=1; }
        }
    }
    let mut prod=1i128; for &l in leaves { prod*=factorial(l as i128); }
    stab*prod
}
fn build_sun(g:usize, leaves:&[usize])->(usize, Vec<u64>){
    let mut edges:Vec<(usize,usize)>=(0..g).map(|i|(i,(i+1)%g)).collect();
    let mut nxt=g;
    for i in 0..g { for _ in 0..leaves[i] { edges.push((i,nxt)); nxt+=1; } }
    let nn=nxt;
    let mut adjm=vec![0u64;nn];
    for (a,b) in edges { adjm[a]|=1u64<<b; adjm[b]|=1u64<<a; }
    (nn, adjm)
}

// ============================ coefficient tables ============================
// cat_coeff[k][I] = sum over connected caterpillar iso-classes on k vertices of
//                   (-1)^{k-1} (k!/|Aut|) Within(cat), grouped by I(cat).
fn build_cat_coeff(nmax:usize)->Vec<HashMap<i128,Big>>{
    let mut tab:Vec<HashMap<i128,Big>>=(0..=nmax).map(|_|HashMap::new()).collect();
    for k in 2..=nmax {
        // dedup caterpillar shapes by tree canonical form; first seq wins.
        let mut seen:HashMap<String,(Vec<usize>,i128,usize,Vec<u64>)>=HashMap::new();
        for seq in gen_seqs(k) {
            let (kk,adj)=build_caterpillar(&seq);
            let (canon,aut)=tree_canon_aut(kk,&adj);
            if seen.contains_key(&canon) { continue; }
            let bm=adj_bitmask(kk,&adj);
            seen.insert(canon,(seq,aut,kk,bm));
        }
        let kf=factorial(k as i128);
        let sign:i128 = if (k-1)%2==0 {1} else {-1};
        for (_c,(seq,aut,kk,bm)) in &seen {
            let ic=indep_count(*kk,bm);
            let seqi:Vec<i128>=seq.iter().map(|&x|x as i128).collect();
            let w=within(&seqi);
            let mult=kf/aut;
            let term=w.mul_i128(sign*mult);
            let e=tab[k].entry(ic).or_insert_with(Big::zero);
            *e=e.add(&term);
        }
    }
    tab
}
// oc_coeff[v][I] = sum over odd-cycle-sun iso-classes on v vertices of
//                  (-1)^{v} (v!/|Aut|) Within, grouped by I.
fn build_oc_coeff(nmax:usize)->Vec<HashMap<i128,Big>>{
    let mut tab:Vec<HashMap<i128,Big>>=(0..=nmax).map(|_|HashMap::new()).collect();
    for v in 3..=nmax {
        let vf=factorial(v as i128);
        let sign:i128 = if v%2==0 {1} else {-1};
        let mut g=3;
        while g<=v {
            let l=v-g;
            let mut seeng:HashSet<Vec<usize>>=HashSet::new();
            for comp in compositions(l,g) {
                let dc=dihedral_canon(&comp);
                if !seeng.insert(dc.clone()) { continue; }
                let (nn,adjm)=build_sun(g,&comp);
                let ic=indep_count(nn,&adjm);
                let aut=aut_oc(g,&comp);
                let mut e=1i128; for &x in &comp { e += (x as i128)*(x as i128 -1)/2; }
                let w=pow2(e as u32);
                let mult=vf/aut;
                let term=w.mul_i128(sign*mult);
                let ent=tab[v].entry(ic).or_insert_with(Big::zero);
                *ent=ent.add(&term);
            }
            g+=2;
        }
    }
    tab
}

// ============================ recurrence & L(n) ============================
fn poly_eval(coeff:&HashMap<i128,Big>, t:usize)->Big{
    let mut s=Big::zero();
    for (i,c) in coeff { s=s.add(&mul_big(c,&ipow(*i,t))); }
    s
}
fn d_rec(vv:usize, j:usize, wcache:&[Big], memo:&mut HashMap<(usize,usize),Big>)->Big{
    if vv==0 { return Big::from_i128(1); }
    if let Some(x)=memo.get(&(vv,j)){ return x.clone(); }
    let mut s=Big::zero();
    for k in 1..=vv {
        if wcache[k].is_zero() { continue; }
        let sub=d_rec(vv-k, j+1, wcache, memo);
        let mut term=mul_big(&wcache[k], &sub);
        term=term.mul_i128(comb((vv-1) as i128,(k-1) as i128));
        term=term.mul_i128(1i128<<j);
        s=s.add(&term);
    }
    memo.insert((vv,j), s.clone());
    s
}
fn a_tree(v:usize, t:usize, cat_coeff:&[HashMap<i128,Big>])->Big{
    let mut wcache=vec![Big::zero(); v+1];
    for k in 1..=v { wcache[k]=poly_eval(&cat_coeff[k], t); }
    let mut memo=HashMap::new();
    d_rec(v,0,&wcache,&mut memo)
}
fn l_count(n:usize, cat_coeff:&[HashMap<i128,Big>], oc_coeff:&[HashMap<i128,Big>])->Big{
    let mut total=pow2(c2(n as i128));
    for v in 2..=n {
        let t=n-v;
        let av=a_tree(v,t,cat_coeff).add(&poly_eval(&oc_coeff[v], t));
        let mut term=av.mul_i128(comb(n as i128, v as i128));
        term=mul_big(&term, &pow2(c2((n-v) as i128)));
        total=total.add(&term);
    }
    total
}

fn main(){
    let args:Vec<String>=std::env::args().collect();
    if args.len()<2 {
        eprintln!("usage: {} N   (1 <= N <= 18)", args[0]);
        std::process::exit(1);
    }
    let n:usize = match args[1].parse() { Ok(x)=>x, Err(_)=>{ eprintln!("bad N"); std::process::exit(1); } };
    if n<1 || n>18 {
        eprintln!("N must satisfy 1 <= N <= 18 (beyond n=18 the gamma=3 label is \
                   invalid: the smallest tournament with domination number 4 has 19 vertices).");
        std::process::exit(1);
    }
    let cat_coeff=build_cat_coeff(n);
    let oc_coeff=build_oc_coeff(n);
    let l=l_count(n,&cat_coeff,&oc_coeff);

    // domination-number distribution
    let g1=pow2(c2((n as i128)-1)).mul_i128(n as i128);
    let g3=l.clone();
    let ttl=pow2(c2(n as i128));
    let g2=ttl.add(&g1.mul_i128(-1)).add(&g3.mul_i128(-1));

    println!("L({}) = {}", n, l.to_string());
    println!("domination-number distribution (n={}):", n);
    println!("  gamma=1: {}", g1.to_string());
    println!("  gamma=2: {}", g2.to_string());
    println!("  gamma=3: {}", g3.to_string());
    println!("  total  : {}", ttl.to_string());
}
