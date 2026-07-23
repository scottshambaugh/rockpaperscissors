// Within(caterpillar) via the 4-state spine TRANSFER MATRIX (exact, polynomial time).
//
// Within(C) = # tournaments on V(C) s.t. every C-edge {x,y} has no common in-neighbor
//           = # tournaments s.t. every out-neighborhood is a C-independent set.
//
// Caterpillar C = spine s_0-s_1-...-s_{m-1} (path), s_i carries l_i extra leaves.
//
// THEORY (derived + validated to zero mismatches vs brute force for all k<=8, and vs
// pruned-DFS for hundreds of caterpillars up to k~13):
//  * The counting series over spine words l_0 l_1 ... has Hankel rank exactly 4
//    (the double-star table [a,b] has rank 4, det 144). So a 4-dim linear (weighted-
//    automaton) representation exists:  Within = row(l_0)·M(l_1)···M(l_{m-2})·col(l_{m-1}).
//  * Per-spine-vertex "leaf factor" has the exact closed form
//        f(l) = 2^C(l,2)*|A|^l + sigma*l*|A*|*2^C(l-1,2)*|A|^(l-1)
//    (each out-neighborhood being C-independent forces |A| = a power of two), which
//    fixes the leaf-count dependence of every channel to combinations of
//        2^C(l,2)*A^l   and   l*2^C(l-1,2)*A^(l-1)  for A in {1,2,4,8}.
//  * We realize the automaton in the basis of "single-symbol suffixes" 0,1,2,3:
//    state R[j] = Within(prefix followed by a spine vertex carrying j leaves), j=0..3.
//    Transitions/read-out solve the 4x4 system G (G[a][d]=DS(a,d)); adj(G)=144 G^{-1}.
//
// Closed-form channels (validated):
//   star(m)   = 2^C(m,2) + m 2^C(m-1,2)                          (single spine vertex)
//   DS(L,0)   = star(L+1);   DS(L,1)=3 2^C(L,2)+L 2^C(L-1,2);
//   DS(L,2)   = 2 DS(L,1);   DS(L,3)=22 2^C(L,2)+8 L 2^C(L-1,2)  (L>=1; L=0 => star(.+1))
//   triples g_{a,c}(l)=Within([a,l,c]) : combinations in the basis above (see TRIP).
//
// Usage: within_tm L0 L1 L2 ...
//
// Arithmetic: a minimal signed big integer keeps everything exact for any k.

#[derive(Clone)]
struct Big { neg: bool, mag: Vec<u64> } // little-endian base 2^64, no trailing zeros
impl Big {
    fn zero() -> Big { Big{neg:false, mag:vec![]} }
    fn norm(mut self) -> Big { while let Some(&0)=self.mag.last(){ self.mag.pop(); } if self.mag.is_empty(){self.neg=false;} self }
    fn from_i128(x: i128) -> Big {
        let neg = x<0; let mut u = (if neg { (x as i128).unsigned_abs() } else { x as u128 }) as u128;
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
        // multiply magnitude (base 2^64) by u128: do it as two 64-bit multipliers
        let lo = (u & 0xFFFFFFFFFFFFFFFF) as u64; u>>=64; let hi = u as u64;
        // result = mag*lo + (mag*hi)<<64
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
            // shift left by 64 bits (one limb)
            p_hi.insert(0,0);
            res = res.add(&Big{neg:false, mag:p_hi});
        }
        res.neg=neg; res.norm()
    }
    fn div_small_exact(&self, d:u64)->Big{ // exact division, d>0
        // long division base 2^64
        let mut rem=0u128; let mut out=vec![0u64;self.mag.len()];
        for i in (0..self.mag.len()).rev(){
            let cur=(rem<<64) | (self.mag[i] as u128);
            out[i]=(cur/(d as u128)) as u64; rem=cur%(d as u128);
        }
        // rem should be 0
        Big{neg:self.neg, mag:out}.norm()
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

fn c2(n:i128)->u32{ (n*(n-1)/2) as u32 } // C(n,2) for n>=0
fn pow2(e:u32)->Big{ // 2^e as Big
    let mut m=vec![0u64; (e/64+1) as usize]; m[(e/64) as usize]=1u64<<(e%64); Big{neg:false,mag:m}.norm()
}
// 2^C(L,2) and L*2^C(L-1,2) as Big
fn f_p(l:i128)->Big{ pow2(c2(l)) }               // 2^C(l,2)
fn f_q(l:i128)->Big{ if l<=0 { Big::zero() } else { pow2(c2(l-1)).mul_i128(l) } } // l*2^C(l-1,2)

// basis eval: coeffs over [p1,p2,p4,p8,q1,q2,q4,q8]
fn chan(l:i128, c:[i128;8])->Big{
    let As=[1i128,2,4,8];
    let mut acc=Big::zero();
    for k in 0..4 {
        // p term: c[k]*2^C(l,2)*A^l
        if c[k]!=0 {
            let mut t=f_p(l); // 2^C(l,2)
            // multiply by A^l
            let mut al=1i128; for _ in 0..l { al*=As[k]; } // A^l (l small<=~13, A<=8 -> fits i128 up to 8^13=2^39)
            t=t.mul_i128(al).mul_i128(c[k]);
            acc=acc.add(&t);
        }
        if c[4+k]!=0 {
            let mut t=f_q(l); // l*2^C(l-1,2)
            let mut al=1i128; for _ in 0..(l-1).max(0) { al*=As[k]; } // A^(l-1)
            if l>=1 { t=t.mul_i128(al).mul_i128(c[4+k]); acc=acc.add(&t); }
        }
    }
    acc
}

fn star(m:i128)->Big{ // 2^C(m,2)+m 2^C(m-1,2)
    f_p(m).add(&f_q(m))
}
// double-star channels d_j(L)=DS(L,j), asymptotic L>=1
const D1:[i128;8]=[3,0,0,0, 1,0,0,0];
const D2:[i128;8]=[6,0,0,0, 2,0,0,0];
const D3:[i128;8]=[22,0,0,0, 8,0,0,0];
fn ds(p:i128,q:i128)->Big{
    if p==0 { return star(q+1); }
    if q==0 { return star(p+1); }
    // both>=1
    let (a,b)= if q<4 {(p,q)} else {(q,p)}; // vary a with fixed small index b
    match b { 1=>chan(a,D1),2=>chan(a,D2),3=>chan(a,D3), _=>unreachable!() }
}
// triple channels g_{a,c}(l), a<=c, over basis
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
// adj(G) = 144 * G^{-1}
const ADJ:[[i128;4];4]=[
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
    // R[j]=DS(seq[0], j)
    let mut r:[Big;4]=[ds(seq[0],0),ds(seq[0],1),ds(seq[0],2),ds(seq[0],3)];
    // transitions for middle symbols
    for &l in &seq[1..m-1] {
        // S[a] = sum_d adj[d][a] * R[d]
        let mut s:[Big;4]=[Big::zero(),Big::zero(),Big::zero(),Big::zero()];
        for a in 0..4 { let mut acc=Big::zero(); for d in 0..4 { acc=acc.add(&r[d].mul_i128(ADJ[d][a])); } s[a]=acc; }
        // R'[j] = (sum_a w_j[a]*S[a]) / 144 ; w_j[a]=triple(a,l,j)
        let mut nr:[Big;4]=[Big::zero(),Big::zero(),Big::zero(),Big::zero()];
        for j in 0..4 {
            let mut acc=Big::zero();
            for a in 0..4 { let w=triple(a as i128,l,j as i128); acc=acc.add(&mul_big(&w,&s[a])); }
            nr[j]=acc.div_small_exact(DETG);
        }
        r=nr;
    }
    // read-out with last symbol x
    let x=seq[m-1];
    // lam = G^{-1} * [DS(a,x)]_a ; Within = sum_d lam[d] R[d] = (1/144) sum_d (sum_a adj[d][a] DS(a,x)) R[d]
    let mut acc=Big::zero();
    for d in 0..4 {
        let mut coeff=Big::zero();
        for a in 0..4 { coeff=coeff.add(&ds(a as i128,x).mul_i128(ADJ[d][a])); }
        acc=acc.add(&mul_big(&coeff,&r[d]));
    }
    acc.div_small_exact(DETG)
}
// big*big multiply (needed for w*S). Implement via repeated mul_i128 over limbs of the smaller.
fn mul_big(a:&Big,b:&Big)->Big{
    if a.is_zero()||b.is_zero(){ return Big::zero(); }
    // multiply a by each limb of b
    let mut res=Big::zero();
    for (i,&limb) in b.mag.iter().enumerate(){
        if limb==0 { continue; }
        let mut part=a.mul_u64(limb);
        for _ in 0..i { part.mag.insert(0,0); }
        res=res.add(&part);
    }
    res.neg = a.neg ^ b.neg; res.norm()
}
impl Big{
    fn mul_u64(&self,f:u64)->Big{
        if f==0||self.is_zero(){return Big::zero();}
        let mut r=vec![]; let mut carry=0u128;
        for &mi in &self.mag { let p=(mi as u128)*(f as u128)+carry; r.push((p&0xFFFFFFFFFFFFFFFF) as u64); carry=p>>64; }
        if carry>0 { r.push(carry as u64); }
        Big{neg:self.neg,mag:r}.norm()
    }
}

fn main(){
    let seq:Vec<i128>=std::env::args().skip(1).map(|s|s.parse().unwrap()).collect();
    let k:i128 = seq.len() as i128 + seq.iter().sum::<i128>();
    let w=within(&seq);
    let ds:Vec<String>=seq.iter().map(|x|x.to_string()).collect();
    println!("[{}] k={} Within={}", ds.join(","), k, w.to_string());
}
