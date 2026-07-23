// Burnside count of tournaments up to isomorphism, in closed form over cycle types.
// Validates against A000568, and shows the framework the S2 version must plug into.
//
// #iso = (1/n!) sum_{lambda |- n, all parts ODD} (n!/z_lambda) * 2^{E(lambda)}
//   E(lambda) = sum_i (l_i - 1)/2  +  sum_{i<j} gcd(l_i, l_j)     (# edge-orbits of pi)
//   (any even part => a self-reversing edge-orbit => 0 fixed tournaments)

fn gcd(a: u64, b: u64) -> u64 { if b == 0 { a } else { gcd(b, a % b) } }

fn factorial(n: u64) -> u128 { (1..=n).map(|x| x as u128).product::<u128>().max(1) }

// z_lambda = prod_k (k^{m_k} * m_k!), m_k = multiplicity of part k
fn z_lambda(parts: &[u64]) -> u128 {
    let mut mult = std::collections::HashMap::new();
    for &p in parts { *mult.entry(p).or_insert(0u64) += 1; }
    let mut z: u128 = 1;
    for (&k, &m) in &mult {
        z *= (k as u128).pow(m as u32);
        z *= factorial(m);
    }
    z
}

fn edge_orbits(parts: &[u64]) -> u64 {
    let mut e = 0u64;
    for &l in parts { e += (l - 1) / 2; }
    for i in 0..parts.len() {
        for j in (i + 1)..parts.len() {
            e += gcd(parts[i], parts[j]);
        }
    }
    e
}

// generate partitions of n (as nonincreasing parts)
fn partitions(n: u64, max: u64, cur: &mut Vec<u64>, out: &mut Vec<Vec<u64>>) {
    if n == 0 { out.push(cur.clone()); return; }
    let hi = max.min(n);
    let mut p = hi;
    while p >= 1 {
        cur.push(p);
        partitions(n - p, p, cur, out);
        cur.pop();
        p -= 1;
    }
}

fn main() {
    for n in 1u64..=16 {
        let mut parts_list = Vec::new();
        partitions(n, n, &mut Vec::new(), &mut parts_list);
        let nfac = factorial(n);
        // sum (n!/z) * 2^E over all-odd partitions; then /n!
        // accumulate as fraction sum of (2^E / z); multiply by n! at end -> integer
        let mut all_terms: Vec<(u128, u64)> = Vec::new(); // (n!/z, E)
        let mut used = 0u64;
        for parts in &parts_list {
            if parts.iter().any(|&l| l % 2 == 0) { continue; }
            used += 1;
            let z = z_lambda(parts);
            let coeff = nfac / z; // n!/z_lambda is an integer
            let e = edge_orbits(parts);
            all_terms.push((coeff, e));
        }
        // total = sum coeff * 2^E ; iso = total / n!
        let mut total: u128 = 0;
        for (coeff, e) in &all_terms {
            total = total.checked_add(coeff.checked_mul(1u128 << e).expect("mul of")).expect("add of");
        }
        let iso = total / nfac;
        let rem = total % nfac;
        println!("n={:2}  iso_tournaments={:<28}  (odd-partitions used={}, exact_div={})",
                 n, iso, used, rem == 0);
    }

    // E(lambda) distribution for n=13 and n=14: which Burnside terms are cheap vs the wall.
    for n in [13u64, 14] {
        println!("\n=== E(lambda) distribution, n={} ===", n);
        let mut parts_list = Vec::new();
        partitions(n, n, &mut Vec::new(), &mut parts_list);
        let nfac = factorial(n);
        let mut rows: Vec<(u64, u128, Vec<u64>)> = Vec::new();
        for parts in &parts_list {
            if parts.iter().any(|&l| l % 2 == 0) { continue; }
            let e = edge_orbits(parts);
            let coeff = nfac / z_lambda(parts);
            rows.push((e, coeff, parts.clone()));
        }
        rows.sort();
        let mut easy = 0; // E <= 24 : brute 2^E feasible
        for (e, coeff, parts) in &rows {
            let tag = if *e <= 24 { easy += 1; "cheap" } else { "HARD " };
            println!("  E={:2} [{}]  coeff(n!/z)={:<20}  lambda={:?}", e, tag, coeff, parts);
        }
        println!("  --> {}/{} terms are cheap (2^E enumerable, E<=24); {} are the wall (E>24)",
                 easy, rows.len(), rows.len() - easy);
    }
}
