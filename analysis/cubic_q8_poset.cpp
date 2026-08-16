// Exhaustive affine-orbit classification and processor search for cubic maps
// F_8^2 -> F_8.  F_8 = F_2[a]/(a^3+a+1), with elements bit-packed 0..7.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using F = uint8_t;

static F mul(F a, F b) {
  F out = 0;
  while (b) {
    if (b & 1) out ^= a;
    b >>= 1;
    bool high = a & 4;
    a <<= 1;
    if (high) a ^= 0b1011;
  }
  return out & 7;
}

static F fpow(F a, int n) {
  F out = 1;
  while (n) {
    if (n & 1) out = mul(out, a);
    a = mul(a, a);
    n >>= 1;
  }
  return out;
}

static F inv(F a) {
  assert(a != 0);
  return fpow(a, 6);
}

struct Poly {
  // Coefficient of x^i y^j lives at c[index(i,j)].
  std::array<F, 28> c{};
};

static int pindex[7][7];
static std::array<std::pair<int, int>, 28> pexp;

static void init_monomials() {
  int at = 0;
  for (int degree = 0; degree <= 6; ++degree) {
    for (int j = 0; j <= degree; ++j) {
      int i = degree - j;
      pindex[i][j] = at;
      pexp[at++] = {i, j};
    }
  }
  assert(at == 28);
}

static Poly padd(const Poly& a, const Poly& b) {
  Poly out;
  for (int i = 0; i < 28; ++i) out.c[i] = a.c[i] ^ b.c[i];
  return out;
}

static Poly pscale(const Poly& a, F s) {
  Poly out;
  for (int i = 0; i < 28; ++i) out.c[i] = mul(a.c[i], s);
  return out;
}

static Poly pmul(const Poly& a, const Poly& b) {
  Poly out;
  for (int ai = 0; ai < 28; ++ai) {
    if (!a.c[ai]) continue;
    auto [ax, ay] = pexp[ai];
    for (int bi = 0; bi < 28; ++bi) {
      if (!b.c[bi]) continue;
      auto [bx, by] = pexp[bi];
      if (ax + ay + bx + by <= 6) {
        int oi = pindex[ax + bx][ay + by];
        out.c[oi] ^= mul(a.c[ai], b.c[bi]);
      }
    }
  }
  return out;
}

static Poly pone() {
  Poly out;
  out.c[pindex[0][0]] = 1;
  return out;
}

static Poly ppow(const Poly& a, int n) {
  Poly out = pone();
  for (int i = 0; i < n; ++i) out = pmul(out, a);
  return out;
}

static const std::array<std::pair<int, int>, 10> BASIS3 = {{
    {0,0}, {1,0}, {0,1}, {2,0}, {1,1}, {0,2},
    {3,0}, {2,1}, {1,2}, {0,3}
}};

static Poly decode10(uint32_t code) {
  Poly out;
  for (auto [i,j] : BASIS3) {
    out.c[pindex[i][j]] = code & 7;
    code >>= 3;
  }
  return out;
}

static uint32_t lower_code(const Poly& p) {
  uint32_t code = 0;
  int shift = 0;
  for (int k = 0; k < 6; ++k) {
    auto [i,j] = BASIS3[k];
    code |= uint32_t(p.c[pindex[i][j]]) << shift;
    shift += 3;
  }
  return code;
}

static uint16_t top_code(const Poly& p) {
  uint16_t code = 0;
  int shift = 0;
  for (int k = 6; k < 10; ++k) {
    auto [i,j] = BASIS3[k];
    code |= uint16_t(p.c[pindex[i][j]]) << shift;
    shift += 3;
  }
  return code;
}

static Poly from_parts(uint16_t top, uint32_t lower) {
  uint32_t code = lower | (uint32_t(top) << 18);
  return decode10(code);
}

struct LinearScale {
  F a, b, c, d, scale;
};

static std::vector<LinearScale> all_linear_scales() {
  std::vector<LinearScale> out;
  for (F a=0; a<8; ++a) for (F b=0; b<8; ++b)
  for (F c=0; c<8; ++c) for (F d=0; d<8; ++d) {
    F det = mul(a,d) ^ mul(b,c);
    if (!det) continue;
    for (F s=1; s<8; ++s) out.push_back({a,b,c,d,s});
  }
  assert(out.size() == 24696);
  return out;
}

static LinearScale inverse_ls(const LinearScale& g) {
  F deti = inv(mul(g.a,g.d) ^ mul(g.b,g.c));
  F si = inv(g.scale);
  return {mul(g.d,deti), mul(g.b,deti), mul(g.c,deti),
          mul(g.a,deti), si};
}

static Poly affine_transform(const Poly& f, const LinearScale& g,
                             F tx=0, F ty=0, F shift=0) {
  Poly x, y;
  x.c[pindex[0][0]] = tx;
  x.c[pindex[1][0]] = g.a;
  x.c[pindex[0][1]] = g.b;
  y.c[pindex[0][0]] = ty;
  y.c[pindex[1][0]] = g.c;
  y.c[pindex[0][1]] = g.d;
  std::array<Poly,4> xp = {pone(), x, ppow(x,2), ppow(x,3)};
  std::array<Poly,4> yp = {pone(), y, ppow(y,2), ppow(y,3)};
  Poly out;
  for (auto [i,j] : BASIS3) {
    F coefficient = f.c[pindex[i][j]];
    if (coefficient) out = padd(out, pscale(pmul(xp[i],yp[j]), coefficient));
  }
  out = pscale(out, g.scale);
  out.c[pindex[0][0]] ^= shift;
  return out;
}

static int image_size(const Poly& f) {
  bool seen[8]{};
  int count = 0;
  for (F x=0; x<8; ++x) for (F y=0; y<8; ++y) {
    F value = 0;
    for (auto [i,j] : BASIS3) {
      F coefficient = f.c[pindex[i][j]];
      if (coefficient) value ^= mul(coefficient, mul(fpow(x,i), fpow(y,j)));
    }
    if (!seen[value]) { seen[value] = true; ++count; }
  }
  return count;
}

struct Orbit {
  int id;
  int top_type;
  uint16_t normal_top;
  uint32_t representative_lower;
  uint32_t affine_size;
  int image;
};

struct TopData {
  uint16_t normal;
  std::vector<LinearScale> stabilizer;
  std::vector<int> lower_to_orbit;
  uint32_t top_orbit_size = 0;
};

struct AlignData {
  int top_type = -1;
  LinearScale align{};
};

class Classifier {
 public:
  std::vector<LinearScale> group;
  std::array<TopData,6> tops;
  std::array<AlignData,4096> align;
  std::vector<Orbit> orbits;

  Classifier() {
    group = all_linear_scales();
    const std::array<std::array<F,4>,6> normals = {{
      {{0,0,0,0}}, {{0,0,0,1}}, {{0,0,1,0}},
      {{0,1,1,0}}, {{0,1,1,1}}, {{1,0,1,2}}
    }};
    for (int t=0; t<6; ++t) {
      uint16_t code = 0;
      for (int k=0; k<4; ++k) code |= uint16_t(normals[t][k]) << (3*k);
      tops[t].normal = code;
    }
    build_top_data();
    build_lower_orbits();
  }

  int classify(const Poly& f) const {
    const auto& a = align[top_code(f)];
    assert(a.top_type >= 0);
    Poly normalized = affine_transform(f, a.align);
    assert(top_code(normalized) == tops[a.top_type].normal);
    int id = tops[a.top_type].lower_to_orbit[lower_code(normalized)];
    assert(id >= 0);
    return id;
  }

 private:
  void build_top_data() {
    // Zero leading part.
    align[0] = {0, {1,0,0,1,1}};
    tops[0].top_orbit_size = 1;
    tops[0].stabilizer = group;
    for (int t=1; t<6; ++t) {
      Poly normal = from_parts(tops[t].normal, 0);
      std::unordered_set<uint16_t> orbit_tops;
      for (const auto& g : group) {
        uint16_t transformed = top_code(affine_transform(normal, g));
        orbit_tops.insert(transformed);
        if (transformed == tops[t].normal) tops[t].stabilizer.push_back(g);
        if (align[transformed].top_type < 0) {
          align[transformed] = {t, inverse_ls(g)};
        }
      }
      tops[t].top_orbit_size = orbit_tops.size();
    }
    for (int code=0; code<4096; ++code) assert(align[code].top_type >= 0);
  }

  void build_lower_orbits() {
    for (int t=0; t<6; ++t) {
      auto& td = tops[t];
      td.lower_to_orbit.assign(1u<<18, -1);
      for (uint32_t lower=0; lower<(1u<<18); ++lower) {
        if (td.lower_to_orbit[lower] >= 0) continue;
        int id = orbits.size();
        std::unordered_set<uint32_t> members;
        Poly representative = from_parts(td.normal, lower);
        for (const auto& g : td.stabilizer) {
          for (F tx=0; tx<8; ++tx) for (F ty=0; ty<8; ++ty) {
            Poly base = affine_transform(representative, g, tx, ty, 0);
            for (F shift=0; shift<8; ++shift) {
              Poly moved = base;
              moved.c[pindex[0][0]] ^= shift;
              assert(top_code(moved) == td.normal);
              members.insert(lower_code(moved));
            }
          }
        }
        for (uint32_t member : members) td.lower_to_orbit[member] = id;
        uint64_t full_size = uint64_t(td.top_orbit_size) * members.size();
        assert(full_size <= UINT32_MAX);
        orbits.push_back({id, t, td.normal, lower,
                          uint32_t(full_size), image_size(representative)});
      }
      std::cerr << "top " << t << ": stabilizer=" << td.stabilizer.size()
                << " top-orbit=" << td.top_orbit_size << " total-orbits="
                << orbits.size() << "\n";
    }
  }
};

static Poly input_poly(uint16_t code) {
  // Six coefficients in order 1,x,y,x^2,xy,y^2.
  Poly p;
  for (int k=0; k<6; ++k) {
    auto [i,j] = BASIS3[k];
    p.c[pindex[i][j]] = code & 7;
    code >>= 3;
  }
  return p;
}

static Poly compose(const Poly& f, const Poly& x, const Poly& y) {
  std::array<Poly,4> xp = {pone(), x, ppow(x,2), ppow(x,3)};
  std::array<Poly,4> yp = {pone(), y, ppow(y,2), ppow(y,3)};
  Poly out;
  for (auto [i,j] : BASIS3) {
    F coefficient = f.c[pindex[i][j]];
    if (coefficient) out = padd(out, pscale(pmul(xp[i],yp[j]), coefficient));
  }
  return out;
}

static bool degree_at_most_3(const Poly& f) {
  for (int k=0; k<28; ++k) {
    auto [i,j] = pexp[k];
    if (i+j > 3 && f.c[k]) return false;
  }
  return true;
}

static Poly affine_scalar(F constant, F x_coefficient=0, F y_coefficient=0) {
  Poly out;
  out.c[pindex[0][0]] = constant;
  out.c[pindex[1][0]] = x_coefficient;
  out.c[pindex[0][1]] = y_coefficient;
  return out;
}

static void add_target(const Classifier& classifier, const Poly& source,
                       const Poly& p, const Poly& q, std::set<int>& targets) {
  Poly target = compose(source, p, q);
  if (degree_at_most_3(target)) targets.insert(classifier.classify(target));
}

static void add_affine_degenerations(const Classifier& classifier,
                                     const Poly& source,
                                     std::set<int>& targets) {
  // Rank two gives the source orbit.  The loops cover every parameterized
  // affine line (with duplication), while rank zero gives a constant.
  targets.insert(classifier.classify(source));
  targets.insert(0);
  for (F p0=0; p0<8; ++p0) for (F q0=0; q0<8; ++q0)
  for (F pu=0; pu<8; ++pu) for (F qu=0; qu<8; ++qu) {
    if (!pu && !qu) continue;
    Poly p = affine_scalar(p0, pu, 0);
    Poly q = affine_scalar(q0, qu, 0);
    add_target(classifier, source, p, q, targets);
  }
}

enum class Branch { Q_FIXED, P_FIXED, SUM_FIXED };

static void add_quadratic_branch(const Classifier& classifier,
                                 const Poly& source, Branch branch,
                                 const Poly& fixed,
                                 std::set<int>& targets) {
  for (uint32_t code=0; code<(1u<<18); ++code) {
    Poly other = input_poly(code);
    Poly p, q;
    if (branch == Branch::Q_FIXED) {
      p = other; q = fixed;
    } else if (branch == Branch::P_FIXED) {
      p = fixed; q = other;
    } else {
      p = other; q = other;
      q.c[pindex[0][0]] ^= fixed.c[pindex[0][0]];
    }
    add_target(classifier, source, p, q, targets);
  }
}

static std::vector<std::set<int>> exact_transition_search(
    const Classifier& classifier) {
  std::vector<std::set<int>> adjacency(classifier.orbits.size());
  for (const auto& orbit : classifier.orbits) {
    Poly source = from_parts(orbit.normal_top, orbit.representative_lower);
    auto& targets = adjacency[orbit.id];
    add_affine_degenerations(classifier, source, targets);

    if (orbit.top_type == 0) {
      // Orbit 1 is linear, orbit 4 parabolic, and orbit 5 split.  The other
      // quadratic types cannot acquire a cubic leading part under an allowed
      // quadratic substitution.
      if (orbit.id == 1) {
        for (int target=0; target<=6; ++target) targets.insert(target);
      } else if (orbit.id == 4) {
        for (int target=0; target<=6; ++target) targets.insert(target);
      } else if (orbit.id == 5) {
        for (F c=0; c<8; ++c) {
          add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                               affine_scalar(c), targets);
          add_quadratic_branch(classifier, source, Branch::P_FIXED,
                               affine_scalar(c), targets);
        }
        add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                             affine_scalar(0,1,0), targets);
        add_quadratic_branch(classifier, source, Branch::P_FIXED,
                             affine_scalar(0,1,0), targets);
      }
    } else if (orbit.top_type == 1) {
      for (F c=0; c<8; ++c)
        add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                             affine_scalar(c), targets);
      add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                           affine_scalar(0,1,0), targets);
    } else if (orbit.top_type == 2) {
      for (F c=0; c<8; ++c) {
        add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                             affine_scalar(c), targets);
        add_quadratic_branch(classifier, source, Branch::P_FIXED,
                             affine_scalar(c), targets);
      }
      add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                           affine_scalar(0,1,0), targets);
    } else if (orbit.top_type == 3) {
      for (F c=0; c<8; ++c) {
        Poly fixed = affine_scalar(c);
        add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                             fixed, targets);
        add_quadratic_branch(classifier, source, Branch::P_FIXED,
                             fixed, targets);
        add_quadratic_branch(classifier, source, Branch::SUM_FIXED,
                             fixed, targets);
      }
    } else if (orbit.top_type == 4) {
      for (F c=0; c<8; ++c)
        add_quadratic_branch(classifier, source, Branch::Q_FIXED,
                             affine_scalar(c), targets);
    }
    std::cerr << "exact orbit " << orbit.id << " edges=" << targets.size()
              << "\n";
  }
  return adjacency;
}

static void print_edges(const std::vector<std::set<int>>& adjacency) {
  for (int source=0; source<(int)adjacency.size(); ++source) {
    std::cout << "EDGES " << source;
    for (int target : adjacency[source]) std::cout << " " << target;
    std::cout << "\n";
  }
}

static void verify_generated_poset(const Classifier& classifier,
                                   const std::vector<std::set<int>>& adjacency) {
  const int n = adjacency.size();
  std::vector<std::set<int>> closure(n);
  for (int start=0; start<n; ++start) {
    closure[start].insert(start);
    std::vector<int> stack{start};
    while (!stack.empty()) {
      int source = stack.back();
      stack.pop_back();
      for (int target : adjacency[source]) {
        if (closure[start].insert(target).second) stack.push_back(target);
      }
    }
  }

  std::vector<int> component_of(n, -1);
  std::vector<std::vector<int>> components;
  for (int source=0; source<n; ++source) {
    if (component_of[source] >= 0) continue;
    int id = components.size();
    components.push_back({});
    for (int target=0; target<n; ++target) {
      if (component_of[target] < 0 && closure[source].count(target) &&
          closure[target].count(source)) {
        component_of[target] = id;
        components.back().push_back(target);
      }
    }
  }
  assert(components.size() == 110);

  std::vector<std::set<int>> component_reach(components.size());
  uint64_t total_size = 0;
  for (int source=0; source<n; ++source) {
    total_size += classifier.orbits[source].affine_size;
    for (int target : closure[source])
      component_reach[component_of[source]].insert(component_of[target]);
  }
  assert(total_size == (uint64_t(1) << 30));

  int cover_count = 0;
  for (int source=0; source<(int)components.size(); ++source) {
    for (int target : component_reach[source]) {
      if (source == target) continue;
      bool covered = false;
      for (int middle : component_reach[source]) {
        if (middle != source && middle != target &&
            component_reach[middle].count(target)) {
          covered = true;
          break;
        }
      }
      if (!covered) ++cover_count;
    }
  }
  assert(cover_count == 206);
  std::cerr << "generated poset: components=" << components.size()
            << " covers=" << cover_count << " maps=" << total_size << "\n";
}

static void print_orbits(const Classifier& classifier) {
  uint64_t total = 0;
  std::array<int,9> images{};
  for (const auto& orbit : classifier.orbits) {
    total += orbit.affine_size;
    images[orbit.image]++;
    std::cout << "ORBIT " << orbit.id << " top=" << orbit.top_type
              << " lower=" << orbit.representative_lower
              << " size=" << orbit.affine_size
              << " image=" << orbit.image << "\n";
  }
  std::cout << "TOTAL orbits=" << classifier.orbits.size()
            << " maps=" << total << " images";
  for (int i=1; i<=8; ++i) if (images[i]) std::cout << " " << i << ":" << images[i];
  std::cout << "\n";
}

static void random_transition_search(const Classifier& classifier,
                                     int trials_per_orbit) {
  std::mt19937 rng(20260816);
  std::uniform_int_distribution<int> coeff(0,7);
  std::vector<std::set<int>> adjacency(classifier.orbits.size());
  for (const auto& orbit : classifier.orbits) adjacency[orbit.id].insert(orbit.id);
  for (const auto& orbit : classifier.orbits) {
    Poly source = from_parts(orbit.normal_top, orbit.representative_lower);
    for (int trial=0; trial<trials_per_orbit; ++trial) {
      uint32_t pc=0, qc=0;
      for (int k=0; k<3; ++k) {
        pc |= uint32_t(coeff(rng)) << (3*k);
        qc |= uint32_t(coeff(rng)) << (3*k);
      }
      std::array<F,3> p2{}, q2{};
      for (F& value : p2) value = coeff(rng);
      for (F& value : q2) value = coeff(rng);
      if (orbit.top_type == 1 || orbit.top_type == 4) {
        q2 = {0,0,0};
      } else if (orbit.top_type == 2) {
        if (coeff(rng)&1) p2 = {0,0,0}; else q2 = {0,0,0};
      } else if (orbit.top_type == 3) {
        int branch = coeff(rng)%3;
        if (branch == 0) p2 = {0,0,0};
        else if (branch == 1) q2 = {0,0,0};
        else q2 = p2;
      } else if (orbit.top_type == 5) {
        p2 = q2 = {0,0,0};
      }
      for (int k=0; k<3; ++k) {
        pc |= uint32_t(p2[k]) << (3*(k+3));
        qc |= uint32_t(q2[k]) << (3*(k+3));
      }
      Poly target = compose(source, input_poly(pc), input_poly(qc));
      if (!degree_at_most_3(target)) continue;
      adjacency[orbit.id].insert(classifier.classify(target));
    }
    std::cerr << "sampled orbit " << orbit.id << " edges="
              << adjacency[orbit.id].size() << "\n";
  }
  for (int source=0; source<(int)adjacency.size(); ++source) {
    (void)source;
  }
  print_edges(adjacency);
}

int main(int argc, char** argv) {
  init_monomials();
  Classifier classifier;
  print_orbits(classifier);
  if (argc > 1 && std::string(argv[1]) == "exact") {
    auto adjacency = exact_transition_search(classifier);
    verify_generated_poset(classifier, adjacency);
    print_edges(adjacency);
  } else {
    int trials = argc > 1 ? std::stoi(argv[1]) : 0;
    if (trials) random_transition_search(classifier, trials);
  }
  return 0;
}
