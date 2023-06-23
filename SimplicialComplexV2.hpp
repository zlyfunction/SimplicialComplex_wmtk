#include <vector>
#include <set>
#include <cassert>
#include <queue>

struct Tuple;

struct Mesh
{
    Tuple sw(const Tuple &t, const int &d) const
    {
        throw std::exception("This is a dummy implementation!");
        return {};
    }
    bool is_boundary(const Tuple &t, const int &d) const
    {
        throw std::exception("This is a dummy implementation!");
        return false;
    }

    int cell_dimension() const
    {
        throw std::exception("This is a dummy implementation!");
        return 3;
    }
};

struct Tuple
{
    Tuple sw(const int &d, const Mesh *m) const
    {
        throw std::exception("This is a dummy implementation!");
        return m->sw(*this, d);
    }

    bool is_boundary(const Mesh *m) const
    {
        throw std::exception("This is a dummy implementation!");
        return false;
    }
};

// note that this code only works for triangle/tet meshes
// TODO: optimize it from n^2 to nlg(n)
class Simplex
{
    int _d;   // dimension
    Tuple _t; // tuple

public:
    Simplex(const int &d, const Tuple &t) : _d{d}, _t{t} {}

    int global_id() const { return -1; }
    int dimension() const { return _d; }
    const Tuple &tuple() const { return _t; }

    bool operator<(const Simplex &rhs) const
    {
        if (_d < rhs._d)
        {
            return true;
        }
        if (_d > rhs._d)
        {
            return false;
        }
        return global_id() < rhs.global_id();
    }

    bool operator==(const Simplex &rhs) const
    {
        return (_d == rhs._d) && (global_id() == rhs.global_id());
    }
};

class SimplicialComplex
{
private:
    std::set<Simplex> simplexes;
    const Mesh *m;

public:
    const std::set<Simplex> &get_simplices() const
    {
        return simplexes;
    }

    std::set<Simplex> get_simplices(const int &dim) const
    {
        std::set<Simplex> ret;
        for (const Simplex &s : simplexes)
        {
            if (s.dimension() == dim)
            {
                ret.insert(s);
            }
        }

        return ret;
    }

    // TODO why do we need that?
    int get_size() const
    {
        int ret = 0;
        for (int d = 0; d < 4; d++)
        {
            ret += simplexes[d].size();
        }
        return ret;
    }

    /**
     * @brief Add simplex to the complex if it is not already in it.
     *
     * @returns false if simplex is already in the complex
     */
    bool add_simplex(const Simplex &s)
    {
        assert(s.dimension() <= 4);
        const auto [it, was_successful] = simplexes.insert(s);
        return was_successful;
    }

    void unify_with_complex(const SimplicialComplex &other)
    {
        // this is N log(N) complexity
        for (const Simplex &s : other.get_simplices())
        {
            add_simplex(s);
        }
    }

    bool operator==(const SimplicialComplex &other) const
    {
        if (simplexes.size() != other.simplexes.size())
        {
            return false;
        }
        // this is N log(N) complexity
        for (const auto &t1 : simplexes)
        {
            const auto it = other.simplexes.find(t1);
            if (it == other.simplexes.end())
            {
                return false;
            }
        }

        return true;
    }

    SimplicialComplex &operator=(const SimplicialComplex &other)
    {
        if (this != &other)
        {
            this->simplexes = other.simplexes;
            this->m = other.m;
        }
        return *this;
    }

    SimplicialComplex(const Mesh *mm) : m{mm} {}

    const Mesh *get_mesh() const { return m; }
};

inline SimplicialComplex get_union(const SimplicialComplex &sc1, const SimplicialComplex &sc2)
{
    SimplicialComplex u = sc1;
    u.unify_with_complex(sc2);
    return u;
}

inline SimplicialComplex get_intersection(const SimplicialComplex &A, const SimplicialComplex &B)
{
    SimplicialComplex sc_union = A;
    SimplicialComplex sc_intersection(A.get_mesh());

    for (const auto &s : B.get_simplices())
    {
        if (!sc_union.add_simplex(s))
        {
            // s is already in A --> s is in the intersection of A and B
            sc_intersection.add_simplex(s);
        }
    }

    return sc_intersection;
}

//////////////////////////////////
// List of Operators
// bd: boundary
// clbd: closed boudnary
// st: start
// clst: closed star
// lnk: link
//////////////////////////////////

// ∂s
/**
 * @brief get the boundary of a simplex
 */
SimplicialComplex bd(const Simplex &s, const Mesh *m)
{
    SimplicialComplex SC(m);

    // exhaustive implementation
    switch (s.dimension())
    {
    case 3:                                                                        // bd(tet) = 4triangles + 6 edges + 4vertices
        SC.add_simplex(Simplex(0, s.tuple()));                                     // A
        SC.add_simplex(Simplex(0, s.tuple().sw(0, m)));                            // B
        SC.add_simplex(Simplex(0, s.tuple().sw(1, m).sw(0, m)));                   // C
        SC.add_simplex(Simplex(0, s.tuple().sw(2, m).sw(0, m)));                   // D
        SC.add_simplex(Simplex(1, s.tuple()));                                     // AB
        SC.add_simplex(Simplex(1, s.tuple().sw(1, m)));                            // AC
        SC.add_simplex(Simplex(1, s.tuple().sw(0, m).sw(1, m)));                   // BC
        SC.add_simplex(Simplex(1, s.tuple().sw(2, m).sw(1, m)));                   // AD
        SC.add_simplex(Simplex(1, s.tuple().sw(0, m).sw(2, m).sw(1, m)));          // BD
        SC.add_simplex(Simplex(1, s.tuple().sw(1, m).sw(0, m).sw(2, m).sw(1, m))); // CD
        SC.add_simplex(Simplex(2, s.tuple()));                                     // ABC
        SC.add_simplex(Simplex(2, s.tuple().sw(2, m)));                            // ABD
        SC.add_simplex(Simplex(2, s.tuple().sw(1, m).sw(2, m)));                   // ACD
        SC.add_simplex(Simplex(2, s.tuple().sw(0, m).sw(1, m).sw(2, m)));          // BCD
        break;
    case 2: // bd(triangle) = 3edges + 3vertices
        SC.add_simplex(Simplex(0, s.tuple()));
        SC.add_simplex(Simplex(0, s.tuple().sw(0, m)));
        SC.add_simplex(Simplex(0, s.tuple().sw(1, m).sw(0, m)));
        SC.add_simplex(Simplex(1, s.tuple()));
        SC.add_simplex(Simplex(1, s.tuple().sw(1, m)));
        SC.add_simplex(Simplex(1, s.tuple().sw(0, m).sw(1, m)));
        /* code */
        break;
    case 1:
        // bd(edge) = 2 vertices
        SC.add_simplex(Simplex(0, s.tuple()));
        SC.add_simplex(Simplex(0, s.tuple().sw(0, m)));
        /* code */
        break;
    case 0:
        break;
    default:
        assert(false);
        break;
    }

    return SC;
}

// ∂s∪{s}
/**
 * @brief get complex of a simplex and its boundary
 */
SimplicialComplex simplex_with_boundary(const Simplex &s, const Mesh *m)
{
    SimplicialComplex sc = bd(s, m);
    sc.add_simplex(s);
    return sc;
}

// Simplex s1,s2, check if A∩B!=∅
// check is intersect(∂s1, ∂s2) has intersections
/**
 * @brief check if simplices with their boundary intersect
 */
inline bool simplices_wbd_intersect(const Simplex &s1, const Simplex &s2, const Mesh *m)
{
    SimplicialComplex s1_bd = simplex_with_boundary(s1, m);
    SimplicialComplex s2_bd = simplex_with_boundary(s2, m);
    SimplicialComplex s1_s2_int = get_intersection(s1_bd, s2_bd);
    return (s1_s2_int.get_size() != 0);
}

SimplicialComplex clst(const Simplex &s, const Mesh *m)
{
    SimplicialComplex sc(m);
    const int &cell_dim = m->cell_dimension(); // TODO: 2 for trimesh, 3 for tetmesh need it in Mesh class

    if (cell_dim == 2)
    {
        switch (s.dimension())
        {
        case 0:
        {
            std::queue<Tuple> q;
            q.push(s.tuple());
            while (!q.empty())
            {
                const Tuple t = q.front();
                q.pop();
                if (sc.add_simplex(Simplex(2, t)))
                {
                    if (!t.is_boundary(m))
                    {
                        q.push(t.sw(2, m));
                    }
                    if (!t.sw(1, m).is_boundary(m))
                    {
                        q.push(t.sw(1, m).sw(2, m));
                    }
                }
            }
            break;
        }
        case 1:
            sc.add_simplex(Simplex(2, s.tuple()));
            if (!s.tuple().is_boundary(m))
            {
                sc.add_simplex(Simplex(2, s.tuple().sw(2, m)));
            }
            break;
        case 2:
            sc.add_simplex(s);
            break;
        default:
            assert(false);
            break;
        }
    }
    else if (cell_dim == 3)
    {
        switch (s.dimension())
        {
        case 0:
        {
            std::queue<Tuple> q;
            q.push(s.tuple());
            while (!q.empty())
            {
                Tuple t = q.front();
                q.pop();
                if (sc.add_simplex(Simplex(3, t)))
                {
                    const Tuple t1 = t;
                    const Tuple t2 = t.sw(2, m);
                    const Tuple t3 = t.sw(1, m).sw(2, m);
                    if (!t1.is_boundary(m))
                    {
                        q.push(t1.sw(3, m));
                    }
                    if (!t2.is_boundary(m))
                    {
                        q.push(t2.sw(3, m));
                    }
                    if (!t3.is_boundary(m))
                    {
                        q.push(t3.sw(3, m));
                    }
                }
            }
            break;
        }
        case 1:
        {
            std::queue<Tuple> q;
            q.push(s.tuple());
            while (!q.empty())
            {
                Tuple t = q.front();
                q.pop();
                if (sc.add_simplex(Simplex(3, t)))
                {
                    if (!t.is_boundary(m))
                    {
                        q.push(t.sw(3, m));
                    }
                    if (!t.sw(2, m).is_boundary(m))
                    {
                        q.push(t.sw(2, m).sw(3, m));
                    }
                }
            }
            break;
        }
        case 2:
        {
            sc.add_simplex(Simplex(3, s.tuple()));
            if (!s.tuple().is_boundary(m))
            {
                sc.add_simplex(Simplex(3, s.tuple().sw(3, m)));
            }
            break;
        }
        case 3:
        {
            sc.add_simplex(s);
            break;
        }
        default:
        {
            assert(false);
            break;
        }
        }
    }

    const auto top_simplices = sc.get_simplices();
    for (const Simplex &ts : top_simplices)
    {
        sc.unify_with_complex(bd(ts, m));
    }
    return sc;
}

SimplicialComplex lnk(const Simplex &s, const Mesh *m)
{
    SimplicialComplex sc_clst = clst(s, m);
    SimplicialComplex sc(m);
    for (const Simplex &ss : sc_clst.get_simplices())
    {
        if (!simplices_wbd_intersect(s, ss, m))
        {
            sc.add_simplex(ss);
        }
    }

    return sc;
}

SimplicialComplex st(const Simplex &s, const Mesh *m)
{
    SimplicialComplex sc_clst = clst(s, m);
    SimplicialComplex sc(m);
    sc.add_simplex(s);
    for (const Simplex &ss : sc_clst.get_simplices())
    {
        if (simplices_wbd_intersect(s, ss, m))
        {
            sc.add_simplex(ss);
        }
    }

    return sc;
}

//////////////////////////////////
// check link condition
// input Tuple t --> edge (a,b)
// check if lnk(a) ∩ lnk(b) == lnk(ab)
//////////////////////////////////
bool link_cond(Tuple t, const Mesh &m)
{
    // TODO: implement this
}

//////////////////////////////////
// k-ring
//////////////////////////////////
std::vector<Tuple> one_ring(Tuple t, const Mesh &m)
{
    Simplex s(0, t);
    SimplicialComplex s_st = st(s, m);
    std::vector<Tuple> Vs = s_st.get_simplices()[1];
    for (int i = 0; i < Vs.size(); i++)
    {
        if (is_same_simplex(t, Vs[i], 0, m))
        {
            Vs[i] = Vs[i].sw(0, m);
        }
    }
    return Vs;
}

std::vector<Tuple> k_ring(Tuple t, const Mesh &m, int k)
{
    assert(k >= 1);

    if (k == 1)
    {
        return one_ring(t, m);
    }
    else
    {
        std::vector<Tuple> t_one_ring = one_ring(t, m);
        std::vector < std::vector < Tuple >>> all_ts;
        for (auto tmp : t_one_ring)
        {
            all_ts.push_back(k_ring(tmp, m, k - 1));
        }

        std::vector<Tuple> ret;
        for (auto t_vec : all_ts)
        {
            for (auto tmp : t_vec)
            {
                bool flag = true;
                for (int i = 0; i < ret.size(); i++)
                {
                    if (is_same_simplex(ret[i], tmp, 0, m))
                    {
                        flag = false;
                        break;
                    }
                    if (flag)
                    {
                        ret.push_back(tmp);
                    }
                }
            }
        }

        return ret;
    }
}
