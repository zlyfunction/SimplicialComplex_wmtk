#include <vector>
#include <assert.h>

using Tuple = double; // random dummy to remove IDE errors

struct Mesh
{
    Tuple sw(const Tuple &t, const int &d) const
    {
        throw std::exception("This is a dummy implementation!");
        return 0;
    }
};

// note that this code only works for triangle/tet meshes
// TODO: optimize it from n^2 to nlg(n)
struct Simplex
{
    int d;
    Tuple t;
};

class SimplicialComplex
{
private:
    std::vector<Simplex> simplexes;
    const Mesh *m;

public:
    const std::vector<Simplex> &get_simplexes() const
    {
        return simplexes;
    }

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
        assert(s.d <= 4);
        for (const auto &tmp : simplexes)
        {
            if (tmp.d != s.d)
            {
                continue;
            }
            // TODO: need implement is_same_simplex for tuples
            // is_same_simplex(t1, t2, dim, m) returns if (dim, t1) and (dim, t2) are the same simplex in Mesh m
            // Need it in Mesh class
            if (is_same_simplex(tmp, s))
            {
                return false;
            }
        }
        simplexes.push_back(s);
        return true;
    }

    void unify_with_complex(const SimplicialComplex &other)
    {
        // TODO this is N^2 complexity!!!
        for (const Simplex &s : other.get_simplexes())
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
        // TODO this is N^2 complexity!!!
        for (const auto &t1 : simplexes)
        {
            bool flag = false;
            for (const auto &t2 : other.simplexes)
            {
                if (is_same_simplex(t1, t2))
                {
                    flag = true;
                    break;
                }
            }
            if (!flag)
                return false;
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

inline bool is_same_simplex(const Simplex s1, const Simplex &s2)
{
    throw std::exception("Implementation required!");
    return false;
}

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

    for (const auto &s : B.get_simplexes())
    {
        if (!sc_union.add_simplex(s))
        {
            // s is already in A --> s is in the intersection of A and B
            sc_intersection.add_simplex(s);
        }
    }

    return sc_intersection;
}

// Simplex s1,s2, check if A∩B!=∅
// check is intersect(∂s1, ∂s2) has intersections
inline bool is_intersect(Simplex s1, Simplex s2)
{
    SimplicialComplex s1_bd = clbd(s1);
    SimplicialComplex s2_bd = cldb(s2);
    SimplicialComplex s1_s2_int = get_intersection(s1_bd, s2_bd);
    return (s1_s2_int.get_size() != 0);
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
    SimplicialComplex sc(m);
    switch (s.d)
    {
    case 0: // vertex
        break;
    case 1: // edge
        sc.add_simplex({0, s.t});
        sc.add_simplex({0, m->sw(s.t, 0)});
        break;

        ////////////////////////////////////////
        //// TODO: GO ON HERE
    case 2: // face
        sc.add_simplex({1, s.t});
        sc.add_simplex({1, m->sw(s.t, 1)});
        sc.add_simplex({1, m->sw(m->sw(s.t, 0), 1)});
        break;

    case 3: // tet
        sc.add_simplex({2, s.t});
        sc.add_simplex({2, s.t.sw(2, m)});
        sc.add_simplex({2, s.t.sw(1, m).sw(2, m)});
        sc.add_simplex({2, s.t.sw(0, m).sw(1, m).sw(2, m)});

        break;

    default:
        assert(false);
        break;
    }

    return sc;
}

// ∂s∪{s}
SimplicialComplex clbd(const Simplex &s, const Mesh &m)
{
    SimplicialComplex SC = bd(s, m);
    SC.AddSimplex(s);
    return SC;
}

SimplicialComplex clst(const Simplex &s, const Mesh &m)
{
    SimplicialComplex SC(m);
    int dim = m.dim; // TODO: 2 for trimesh, 3 for tetmesh need it in Mesh class

    if (dim == 2)
    {
        switch (s.d)
        {
        case 0:
            std::queue<Tuple> q;
            q.push(s.t);
            while (!q.empty())
            {
                auto t = q.front();
                q.pop();
                if (SC.AddSimplex(Simplex(2, t)))
                {
                    if (!t.boundary(m))
                    {
                        q.push(t.sw(2, m));
                    }
                    if (!t.sw(1, m).boundary(m))
                    {
                        q.push(t.sw(1, m).sw(2, m));
                    }
                }
            }
            break;
        case 1:
            SC.add_simplex(Simplex(2, s.t));
            if (!s.t.boundary(m))
            {
                SC.add_simplex(Simplex(2, s.t.sw(2, m)));
            }
            break;
        case 2:
            SC.add_simplex(s);
            break;
        default:
            assert(false);
            break;
        }
    }
    else if (dim == 3)
    {
        switch (s.d)
        {
        case 0:
            std::queue<Tuple> q;
            q.push(s.t);
            while (!q.empty())
            {
                auto t = q.front();
                q.pop();
                if (SC.AddSimplex(Simplex(3, t)))
                {
                    if (!t.boundary(m))
                    {
                        q.push(t.sw(3, m));
                    }
                    if (!t.sw(2, m).boundary(m))
                    {
                        q.push(t.sw(2, m).sw(3, m));
                    }
                    if (!t.sw(1, m).sw(2, m).boundary(m))
                    {
                        q.push(t.sw(1, m).sw(2, m).sw(3, m));
                    }
                }
            }
            break;
        case 1:
            std::queue<Tuple> q;
            q.push(s.t);
            while (!q.empty())
            {
                auto t = q.front();
                q.pop();
                if (SC.AddSimplex(Simplex(3, t)))
                {
                    if (!t.boundary(m))
                    {
                        q.push(t.sw(3, m));
                    }
                    if (!t.sw(2, m).boundary(m))
                    {
                        q.push(t.sw(2, m).sw(3, m));
                    }
                }
            }
            break;
        case 2:
            SC.add_simplex(Simplex(3, s.t));
            if (!s.t.boundary(m))
            {
                SC.add_simplex(Simplex(3, s.t.sw(3, m)));
            }
            break;
        case 3:
            SC.add_simplex(s);
            break;
        default:
            assert(false);
            break;
        }
    }

    auto top_tuples = SC.get_simplexes[dim];
    for (Tuple t : top_tuples)
    {
        SC.union(bd(Simplex(dim, t), m));
    }
    return SC;
}

SimplicialComplex lnk(const Simplex &s, const Mesh &m)
{
    SimplicialComplex SC_clst = clst(s, m);
    SimplicialComplex SC(m);
    auto simplexes = SC_clst.get_simplexes();
    for (int d = 0; d < 4; d++)
    {
        for (auto t : simplexes[d])
        {
            if (!is_intersect(s, Simplex(d, t)))
            {
                SC.add_simplex(Simplex(d, t));
            }
        }
    }
    return SC;
}

SimplicialComplex st(const Simplex &s, const Mesh &m)
{
    SimplicialComplex SC_clst = clst(s, m);
    SimplicialComplex SC;
    SC.add_simplex(s);
    auto simplexes = SC_clst.get_simplexes();
    for (int d = s.d + 1; d < 4; d++)
    {
        for (auto t : simplexes[d])
        {
            if (is_intersect(s, Simplex(d, t)))
            {
                SC.add_simplex(Simplex(d, t));
            }
        }
    }
    return SC;
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
    std::vector<Tuple> Vs = s_st.get_simplexes()[1];
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